// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include <atomic>
#include <algorithm>

#include "StoredFragments.h"
#include "../simd/math.h"

namespace BPCells {

StoredFragmentsBase::StoredFragmentsBase(
    UIntReader &&cell,
    UIntReader &&start,
    UIntReader &&end,
    UIntReader &&end_max,
    ULongReader &&chr_ptr,
    std::unique_ptr<StringReader> &&chr_names,
    std::unique_ptr<StringReader> &&cell_names
)
    : cell(std::move(cell))
    , start(std::move(start))
    , end(std::move(end))
    , end_max(std::move(end_max))
    , chr_ptr(std::move(chr_ptr))
    , chr_names(std::move(chr_names))
    , cell_names(std::move(cell_names)) {}

// Read end_max_buf from end_max iterator, making it equal to the values between
// start_idx and end_idx
void StoredFragmentsBase::readEndMaxBuf(uint64_t start_idx, uint64_t end_idx) {
    if (start_idx == end_idx) {
        end_max_buf.resize(0);
        return;
    }
    // Read the end_max buffer
    end_max_buf.resize((end_idx - 1) / 128 - start_idx / 128 + 1);
    end_max.seek(start_idx / 128);
    uint64_t i = 0;
    while (true) {
        end_max.ensureCapacity(1);
        uint64_t load_amount = std::min(end_max.capacity(), (uint64_t) end_max_buf.size() - i);

        std::memmove(&end_max_buf[i], end_max.data(), load_amount * sizeof(uint32_t));
        i += load_amount;
        if (i >= end_max_buf.size()) break;
        end_max.advance(load_amount);
    }
    // Handle starting a chromosome on a non-multiple of 128, such end_max[0]
    // may be the max of the prior chromosome rather than max of the new one
    if (start_idx % 128 != 0 && end_max_buf.size() > 1) {
        end_max_buf[0] = std::min(end_max_buf[0], end_max_buf[1]);
    }
}

bool StoredFragmentsBase::isSeekable() const { return true; }
void StoredFragmentsBase::seek(uint32_t chr_id, uint32_t base) {
    if (chr_id != current_chr) {
        if ((int64_t)chr_id >= chrCount()) {
            // Seeking to a chromosome larger than exists in the fragments.
            // Make it so next load and nextChr calls will return false
            current_chr = chr_id;
            current_idx = UINT64_MAX;
            chr_start_ptr = 0;
            chr_end_ptr = 0;
            return;
        }
        current_chr = chr_id;
        chr_ptr.seek(chr_id * 2);
        chr_start_ptr = chr_ptr.read_one();
        chr_end_ptr = chr_ptr.read_one();

        readEndMaxBuf(chr_start_ptr, chr_end_ptr);
    }

    // Binary search for base in end_max
    uint64_t current_block = chr_start_ptr / 128 +
                             std::upper_bound(end_max_buf.begin(), end_max_buf.end(), base) -
                             end_max_buf.begin();

    // Add this max in case we're seeking to the first block in a chromosome that
    // doesn't start on a multiple of 128
    current_idx = std::max(chr_start_ptr, current_block * 128);

    start.seek(current_idx);
    end.seek(current_idx);
    cell.seek(current_idx);
    current_capacity = 0;
}

void StoredFragmentsBase::restart() {
    current_chr = UINT32_MAX;
    current_idx = UINT64_MAX;
    chr_ptr.seek(0);
}

int StoredFragmentsBase::chrCount() const { return chr_names->size(); }
int StoredFragmentsBase::cellCount() const { return cell_names->size(); }

const char *StoredFragmentsBase::chrNames(uint32_t chr_id) { return chr_names->get(chr_id); }
const char *StoredFragmentsBase::cellNames(uint32_t cell_id) { return cell_names->get(cell_id); }

bool StoredFragmentsBase::nextChr() {
    current_chr += 1;
    if ((int64_t)current_chr >= chrCount()) {
        current_chr -= 1;
        current_idx = UINT64_MAX;
        return false;
    }

    chr_start_ptr = chr_ptr.read_one();
    chr_end_ptr = chr_ptr.read_one();
    // Check if we need to perform any seeks, or if we're all good
    // since we just finished a chromosome
    if (current_idx != chr_start_ptr) {
        // We need to perform seeks to get to the right data
        // reading location
        start.seek(chr_start_ptr);
        end.seek(chr_start_ptr);
        cell.seek(chr_start_ptr);
        current_capacity = 0;
    }
    current_idx = chr_start_ptr;
    readEndMaxBuf(chr_start_ptr, chr_end_ptr);
    return true;
}
uint32_t StoredFragmentsBase::currentChr() const { return current_chr; }

uint32_t StoredFragmentsBase::capacity() const { return current_capacity; }

uint32_t *StoredFragmentsBase::cellData() { return cell.data(); }
uint32_t *StoredFragmentsBase::startData() { return start.data(); }
uint32_t *StoredFragmentsBase::endData() { return end.data(); }

StoredFragments StoredFragments::openUnpacked(
    ReaderBuilder &rb,
    std::unique_ptr<StringReader> &&chr_names,
    std::unique_ptr<StringReader> &&cell_names
) {

    ULongReader chr_ptr;
    if (rb.readVersion() == "unpacked-fragments-v1") {
        chr_ptr = rb.openUIntReader("chr_ptr").convert<uint64_t>();
    } else if (rb.readVersion() == "unpacked-fragments-v2") {
        chr_ptr = rb.openULongReader("chr_ptr");
    } else {
        throw std::runtime_error(
            std::string("Version does not match unpacked-fragments-v2: ") + rb.readVersion()
        );
    }

    if (!chr_names) chr_names = rb.openStringReader("chr_names");
    if (!cell_names) cell_names = rb.openStringReader("cell_names");

    return StoredFragments(
        rb.openUIntReader("cell"),
        rb.openUIntReader("start"),
        rb.openUIntReader("end"),
        rb.openUIntReader("end_max"),
        std::move(chr_ptr),
        std::move(chr_names),
        std::move(cell_names)
    );
}

bool StoredFragments::load() {
    if (current_idx >= chr_end_ptr) {
        return false;
    }

    cell.advance(current_capacity);
    start.advance(current_capacity);
    end.advance(current_capacity);

    // Load cell, start, or end if necessary
    if (cell.capacity() == 0) cell.ensureCapacity(1);
    if (start.capacity() == 0) start.ensureCapacity(1);
    if (end.capacity() == 0) end.ensureCapacity(1);

    current_capacity =
        std::min({cell.capacity(), start.capacity(), end.capacity(), chr_end_ptr - current_idx});
    current_idx += current_capacity;
    return true;
}

StoredFragmentsPacked StoredFragmentsPacked::openPacked(
    ReaderBuilder &rb,
    uint32_t load_size,
    std::unique_ptr<StringReader> &&chr_names,
    std::unique_ptr<StringReader> &&cell_names
) {
    ULongReader chr_ptr, start_idx_offsets, end_idx_offsets, cell_idx_offsets;
    if (rb.readVersion() == "packed-fragments-v1") {
        chr_ptr = rb.openUIntReader("chr_ptr").convert<uint64_t>();
        start_idx_offsets = ConstNumReader<uint64_t>::create({0, UINT64_MAX});
        end_idx_offsets = ConstNumReader<uint64_t>::create({0, UINT64_MAX});
        cell_idx_offsets = ConstNumReader<uint64_t>::create({0, UINT64_MAX});
    } else if (rb.readVersion() == "packed-fragments-v2") {
        chr_ptr = rb.openULongReader("chr_ptr");
        start_idx_offsets = rb.openULongReader("start_idx_offsets");
        end_idx_offsets = rb.openULongReader("end_idx_offsets");
        cell_idx_offsets = rb.openULongReader("cell_idx_offsets");
    } else {
        throw std::runtime_error(
            std::string("Version does not match packed-fragments-v2: ") + rb.readVersion()
        );
    }

    uint64_t count = 0;
    for (uint64_t i = 0; i < chr_ptr.size(); i++) {
        count = std::max(count, chr_ptr.read_one());
    }
    chr_ptr.seek(0);

    if (!chr_names) chr_names = rb.openStringReader("chr_names");
    if (!cell_names) cell_names = rb.openStringReader("cell_names");

    return StoredFragmentsPacked(
        UIntReader(
            std::make_unique<BP128UIntReader>(
                rb.openUIntReader("cell_data"),
                rb.openUIntReader("cell_idx"),
                std::move(cell_idx_offsets),
                count
            ),
            load_size,
            load_size
        ),
        UIntReader(
            std::make_unique<BP128_D1_UIntReader>(
                rb.openUIntReader("start_data"),
                rb.openUIntReader("start_idx"),
                std::move(start_idx_offsets),
                rb.openUIntReader("start_starts"),
                count
            ),
            load_size,
            load_size
        ),
        UIntReader(
            std::make_unique<BP128UIntReader>(
                rb.openUIntReader("end_data"),
                rb.openUIntReader("end_idx"),
                std::move(end_idx_offsets),
                count
            ),
            load_size,
            load_size
        ),
        rb.openUIntReader("end_max"),
        std::move(chr_ptr),
        std::move(chr_names),
        std::move(cell_names)
    );
}

bool StoredFragmentsPacked::load() {
    if (current_idx >= chr_end_ptr) {
        return false;
    }

    cell.advance(current_capacity);
    start.advance(current_capacity);
    end.advance(current_capacity);

    // Load cell, start, or end if necessary
    if (cell.capacity() == 0) cell.ensureCapacity(1);
    if (start.capacity() == 0) start.ensureCapacity(1);
    if (end.capacity() == 0) end.ensureCapacity(1);

    current_capacity =
        std::min({cell.capacity(), start.capacity(), end.capacity(), chr_end_ptr - current_idx});
    current_idx += current_capacity;

    uint32_t *start = startData();
    uint32_t *end = endData();
    simd::add(end, start, current_capacity);
    return true;
}

StoredFragmentsWriter::StoredFragmentsWriter(
    UIntWriter &&cell,
    UIntWriter &&start,
    UIntWriter &&end,
    UIntWriter &&end_max,
    ULongWriter &&chr_ptr,
    std::unique_ptr<StringWriter> &&chr_names,
    std::unique_ptr<StringWriter> &&cell_names,
    bool subtract_start_from_end
)
    : cell(std::move(cell))
    , start(std::move(start))
    , end(std::move(end))
    , end_max(std::move(end_max))
    , chr_ptr(std::move(chr_ptr))
    , chr_names(std::move(chr_names))
    , cell_names(std::move(cell_names))
    , subtract_start_from_end(subtract_start_from_end) {}

StoredFragmentsWriter StoredFragmentsWriter::createUnpacked(WriterBuilder &wb) {
    wb.writeVersion("unpacked-fragments-v2");
    return StoredFragmentsWriter(
        wb.createUIntWriter("cell"),
        wb.createUIntWriter("start"),
        wb.createUIntWriter("end"),
        wb.createUIntWriter("end_max"),
        wb.createULongWriter("chr_ptr"),
        wb.createStringWriter("chr_names"),
        wb.createStringWriter("cell_names"),
        false
    );
}

StoredFragmentsWriter StoredFragmentsWriter::createPacked(WriterBuilder &wb, uint32_t buffer_size) {
    wb.writeVersion("packed-fragments-v2");

    return StoredFragmentsWriter(
        UIntWriter(
            std::make_unique<BP128UIntWriter>(
                wb.createUIntWriter("cell_data"),
                wb.createUIntWriter("cell_idx"),
                wb.createULongWriter("cell_idx_offsets")
            ),
            buffer_size
        ),
        UIntWriter(
            std::make_unique<BP128_D1_UIntWriter>(
                wb.createUIntWriter("start_data"),
                wb.createUIntWriter("start_idx"),
                wb.createULongWriter("start_idx_offsets"),
                wb.createUIntWriter("start_starts")
            ),
            buffer_size
        ),
        UIntWriter(
            std::make_unique<BP128UIntWriter>(
                wb.createUIntWriter("end_data"),
                wb.createUIntWriter("end_idx"),
                wb.createULongWriter("end_idx_offsets")
            ),
            buffer_size
        ),
        wb.createUIntWriter("end_max"),
        wb.createULongWriter("chr_ptr"),
        wb.createStringWriter("chr_names"),
        wb.createStringWriter("cell_names"),
        true
    );
}

void StoredFragmentsWriter::write(FragmentLoader &fragments, std::atomic<bool> *user_interrupt) {
    uint32_t cur_end_max = 0;
    uint32_t prev_end_max = 0;
    uint64_t idx = 0;

    std::vector<uint64_t> chr_ptr_buf;

    uint64_t write_capacity =
        std::min({cell.maxCapacity(), start.maxCapacity(), end.maxCapacity()});

    fragments.restart();
    while (fragments.nextChr()) {
        uint32_t chr_id = fragments.currentChr();

        if (chr_id * 2 + 2 >= chr_ptr_buf.size()) {
            chr_ptr_buf.resize(chr_id * 2 + 2);
        }
        chr_ptr_buf[chr_id * 2] = idx;

        // This helps handle the case of multiple chromosomes ending within
        // the same block of 128 fragments
        prev_end_max = std::max(prev_end_max, cur_end_max);
        cur_end_max = 0;
        while (fragments.load()) {
            uint64_t capacity = fragments.capacity();

            const uint32_t *in_cell_data = fragments.cellData();
            const uint32_t *in_start_data = fragments.startData();
            uint32_t *in_end_data = fragments.endData();

            uint64_t i = 0;

            // Finish calculating end_max for first partial chunk
            if (idx % 128 != 0) {
                const uint64_t batch = std::min(capacity, 128 - idx % 128);
                cur_end_max = std::max(cur_end_max, simd::max(in_end_data, batch));
                i += batch;

                if ((idx + i) % 128 == 0) {
                    end_max.write_one(std::max(cur_end_max, prev_end_max));
                    prev_end_max = 0;
                }
            }
            // Calculate end_max for middle chunks
            for (; i + 128 <= capacity; i += 128) {
                cur_end_max = std::max(cur_end_max, simd::max(in_end_data + i, 128));
                end_max.write_one(std::max(cur_end_max, prev_end_max));
                prev_end_max = 0;
            }
            // Calculate end_max for trailing partial chunk
            if (i < capacity) {
                cur_end_max = std::max(cur_end_max, simd::max(in_end_data + i, capacity - i));
            }

            if (subtract_start_from_end) {
                simd::sub(in_end_data, in_start_data, capacity);
            }

            // Write data in chunks in case our output read capacity is less than our
            // input read amount
            for (uint64_t written = 0; written < capacity;) {
                uint64_t write_amount = std::min(write_capacity, capacity - written);
                cell.ensureCapacity(write_amount);
                start.ensureCapacity(write_amount);
                end.ensureCapacity(write_amount);

                std::memmove(cell.data(), in_cell_data + written, write_amount * sizeof(uint32_t));
                std::memmove(
                    start.data(), in_start_data + written, write_amount * sizeof(uint32_t)
                );
                std::memmove(end.data(), in_end_data + written, write_amount * sizeof(uint32_t));

                cell.advance(write_amount);
                start.advance(write_amount);
                end.advance(write_amount);
                written += write_amount;
            }

            idx += capacity;

            if (user_interrupt != NULL && *user_interrupt) return;
        }
        chr_ptr_buf[chr_id * 2 + 1] = idx;
    }
    if (idx % 128 != 0) {
        end_max.write_one(std::max(cur_end_max, prev_end_max));
    }

    if (fragments.chrCount() * 2 > (int64_t)chr_ptr_buf.size()) {
        chr_ptr_buf.resize(fragments.chrCount() * 2, idx);
    }

    for (auto &ptr : chr_ptr_buf) {
        chr_ptr.write_one(ptr);
    }

    cell.finalize();
    start.finalize();
    end.finalize();
    end_max.finalize();
    chr_ptr.finalize();

    // Get cell and chromosome names. This probably incurs a few extra copies,
    // but it shouldn't matter since writing the actual fragments should dominate cost
    std::vector<std::string> cell_names;
    for (int i = 0;; i++) {
        const char *cell_name = fragments.cellNames(i);
        if (cell_name == NULL) break;
        cell_names.push_back(std::string(cell_name));
    }
    this->cell_names->write(VecStringReader(cell_names));

    std::vector<std::string> chr_names;
    for (int i = 0;; i++) {
        const char *chr_name = fragments.chrNames(i);
        if (chr_name == NULL) break;
        chr_names.push_back(std::string(chr_name));
    }
    this->chr_names->write(VecStringReader(chr_names));
}

} // end namespace BPCells
