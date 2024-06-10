// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "bp128.h"

#include "../simd/bp128.h"

namespace BPCells {

// ################### BP128 (Vanilla + Base class) #########################################
BP128UIntReader::BP128UIntReader(
    UIntReader &&data, UIntReader &&idx, ULongReader &&idx_offsets, uint64_t count
)
    : data(std::move(data))
    , idx(std::move(idx))
    , idx_offsets(std::move(idx_offsets))
    , count(count)
    , prev_idx(this->idx.read_one())
    , prev_offset_boundary(this->idx_offsets.read_one())
    , next_offset_boundary(this->idx_offsets.read_one()) {}

uint64_t BP128UIntReader::size() const { return count; }

void BP128UIntReader::seek(uint64_t new_pos) {
    if (pos % 128 != 0 && new_pos % 128 != 0 && new_pos / 128 == pos / 128) {
        pos = new_pos;
        return;
    }
    pos = new_pos;

    // Reset the offset boundaries if we're seeking too far backwards
    if (pos / 128 < prev_offset_boundary) {
        idx_offsets.seek(0);
        idx_offset = 0;
        prev_offset_boundary = idx_offsets.read_one();
        next_offset_boundary = idx_offsets.read_one();
    }
    // Handle if we have passed an offset boundary, and adjust idx_offset
    while (pos / 128 >= next_offset_boundary) {
        prev_offset_boundary = next_offset_boundary;
        next_offset_boundary = idx_offsets.read_one();
        idx_offset += OFFSET_INCREMENT;
    }

    seekLoaders();

    if (pos % 128 != 0) load128(buf);
}

uint64_t BP128UIntReader::load(uint32_t *out, uint64_t count) {
    count = std::min(this->count - pos, count);
    uint64_t i = 0;
    if (pos % 128 != 0) {
        i = std::min(count, 128 - pos % 128);
        std::memmove(out, buf + pos % 128, i * sizeof(uint32_t));
        pos += i;
    }
    for (; i + 128 <= count; i += 128) {
        load128(out + i);
        pos += 128;
    }
    if (i == 0 && count > 0) {
        // Handle leftovers. If we get here we know count < 128 and
        // we started aligned to a boundary of 128
        load128(buf);
        std::memmove(out, buf, count * sizeof(uint32_t));
        pos += count;
        i += count;
    }
    return i;
}

void BP128UIntReader::setOffsetIncrement(uint64_t val) {
    OFFSET_INCREMENT = val;
}

void BP128UIntReader::load128(uint32_t *out) {
    if (1 + pos / 128 >= next_offset_boundary) {
        // Check if the next_idx needs a new offset increment
        prev_offset_boundary = next_offset_boundary;
        next_offset_boundary = idx_offsets.read_one();
        idx_offset += OFFSET_INCREMENT;
    }
    uint64_t next_idx = idx.read_one() + idx_offset;
    uint32_t bits = (next_idx - prev_idx) / 4;
    data.ensureCapacity(bits * 4);
    load128(data.data(), out, bits);
    data.advance(bits * 4);
    prev_idx = next_idx;
}

void BP128UIntReader::load128(uint32_t *in, uint32_t *out, uint32_t bits) {
    simd::bp128::unpack(in, out, bits);
}

void BP128UIntReader::seekLoaders() {
    idx.seek(pos / 128);
    prev_idx = idx.read_one() + idx_offset;
    data.seek(prev_idx);
}

BP128UIntWriter::BP128UIntWriter(UIntWriter &&data, UIntWriter &&idx, ULongWriter &&idx_offsets)
    : data(std::move(data))
    , idx(std::move(idx))
    , idx_offsets(std::move(idx_offsets)) {

    this->idx.write_one(0);
    this->idx_offsets.write_one(0);
}

uint64_t BP128UIntWriter::write(uint32_t *in, uint64_t count) {
    uint64_t i = 0;
    if (buf_pos != 0 || count < 128) {
        i = std::min(count, 128 - buf_pos);
        std::memmove(buf + buf_pos, in, i * sizeof(uint32_t));
        buf_pos = (buf_pos + i) % 128;
        if (buf_pos == 0) pack128(buf);
    }
    for (; i + 128 <= count; i += 128) {
        pack128(in + i);
    }
    return i;
}

void BP128UIntWriter::finalize() {
    if (buf_pos != 0) {
        for (; buf_pos < 128; buf_pos++) {
            buf[buf_pos] = buf[buf_pos - 1];
        }
        pack128(buf);
    }
    idx_offsets.write_one(1 + pos/128);
    data.finalize();
    idx.finalize();
    idx_offsets.finalize();
}

void BP128UIntWriter::setOffsetIncrement(uint64_t val) {
    OFFSET_INCREMENT = val;
}

void BP128UIntWriter::pack128(uint32_t *in) {
    pos += 128;
    uint32_t bits = this->bits(in);
    data.ensureCapacity(bits * 4);
    pack128(in, data.data(), bits);
    data.advance(bits * 4);
    cur_idx += bits * 4;
    if (cur_idx >= OFFSET_INCREMENT) {
        cur_idx -= OFFSET_INCREMENT;
        idx_offsets.write_one(pos / 128);
    }
    idx.write_one(cur_idx);
}

void BP128UIntWriter::pack128(uint32_t *in, uint32_t *out, uint32_t bits) {
    simd::bp128::pack(in, out, bits);
}

uint32_t BP128UIntWriter::bits(const uint32_t *in) const { return simd::bp128::maxbits(in); }

// ######################## BP128 (D1) #########################################

BP128_D1_UIntReader::BP128_D1_UIntReader(
    UIntReader &&data,
    UIntReader &&idx,
    ULongReader &&idx_offsets,
    UIntReader &&starts,
    uint64_t count
)
    : BP128UIntReader(std::move(data), std::move(idx), std::move(idx_offsets), count)
    , starts(std::move(starts)) {}

void BP128_D1_UIntReader::seekLoaders() { starts.seek(pos / 128); BP128UIntReader::seekLoaders(); }

void BP128_D1_UIntReader::load128(uint32_t *in, uint32_t *out, uint32_t bits) {
    uint32_t start = starts.read_one();
    simd::bp128::unpack_d1(start, in, out, bits);
}

BP128_D1_UIntWriter::BP128_D1_UIntWriter(
    UIntWriter &&data, UIntWriter &&idx, ULongWriter &&idx_offsets, UIntWriter &&starts
)
    : BP128UIntWriter(std::move(data), std::move(idx), std::move(idx_offsets))
    , starts(std::move(starts)) {}

void BP128_D1_UIntWriter::pack128(uint32_t *in, uint32_t *out, uint32_t bits) {
    uint32_t start = in[0];
    simd::bp128::pack_d1(start, in, out, bits);
    starts.write_one(start);
}

uint32_t BP128_D1_UIntWriter::bits(const uint32_t *in) const {
    return simd::bp128::maxbits_d1(in[0], in);
}

void BP128_D1_UIntWriter::finalize() {
    BP128UIntWriter::finalize();
    starts.finalize();
}

// ######################## BP128 (D1Z) #########################################

BP128_D1Z_UIntReader::BP128_D1Z_UIntReader(
    UIntReader &&data, UIntReader &&idx, ULongReader &&idx_offsets, UIntReader &&starts, uint64_t count
)
    : BP128UIntReader(std::move(data), std::move(idx), std::move(idx_offsets), count)
    , starts(std::move(starts)) {}

void BP128_D1Z_UIntReader::seekLoaders() { starts.seek(pos / 128); BP128UIntReader::seekLoaders();}

void BP128_D1Z_UIntReader::load128(uint32_t *in, uint32_t *out, uint32_t bits) {
    uint32_t start = starts.read_one();
    simd::bp128::unpack_d1z(start, in, out, bits);
}

BP128_D1Z_UIntWriter::BP128_D1Z_UIntWriter(
    UIntWriter &&data, UIntWriter &&idx, ULongWriter &&idx_offsets, UIntWriter &&starts
)
    : BP128UIntWriter(std::move(data), std::move(idx), std::move(idx_offsets))
    , starts(std::move(starts)) {}

void BP128_D1Z_UIntWriter::pack128(uint32_t *in, uint32_t *out, uint32_t bits) {
    uint32_t start = in[0];
    simd::bp128::pack_d1z(start, in, out, bits);
    starts.write_one(start);
}

uint32_t BP128_D1Z_UIntWriter::bits(const uint32_t *in) const {
    return simd::bp128::maxbits_d1z(in[0], in);
}

void BP128_D1Z_UIntWriter::finalize() {
    BP128UIntWriter::finalize();
    starts.finalize();
}

// ######################## BP128 (FOR) #########################################
//  This is just FOR encoding using a constant 1 as the frame of reference.
//  (i.e. to encode values of an integer sparse matrix)

BP128_FOR_UIntReader::BP128_FOR_UIntReader(UIntReader &&data, UIntReader &&idx, ULongReader &&idx_offsets, uint64_t count)
    : BP128UIntReader(std::move(data), std::move(idx), std::move(idx_offsets), count) {}



void BP128_FOR_UIntReader::load128(uint32_t *in, uint32_t *out, uint32_t bits) {
    simd::bp128::unpack_FOR(1, in, out, bits);
}

BP128_FOR_UIntWriter::BP128_FOR_UIntWriter(UIntWriter &&data, UIntWriter &&idx, ULongWriter &&idx_offsets)
    : BP128UIntWriter(std::move(data), std::move(idx), std::move(idx_offsets)) {}

void BP128_FOR_UIntWriter::pack128(uint32_t *in, uint32_t *out, uint32_t bits) {
    simd::bp128::pack_FOR(1, in, out, bits);
}

uint32_t BP128_FOR_UIntWriter::bits(const uint32_t *in) const {
    return simd::bp128::maxbits_FOR(1, in);
}

} // end namespace BPCells
