// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "MergeFragments.h"
#include <algorithm>
#include <cstring>
#include <functional>
#include <unordered_map>

#include <dary_heap/dary_heap.hpp>
#include "../utils/radix_sort.h"

namespace BPCells {

MergeFragments::ChunkedLoader::ChunkedLoader(
    std::unique_ptr<FragmentLoader> &&loader, uint32_t cell_offset, uint32_t chunk_size
)
    : frags(std::move(loader))
    , cell_offset(cell_offset)
    , chunk_size(chunk_size) {}

uint32_t
MergeFragments::ChunkedLoader::load_chunk(uint32_t *start, uint32_t *end, uint32_t *cell) {
    uint32_t copied = 0;
    // Use up any remaining loaded data
    if (loaded > 0) {
        uint32_t copy_amount = std::min<uint32_t>(chunk_size, loaded);
        uint32_t offset = frags->capacity() - loaded;
        std::memmove(start, frags->startData() + offset, sizeof(uint32_t) * copy_amount);
        std::memmove(end, frags->endData() + offset, sizeof(uint32_t) * copy_amount);
        std::memmove(cell, frags->cellData() + offset, sizeof(uint32_t) * copy_amount);
        copied += copy_amount;
        loaded -= copy_amount;
    }

    // Call load repeatedly until we fill the chunk or run out of data
    while (copied < chunk_size) {
        if (!frags->load()) break;
        loaded = frags->capacity();
        uint32_t copy_amount = std::min<uint32_t>(chunk_size - copied, loaded);
        std::memmove(start + copied, frags->startData(), sizeof(uint32_t) * copy_amount);
        std::memmove(end + copied, frags->endData(), sizeof(uint32_t) * copy_amount);
        std::memmove(cell + copied, frags->cellData(), sizeof(uint32_t) * copy_amount);
        copied += copy_amount;
        loaded -= copy_amount;
    }
    for (uint32_t i = 0; i < copied; i++) {
        cell[i] += cell_offset;
    }
    return copied;
}

uint32_t MergeFragments::ChunkedLoader::peek_start() {
    if (loaded == 0) {
        if (!frags->load()) return UINT32_MAX;
        loaded = frags->capacity();
    }
    uint32_t offset = frags->capacity() - loaded;
    return frags->startData()[offset];
}

void MergeFragments::ChunkedLoader::seek(uint32_t chr_id, uint32_t base) {
    loaded = 0;
    frags->seek(chr_id, base);
}

void MergeFragments::ChunkedLoader::restart() {
    loaded = 0;
    frags->restart();
}

int MergeFragments::ChunkedLoader::chrCount() const { return frags->chrCount(); }
int MergeFragments::ChunkedLoader::cellCount() const { return frags->cellCount(); }

const char *MergeFragments::ChunkedLoader::chrNames(uint32_t chr_id) {
    return frags->chrNames(chr_id);
}
const char *MergeFragments::ChunkedLoader::cellNames(uint32_t cell_id) {
    return frags->cellNames(cell_id);
}

bool MergeFragments::ChunkedLoader::nextChr() {
    loaded = 0;
    return frags->nextChr();
}

uint32_t MergeFragments::ChunkedLoader::currentChr() const { return frags->currentChr(); }

MergeFragments::MergeFragments(
    std::vector<std::unique_ptr<FragmentLoader>> &&fragments,
    const std::vector<std::string> &chr_order,
    uint32_t load_size,
    uint32_t chunk_size
)
    : load_size(load_size)
    , chunk_size(chunk_size)
    // load_size <= chunk_size * (chunks_per_input - 1) * fragments.size()
    , chunks_per_input(
          1 + (load_size + chunk_size * fragments.size() - 1) / (chunk_size * fragments.size())
      )
    , chr_order(chr_order)
    , start((chunks_per_input + 1) * chunk_size * fragments.size())
    , end((chunks_per_input + 1) * chunk_size * fragments.size())
    , cell((chunks_per_input + 1) * chunk_size * fragments.size())
    , start_buf((chunks_per_input + 1) * chunk_size * fragments.size())
    , end_buf((chunks_per_input + 1) * chunk_size * fragments.size())
    , cell_buf((chunks_per_input + 1) * chunk_size * fragments.size())
    , heap(fragments.size()) {

    if (fragments.size() < 2) throw std::runtime_error("Must have >= 2 fragments to merge");

    // Find the cell_id offsets to apply to each fragment
    cell_id_offset.push_back(0);
    for (uint32_t i = 0; i < fragments.size(); i++) {
        if (fragments[i]->cellCount() == -1 || fragments[i]->chrCount() == -1 ||
            !fragments[i]->isSeekable())
            throw std::runtime_error(
                "MergeFragments Error: all input fragments to merge must have known cell + chr "
                "counts and be seekable. Convert inputs to BPCells format first if needed."
            );
        cell_id_offset.push_back(cell_id_offset.back() + fragments[i]->cellCount());
    }

    // Find the chromosome ordering for each input fragment
    std::unordered_map<std::string, uint32_t> chr_order_lookup;
    for (uint32_t j = 0; j < chr_order.size(); j++) {
        chr_order_lookup[chr_order[j]] = j;
    }
    source_chr.resize(fragments.size());
    for (uint32_t i = 0; i < fragments.size(); i++) {
        auto &f = fragments[i];
        source_chr[i].resize(chr_order.size());
        for (auto &x : source_chr[i]) {
            x = UINT32_MAX;
        }
        int32_t chr_count = f->chrCount();
        for (int32_t j = 0; j < chr_count; j++) {
            if (chr_order_lookup.find(f->chrNames(j)) != chr_order_lookup.end()) {
                source_chr[i][chr_order_lookup[f->chrNames(j)]] = j;
            } else {
                throw std::runtime_error(
                    "MergeFragments Error: Input index " + std::to_string(i) + " has chromosome " +
                    std::string(f->chrNames(j)) + " which is not included in the output ordering."
                );
            }
        }
    }

    // Wrap input fragments in ChunkedLoader
    for (uint32_t i = 0; i < fragments.size(); i++) {
        frags.push_back(ChunkedLoader(std::move(fragments[i]), cell_id_offset[i], chunk_size));
    }
}

bool MergeFragments::isSeekable() const { return true; }

void MergeFragments::seek(uint32_t chr_id, uint32_t base) {
    for (uint32_t i = 0; i < frags.size(); i++) {
        if (chr_id < chr_order.size()) {
            frags[i].seek(source_chr[i][chr_id], base);
        } else {
            frags[i].seek(UINT32_MAX, base);
        }
    }

    loaded = 0;
    available = 0;
    current_chr = chr_id;
}

void MergeFragments::restart() {
    for (auto &&f : frags) {
        f.restart();
    }
    loaded = 0;
    available = 0;
    current_chr = UINT32_MAX;
}

int MergeFragments::chrCount() const { return chr_order.size(); }

int MergeFragments::cellCount() const { return cell_id_offset.back(); }

const char *MergeFragments::chrNames(uint32_t chr_id) {
    if (chr_id >= chr_order.size()) return NULL;
    return chr_order[chr_id].c_str();
}

const char *MergeFragments::cellNames(uint32_t cell_id) {
    auto it = std::upper_bound(cell_id_offset.begin(), cell_id_offset.end(), cell_id);
    uint32_t idx = it - cell_id_offset.begin() - 1;

    if (idx == frags.size()) idx--;

    return frags[idx].cellNames(cell_id - cell_id_offset[idx]);
}

bool MergeFragments::nextChr() {
    loaded = 0;
    available = 0;
    current_chr += 1;
    if ((int64_t)current_chr >= chrCount()) {
        current_chr -= 1;
        return false;
    }

    for (uint32_t i = 0; i < frags.size(); i++) {
        if (source_chr[i][current_chr] != UINT32_MAX) {
            frags[i].seek(source_chr[i][current_chr], 0);
        }
    }
    return true;
}

uint32_t MergeFragments::currentChr() const { return current_chr; }

bool MergeFragments::load() {
    // Try loading data from our existing available set
    cur_idx += capacity();

    // If we've run out of data, load and sort some more
    if (cur_idx >= available) {
        // Move leftover data to start of buffer
        uint32_t leftovers = loaded - available;
        std::memmove(start.data(), start.data() + available, sizeof(uint32_t) * leftovers);
        std::memmove(end.data(), end.data() + available, sizeof(uint32_t) * leftovers);
        std::memmove(cell.data(), cell.data() + available, sizeof(uint32_t) * leftovers);

        loaded = leftovers;
        available = 0;
        // Use min-heap of (start_coord, loader_idx) to load 2N chunks
        for (uint32_t i = 0; i < frags.size(); i++) {
            heap[i] = {frags[i].peek_start(), i};
        }
        dary_heap::make_heap<4>(heap.begin(), heap.end(), std::greater<>{});

        for (uint32_t i = 0; i < chunks_per_input * frags.size(); i++) {
            if (heap.front().first == UINT32_MAX) break; // No more data to load in this chromosome
            uint32_t idx = heap.front().second;
            if (loaded + chunk_size > start.size()) {
                throw std::runtime_error("MergeFragments: Not enough space to load input chunk");
            }
            loaded += frags[idx].load_chunk(
                start.data() + loaded, end.data() + loaded, cell.data() + loaded
            );
            // Update the chunk info and put back in heap
            heap.front().first = frags[idx].peek_start();
            dary_heap::pop_heap<4>(heap.begin(), heap.end(), std::greater<>{});
            dary_heap::push_heap<4>(heap.begin(), heap.end(), std::greater<>{});
        }
        // Sort the input data
        lsdRadixSortArrays<uint32_t, uint32_t, uint32_t>(
            loaded, start, end, cell, start_buf, end_buf, cell_buf
        );

        // Find how much of the data is guaranteed to be in the final output ordering
        // We can use the heap rather than calling peek_start() again for each fragment
        uint32_t max_sorted_start = UINT32_MAX;
        for (uint32_t i = 0; i < frags.size(); i++) {
            max_sorted_start = std::min(max_sorted_start, heap[i].first);
        }
        available = std::upper_bound(start.begin(), start.begin() + loaded, max_sorted_start) -
                    start.begin();

        cur_idx = 0;
    }

    return cur_idx < available;
}

uint32_t MergeFragments::capacity() const { return std::min(available - cur_idx, load_size); }

uint32_t *MergeFragments::cellData() { return cell.data() + cur_idx; }
uint32_t *MergeFragments::startData() { return start.data() + cur_idx; }
uint32_t *MergeFragments::endData() { return end.data() + cur_idx; }
} // namespace BPCells