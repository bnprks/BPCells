// Copyright 2024 BPCells contributors
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once
#include "../fragmentIterators/FragmentIterator.h"
#include "../simd/math.h"
#include "../utils/radix_sort.h"
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

namespace BPCells {

// Transform sorted fragments into sorted insertions -- i.e. merge
// the starts + ends into a sorted stream of insertions. Note that since
// fragments store end coordinates non-inclusive, end coordinates will be shifted
// down by 1bp
class InsertionIterator {
  private:
    FragmentLoader &frags;

    // Pointers for accessing sorted start coordinate data
    uint32_t start_idx = 0;
    std::vector<uint32_t> start_data;
    std::vector<uint32_t> start_cell;

    // Space for storing + sorting end coordinate data (shifted by 1)
    size_t end_idx = 0;
    // Current index of the already loaded data; set to max for no data loaded
    size_t load_idx = std::numeric_limits<size_t>::max();
    // Vectors for end cordinate data
    std::vector<uint32_t> end_data, end_data_buf;
    // Vectors for end cell data
    std::vector<uint32_t> end_cell, end_cell_buf;

    uint32_t next_cell;
    uint32_t next_coord;
    uint32_t current_chr;
    bool use_start;

    // Load new fragments, shift end coordinates, and add to the sorted end_data buf
    // Pre-condition:
    // - We've iterated through all the insertions in start_data + start_cell, and
    //   there are some insertions left in end_data
    // Post-condition:
    // - start_data + start_cell have new insertions, and end_data + end_cell have
    //   the corresponding end insertions added + sorted
    // - Unless we have run out of insertions to load, end_data has been filled all
    //   the way to end_data.capacity(), and same for end_cell
    inline void loadFragments() {
        uint32_t leftover_size = end_data.size() - end_idx;

        // Increase buffer size if ends remaining is a large fraction of buffer size
        // This avoids having to re-sort data too frequently
        if (leftover_size * 2 > end_data.size()) {
            // Increase capacity by at least 1.5x the current size (may be more depending on
            // std::vector implementation)
            uint32_t new_size = end_data.size() + 1 + end_data.size() / 2;
            end_data.resize(new_size);
            end_cell.resize(end_data.capacity());
            end_data_buf.resize(end_data.capacity());
            end_cell_buf.resize(end_cell.capacity());
        }

        // Copy data to remove our used insertions the start of end_data
        std::memmove(end_data.data(), end_data.data() + end_idx, sizeof(uint32_t) * leftover_size);
        std::memmove(end_cell.data(), end_cell.data() + end_idx, sizeof(uint32_t) * leftover_size);

        // Reset our vectors
        end_idx = 0;
        end_data.resize(leftover_size);
        end_cell.resize(leftover_size);
        start_idx = 0;
        start_data.resize(0);
        start_cell.resize(0);

        bool chr_end = false;

        // Copy any leftovers from the previous `frags.load()`
        if (load_idx != std::numeric_limits<size_t>::max() && load_idx < frags.capacity()) {
            uint32_t load_size =
                std::min(end_data.capacity() - end_data.size(), frags.capacity() - load_idx);
            uint32_t *starts = frags.startData() + load_idx;
            uint32_t *cells = frags.cellData() + load_idx;
            uint32_t *ends = frags.endData() + load_idx;
            start_data.insert(start_data.end(), starts, starts + load_size);
            start_cell.insert(start_cell.end(), cells, cells + load_size);
            end_data.insert(end_data.end(), ends, ends + load_size);
            end_cell.insert(end_cell.end(), cells, cells + load_size);
            load_idx += load_size;
        }

        // Load data until we have exactly hit end_data.capacity()
        while (end_data.size() < end_data.capacity()) {
            if (!frags.load()) {
                chr_end = true;
                load_idx = std::numeric_limits<size_t>::max();
                break;
            }
            load_idx = std::min((size_t)frags.capacity(), end_data.capacity() - end_data.size());
            uint32_t *starts = frags.startData();
            uint32_t *cells = frags.cellData();
            uint32_t *ends = frags.endData();

            // Subtract 1 form the end coordinates
            simd::add(ends, -1, frags.capacity());

            // Copy data into our vectors
            start_data.insert(start_data.end(), starts, starts + load_idx);
            start_cell.insert(start_cell.end(), cells, cells + load_idx);
            end_data.insert(end_data.end(), ends, ends + load_idx);
            end_cell.insert(end_cell.end(), cells, cells + load_idx);
        }

        // While start sites are guaranteed to be sorted, end sites are not, so sort with a buffer
        lsdRadixSortArrays<uint32_t, uint32_t>(
            end_data.size(), end_data, end_cell, end_data_buf, end_cell_buf
        );

        if (chr_end) {
            // Use a sentinel value of UINT32_MAX for starts, so we just need to check
            // if end_idx >= end_capacity
            start_data.push_back(UINT32_MAX);
            start_cell.push_back(UINT32_MAX);
        }
    }

  public:
    InsertionIterator(FragmentLoader &loader);

    inline bool nextChr() {
        bool ret = frags.nextChr();
        if (ret) current_chr = frags.currentChr();
        start_idx = 0;
        end_idx = 0;
        start_data.resize(0);
        end_data.resize(0);
        load_idx = std::numeric_limits<size_t>::max();
        return ret;
    }

    inline void restart() {
        frags.restart();
        start_idx = 0;
        end_idx = 0;
        start_data.resize(0);
        end_data.resize(0);
        load_idx = std::numeric_limits<size_t>::max();
    }

    inline bool nextInsertion() {
        if (start_idx >= start_data.size()) loadFragments();
        if (end_idx >= end_data.size()) return false;
        // Prioritize clearing start data when end data is equal (pragmatically doesn't change
        // bedgraph output)
        use_start = end_data[end_idx] >= start_data[start_idx];
        next_coord = use_start ? start_data[start_idx] : end_data[end_idx];
        next_cell = use_start ? start_cell[start_idx] : end_cell[end_idx];
        // advance to next start idx if we use a start coord, else advance next end idx
        start_idx += use_start;
        end_idx += !use_start;

        return true;
    }

    inline void seek(uint32_t chr_id, uint32_t base) {
        frags.seek(chr_id, base);
        current_chr = frags.currentChr();
        start_idx = 0;
        end_idx = 0;
        start_data.resize(0);
        end_data.resize(0);
        load_idx = std::numeric_limits<size_t>::max();
    }

    inline uint32_t chr() const { return current_chr; };
    inline uint32_t cell() const { return next_cell; };
    inline uint32_t coord() const { return next_coord; };
    inline bool isStart() const { return use_start; };
};

} // end namespace BPCells