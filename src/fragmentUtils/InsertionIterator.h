#pragma once

#include <cstring>
#include <vector>

#include "../bitpacking/simd_vec.h"
#include "../fragmentIterators/FragmentIterator.h"
#include "../utils/radix_sort.h"

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
    uint32_t end_idx = 0;
    uint32_t end_capacity = 0;
    std::vector<uint32_t> end_data, end_data_buf;
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
    inline void loadFragments() {
        // Increase buffer size if ends remaining is a large fraction of buffer size
        // This avoids having to re-sort data too frequently
        if (end_idx * 2 < end_data.size()) {
            uint32_t new_size = end_data.size() + 1 + end_data.size() / 2;
            end_data.resize(new_size);
            end_data_buf.resize(new_size);
            end_cell.resize(new_size);
            end_cell_buf.resize(new_size);
        }
        // Copy data to remove our used insertions the start of end_data
        std::memmove(&end_data[0], &end_data[end_idx], sizeof(uint32_t) * (end_capacity - end_idx));
        std::memmove(&end_cell[0], &end_cell[end_idx], sizeof(uint32_t) * (end_capacity - end_idx));
        end_capacity -= end_idx;
        end_idx = 0;

        uint32_t orig_size = end_data.size();
        uint32_t min_load = std::max(1U, (uint32_t)end_data.size() / 2);

        // Reset our vectors
        start_idx = 0;
        start_data.resize(0);
        start_cell.resize(0);
        end_data.resize(end_capacity);
        end_cell.resize(end_capacity);

        bool chr_end = false;

        while (start_data.size() < min_load) {
            if (!frags.load()) {
                chr_end = true;
                break;
            }
            uint32_t capacity = frags.capacity();
            uint32_t *starts = frags.startData();
            uint32_t *cells = frags.cellData();
            uint32_t *ends = frags.endData();
            vec one = splat(1);

            // Subtract 1 form the end coordinates
            uint32_t i;
            for (i = 0; i + 4 <= capacity; i += 4) {
                store((vec *)(ends + i), sub(load((vec *)(ends + i)), one));
            }
            for (; i < capacity; i++)
                ends[i] = ends[i] - 1;

            // Copy data into our vectors
            start_data.insert(start_data.end(), starts, starts + capacity);
            start_cell.insert(start_cell.end(), cells, cells + capacity);
            end_data.insert(end_data.end(), ends, ends + capacity);
            end_cell.insert(end_cell.end(), cells, cells + capacity);
        }

        // Resize & sort end_data
        uint32_t end_size = std::max((uint32_t)end_data.size(), orig_size);
        end_data.resize(end_size);
        end_cell.resize(end_size);
        end_data_buf.resize(end_size);
        end_cell_buf.resize(end_size);
        lsdRadixSortArrays<uint32_t, uint32_t>(
            end_capacity + start_data.size(), end_data, end_cell, end_data_buf, end_cell_buf
        );

        end_capacity += start_data.size();

        if (chr_end) {
            // Use a sentinel value of UINT32_MAX for starts, so we just need to check
            // if end_idx >= end_capacity
            start_data.push_back(UINT32_MAX);
            start_cell.push_back(UINT32_MAX);
        }
    }

  public:
    InsertionIterator(FragmentLoader &loader);

    inline void restart() {
        frags.restart();
        start_idx = 0;
        end_idx = 0;
        start_data.resize(0);
        end_capacity = 0;
    }
    inline bool nextChr() {
        start_idx = 0;
        end_idx = 0;
        bool ret = frags.nextChr();
        if (ret) current_chr = frags.currentChr();
        start_data.resize(0);
        end_capacity = 0;
        return ret;
    }
    inline bool nextInsertion() {
        if (start_idx >= start_data.size()) loadFragments();
        if (end_idx >= end_capacity) return false;
        use_start = end_data[end_idx] >= start_data[start_idx];
        next_coord = use_start ? start_data[start_idx] : end_data[end_idx];
        next_cell = use_start ? start_cell[start_idx] : end_cell[end_idx];
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
        end_capacity = 0;
    }

    inline uint32_t chr() const { return current_chr; };
    inline uint32_t cell() const { return next_cell; };
    inline uint32_t coord() const { return next_coord; };
    inline bool isStart() const { return use_start; };
};

} // end namespace BPCells