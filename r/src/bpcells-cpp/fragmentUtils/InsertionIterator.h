// // Copyright 2022 BPCells contributors
// // 
// // Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// // https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// // <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// // option. This file may not be copied, modified, or distributed
// // except according to those terms.
// 
// #pragma once
// #include <cstring>
// #include <vector>
// #include <queue>
// #include <cstdint>
// #include "../fragmentIterators/FragmentIterator.h"
// #include "../utils/radix_sort.h"
// #include "../simd/math.h"
// 
// struct Insertion {
//     uint32_t coord;
//     uint32_t cell;
//     bool is_start; // True if this is a start coordinate, false if it's an end coordinate
//     uint32_t chrom;
//     Insertion(uint32_t coord, uint32_t cell, bool is_start, uint32_t chrom) : coord(coord), cell(cell), is_start(is_start), chrom(chrom) {}
// };
// struct CompareInsertion {
//     bool operator()(const Insertion &a, const Insertion &b) {
//         return a.coord >= b.coord;
//     }
// };
// 
// namespace BPCells {
// 
// // Transform sorted fragments into sorted insertions -- i.e. merge
// // the starts + ends into a sorted stream of insertions. Note that since
// // fragments store end coordinates non-inclusive, end coordinates will be shifted
// // down by 1bp
// class InsertionIterator {
//   private:
//     FragmentLoader &frags;
//     std::priority_queue<Insertion, std::vector<Insertion>, CompareInsertion> frag_queue;
//     uint32_t highest_start;
// 
//     uint32_t next_cell;
//     uint32_t next_coord;
//     uint32_t current_chr;
//     bool use_start;
// 
//     // Load new fragments, shift end coordinates, and add to the sorted end_data buf
//     // Pre-condition:
//     // - We've iterated through all the insertions in start_data + start_cell, and
//     //   there are some insertions left in end_data
//     // Post-condition:
//     // - start_data + start_cell have new insertions, and end_data + end_cell have
//     //   the corresponding end insertions added + sorted
//     inline bool loadFragments() {
//         if (!frags.load()) {
//             return false;
//         }
//         uint32_t* starts = frags.startData();
//         uint32_t* ends = frags.endData();
//         uint32_t* cells = frags.cellData();
// 
//         // Insert start and end fragments into the priority queue
//         for (uint32_t i = 0; i < frags.capacity(); ++i) {
//             frag_queue.push(Insertion(starts[i], cells[i], true, current_chr));
//             highest_start = std::max(highest_start, starts[i]);
//             frag_queue.push(Insertion(ends[i] - 1, cells[i], false, current_chr));
//         }
//         return true;
//     }
//   public:
// 
//     InsertionIterator(FragmentLoader &loader);
// 
// 
//     inline bool nextChr() {
//         bool ret = frags.nextChr();
//         if (ret) {
//             current_chr = frags.currentChr();
//             std::priority_queue<Insertion, std::vector<Insertion>, CompareInsertion> emp;
//             std::swap(frag_queue, emp);
//             highest_start = 0;
//         }
//         return ret;
//     }
// 
//     inline void restart() {
//         frags.restart();
//         std::priority_queue<Insertion, std::vector<Insertion>, CompareInsertion> emp;
//         std::swap(frag_queue, emp);
//     }
// 
//     inline bool nextInsertion(bool rec = false) {
//         if (frag_queue.empty() && !rec) {
//             bool ret = loadFragments();
//             if (!ret) return false;
//         }
//         Insertion next_frag = frag_queue.top();
//         // Load fragments when at the start sites with the highest coord remaining in the queue
//         // This is to ensure that we don't miss any insertions that start before the largest end coord in
//         // the queue.
//         while ((next_frag.is_start) && (next_frag.coord >= highest_start)) {
//             bool ret = loadFragments();
//             if (ret) {
//                 next_frag = frag_queue.top();
//             } else {
//                 break;
//             }
//         }
//         next_frag = frag_queue.top();
//         frag_queue.pop();
//         next_cell = next_frag.cell;
//         next_coord = next_frag.coord;
//         use_start = next_frag.is_start;
//         return true;
//     }
// 
//     inline void seek(uint32_t chr_id, uint32_t base) {
//         frags.seek(chr_id, base);
//         current_chr = frags.currentChr();
//         highest_start = 0;
//         std::priority_queue<Insertion, std::vector<Insertion>, CompareInsertion> emp;
//         std::swap(frag_queue, emp);
//     }
// 
//     inline uint32_t chr() const { return current_chr; };
//     inline uint32_t cell() const { return next_cell; };
//     inline uint32_t coord() const { return next_coord; };
//     inline bool isStart() const { return use_start; };
// };
// 
// } // end namespace BPCells

// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <cstring>
#include <vector>

#include "../fragmentIterators/FragmentIterator.h"
#include "../utils/radix_sort.h"
#include "../simd/math.h"

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

            // Subtract 1 form the end coordinates
            simd::add(ends, -1, capacity);

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