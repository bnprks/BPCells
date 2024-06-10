// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <algorithm>
#include <map>
#include <vector>

#include "../utils/radix_sort.h"

// Summary:
// These accumulator classes help during construction of a sparse matrix e.g. in PeakMatrix.
// You pass in non-zero matrix entries via add_one calls, then the accumulator adds together
// values at the same matrix entry and returns them via the load, and *Data functions.

// Algorithm descriptions:
// MatrixAccumulator:
//  - Store a 3 data vectors of (col, row, val). When we need to load entries, sort
//    these vectors in order of (col, row), then add together the values from identical
//    entries.
//  - Memory should not be much more than 2x the max number of non-zero entries in a column
//  - time is about 2.5 seconds to calculate peak matrix on 10k PBMCs
// MatrixAccumulator_vec
//  - For each active column, store a vector of rows containing the non-zero rows, and a
//    dense vector of values (containing zeros for unobserved rows)
//  - Memory should be approx # of cells * max # of peaks that a single fragment can overlap
//    (empirically about 10 on my 10k PBMC dataset). This makes it likely unsuitable for
//    calculating tile matrices
//  - time is about 1.6 seconds to calculate peak matrix on 10k PBMCs

namespace BPCells {

template <typename T> class MatrixAccumulator {
  private:
    std::vector<uint32_t> row_data, col_data, row_buf, col_buf;
    std::vector<T> val_data, val_buf;

    uint32_t entries_stored = 0;
    uint32_t output_idx = UINT32_MAX;
    uint32_t load_capacity = 0;

    void compactData() {
        // Sort by (col, row) by sorting stably first by row then by col
        lsdRadixSortArrays<uint32_t, uint32_t, T>(
            entries_stored, row_data, col_data, val_data, row_buf, col_buf, val_buf
        );
        lsdRadixSortArrays<uint32_t, uint32_t, T>(
            entries_stored, col_data, row_data, val_data, col_buf, row_buf, val_buf
        );

        // Add together elements belonging to the same matrix index
        uint32_t out_idx = 0;
        for (uint32_t i = 1; i < entries_stored; i++) {
            bool match = row_data[i] == row_data[out_idx] && col_data[i] == col_data[out_idx];
            out_idx += !match;
            row_data[out_idx] = row_data[i];
            col_data[out_idx] = col_data[i];
            val_data[out_idx] = match ? val_data[out_idx] + val_data[i] : val_data[i];
        }
        // If we compacted less than 50% of the data, then grow our buffer size
        if ((out_idx + 1) * 2 > row_data.size()) {
            uint32_t new_size = 1 + row_data.size() + row_data.size() / 2;
            row_data.resize(new_size);
            val_data.resize(new_size);
            col_data.resize(new_size);
            row_buf.resize(new_size);
            val_buf.resize(new_size);
            col_buf.resize(new_size);
        }
        if (entries_stored != 0) entries_stored = out_idx + 1;
    }

  public:
    void clear() {
        entries_stored = 0;
        output_idx = UINT32_MAX;
        load_capacity = 0;

        // Reset capacity down to some low but manageable size
        row_data.resize(64);
        val_data.resize(64);
        col_data.resize(64);
        row_buf.resize(64);
        val_buf.resize(64);
        col_buf.resize(64);
    }

    // Add an item to the accumulator, or no-op if val == 0
    inline void add_one(uint32_t col, uint32_t row, T val) {
        if (output_idx != UINT32_MAX) {
            std::memmove(
                row_data.data(),
                row_data.data() + output_idx,
                sizeof(uint32_t) * (entries_stored - output_idx)
            );
            std::memmove(
                col_data.data(),
                col_data.data() + output_idx,
                sizeof(uint32_t) * (entries_stored - output_idx)
            );
            std::memmove(
                val_data.data(),
                val_data.data() + output_idx,
                sizeof(T) * (entries_stored - output_idx)
            );
            entries_stored = entries_stored - output_idx;
        }
        if (entries_stored >= row_data.size()) {
            compactData();
        }
        row_data[entries_stored] = row;
        col_data[entries_stored] = col;
        val_data[entries_stored] = val;
        entries_stored += (val != 0);
        output_idx = UINT32_MAX;
        load_capacity = 0;
    }

    // Say whether the buffer is full enough to be ready for loading
    bool ready_for_loading() const {
        // Ready for loading if we have under 25% capacity for additional entries
        // The idea is to prevents repeated compactions for many small, overlapping outputs.
        // I'm not sure how important this really is in practical use cases
        return entries_stored * 4 >= row_data.size();
    }

    // Discard entries until `col`, and return whether or not the next available entry
    // is at least column `col` (or greater)
    bool discard_until(uint32_t min_col) {
        if (output_idx == UINT32_MAX) {
            compactData();
            output_idx = 0;
        } else {
            output_idx += load_capacity;
            if (output_idx == entries_stored) return false;
        }
        while (output_idx < entries_stored && col_data[output_idx] < min_col)
            output_idx += 1;
        load_capacity = 0;
        return output_idx < entries_stored && col_data[output_idx] >= min_col;
    }

    // Load accumulated entries up to max_col.
    bool load(uint32_t max_col, uint32_t max_entries) {
        if (output_idx == UINT32_MAX) {
            compactData();
            output_idx = 0;
        } else {
            output_idx += load_capacity;
            load_capacity = 0;
            if (output_idx == entries_stored) return false;
            if (col_data[output_idx] > max_col) return false;
        }
        load_capacity = std::min(max_entries, entries_stored - output_idx);
        if (col_data[output_idx] == col_data[output_idx + load_capacity - 1]) return true;
        load_capacity = std::upper_bound(col_data.begin() + output_idx, col_data.begin() + output_idx + load_capacity, col_data[output_idx]) -
            (col_data.begin() + output_idx);
        return true;
    }

    // Functions to access the loaded values
    uint32_t capacity() const { return load_capacity; }
    uint32_t currentCol() const { return col_data[output_idx]; }
    uint32_t *rowData() { return row_data.data() + output_idx; }
    T *valData() { return val_data.data() + output_idx; }
};

template <typename T> class MatrixAccumulator_vec {
  private:
    class Counter {
      public:
        uint32_t col;
        std::vector<uint32_t> rows;
        std::vector<T> vals;
    };
    std::vector<Counter> active_counters;
    std::vector<Counter> counter_pool;

    std::vector<T> val_data;
    std::vector<uint32_t> row_data;

    uint32_t current_col;
    uint32_t output_idx = UINT32_MAX;
    uint32_t load_capacity;

    uint32_t max_active_counters = 0;

  public:
    void clear() {
        output_idx = UINT32_MAX;
        active_counters.resize(0);
        counter_pool.resize(0);
    }

    // Add an item to the accumulator, or no-op if val == 0
    inline void add_one(uint32_t col, uint32_t row, T val) {
        if (val == 0) return;
        output_idx = UINT32_MAX;

        for (auto &counter : active_counters) {
            if (counter.col != col) continue;
            if (counter.vals.size() <= row) counter.vals.resize(row + 1);
            if (counter.vals[row] == 0) counter.rows.push_back(row);
            counter.vals[row] += val;
            return;
        }

        // Start a new counter
        if (counter_pool.size() == 0) counter_pool.resize(1);
        Counter &counter = counter_pool.back();
        counter.col = col;
        counter.rows.resize(0);
        counter.rows.push_back(row);
        if (counter.vals.size() <= row) counter.vals.resize(row + 1);
        counter.vals[row] += val;
        active_counters.push_back(std::move(counter));
        counter_pool.pop_back();
    }

    // Say whether the buffer is full enough to be ready for loading
    bool ready_for_loading() const { return true; }

    // Discard entries until `col`, and return whether or not the next available entry
    // is at least column `col` (or greater)
    bool discard_until(uint32_t min_col) {
        output_idx = UINT32_MAX;
        bool has_entries = false;
        for (int i = 0; i < active_counters.size(); i++) {
            if (active_counters[i].col >= min_col) {
                has_entries |= active_counters[i].rows.size() > 0;
                continue;
            }
            auto &counter = active_counters[i];
            for (auto &row : active_counters[i].rows)
                counter.vals[row] = 0;
            counter_pool.push_back(std::move(counter));
            active_counters[i] = std::move(active_counters.back());
            active_counters.pop_back();
            i--;
        }
        return has_entries;
    }

    // Load accumulated entries up to max_col (inclusive)
    bool load(uint32_t max_col, uint32_t max_entries) {
        // Handle if we have leftover entries from a completed column
        if (output_idx != UINT32_MAX && max_col >= current_col) {
            output_idx += load_capacity;
            if (output_idx < row_data.size()) {
                load_capacity = std::min((uint32_t)row_data.size() - output_idx, max_entries);
                return true;
            }
        }

        // Load a new column of entries
        // Sort active cols in decreasing order
        if (active_counters.size() == 0) return false;
        for (int i = 0; i < active_counters.size(); i++) {
            if (active_counters[i].col < active_counters.back().col)
                std::swap(active_counters[i], active_counters.back());
        }
        current_col = active_counters.back().col;
        if (current_col > max_col) return false;

        val_data.resize(0);
        std::vector<uint32_t> &row_vec = active_counters.back().rows;
        row_data.resize(row_vec.size());
        lsdRadixSortArrays(row_vec.size(), row_vec, row_data);
        std::swap(row_data, row_vec);
        for (auto &row : row_data) {
            val_data.push_back(active_counters.back().vals[row]);
            active_counters.back().vals[row] = 0;
        }
        counter_pool.push_back(std::move(active_counters.back()));
        active_counters.pop_back();

        load_capacity = std::min((uint32_t)row_data.size(), max_entries);
        output_idx = 0;
        return row_data.size() > 0;
    }

    // Functions to access the loaded values
    uint32_t capacity() const { return load_capacity; }
    uint32_t currentCol() const { return current_col; }
    uint32_t *rowData() { return row_data.data() + output_idx; }
    T *valData() { return val_data.data() + output_idx; }
};

} // end namespace BPCells