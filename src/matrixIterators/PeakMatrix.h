#pragma once

#include <algorithm>

#include "MatrixIterator.h"
#include "../fragmentIterators/FragmentIterator.h"
#include "../radix_sort.h"
#include "../arrayIO/array_interfaces.h"
#include "../bitpacking/simd_vec.h"

namespace BPCells {

class Peak {
public:
    uint32_t chr, start, end;
};

template<typename T>
class MatrixAccumulator {
private:
    std::vector<uint32_t> row_data, col_data, row_buf, col_buf;
    std::vector<T> val_data, val_buf;

    uint32_t entries_stored = 0;
    uint32_t output_idx = UINT32_MAX;
    uint32_t load_capacity = 0;

    void compactData() {
        // Sort by (col, row) by sorting stably first by row then by col
        lsdRadixSortArrays<uint32_t, T>(
            entries_stored, 
            row_data, col_data, val_data,
            row_buf, col_buf, val_buf
        );
        lsdRadixSortArrays<uint32_t, T>(
            entries_stored, 
            col_data, row_data, val_data,
            col_buf, row_buf, val_buf
        );
        // Add together elements belonging to the same matrix index
        uint32_t out_idx = 0;
        for (uint32_t i = 1; i < entries_stored; i++) {
            bool match = row_data[i] == row_data[out_idx] && 
                         col_data[i] == col_data[out_idx];
            out_idx += !match;
            row_data[out_idx] = row_data[i];
            col_data[out_idx] = col_data[i];
            val_data[out_idx] = match ? val_data[out_idx] + val_data[i] : val_data[i];
        }
        // If we compacted less than 50% of the data, then grow our buffer size
        if ((out_idx + 1) * 2 > entries_stored) {
            uint32_t new_size = 1 + entries_stored + entries_stored/2;
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
    }

    // Add an item to the accumulator, or no-op if val == 0
    inline void add_one(uint32_t col, uint32_t row, T val) {
        if (output_idx != UINT32_MAX) {
            std::memmove(row_data.data(), row_data.data() + output_idx, sizeof(uint32_t) * (entries_stored - output_idx));
            std::memmove(col_data.data(), col_data.data() + output_idx, sizeof(uint32_t) * (entries_stored - output_idx));
            std::memmove(val_data.data(), val_data.data() + output_idx, sizeof(T) * (entries_stored - output_idx));
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
        // Ready for loading if we have under 50% capacity for additional entries
        // The idea is to prevents repeated compactions for many small, overlapping outputs.
        // I'm not sure how important this really is in practical use cases
        return entries_stored*2 >= row_data.size();
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
        while(output_idx < entries_stored && col_data[output_idx] < min_col) output_idx += 1;
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
        while (col_data[output_idx] != col_data[output_idx + load_capacity - 1]) load_capacity -= 1;
        return true;
    }
    
    // Functions to access the loaded values
    uint32_t capacity() const {return load_capacity;}
    uint32_t currentCol() const {return col_data[output_idx];}
    uint32_t* rowData() {return row_data.data() + output_idx;}
    T* valData() {return val_data.data() + output_idx;}
};

// Output cell x peak matrix (rows = cell_id, col = peak_id)
// Regions are given as half-open format
// Peaks can overlap, and columns are ordered by (end, start) coordinate
class PeakMatrix : public MatrixLoader<uint32_t> {
private:
    FragmentLoader &frags;
    std::unique_ptr<StringReader> chr_levels;
    MatrixAccumulator<uint32_t> accumulator;
    std::vector<uint32_t> end_sorted_lookup; // end_sorted_lookup[i] gives the index of end-sorted peak i in the start-sorted list
    std::vector<Peak> sorted_peaks;
    std::vector<Peak> active_peaks;
    uint32_t next_completed_peak = 0;
    uint32_t current_output_peak = UINT32_MAX;
    uint32_t next_active_peak = 0;
    uint32_t n_peaks;

    void loadFragments();
public:
    // Note: It's the caller's responsibility to make sure that
    // the FragmentLoader will not be deleted while this PeakMatrix is still alive
    PeakMatrix(FragmentLoader &frags, 
        const std::vector<uint32_t> &chr, const std::vector<uint32_t> &start,
        const std::vector<uint32_t> &end, std::unique_ptr<StringReader> &&chr_levels);
    
    uint32_t rows() const override;
    uint32_t cols() const override;

    const char* rowNames(uint32_t row) const override;
    const char* colNames(uint32_t col) const override;

    void restart() override;
    void seekCol(uint32_t col) override;

    bool nextCol() override;
    
    uint32_t currentCol() const override;

    bool load() override;
    
    uint32_t capacity() const override;

    uint32_t* rowData() override;
    uint32_t* valData() override;

};

} // end namespace BPCells