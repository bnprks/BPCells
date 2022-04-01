#pragma once

#include <algorithm>
#include <map>
#include <vector>

#include "../radix_sort.h"
// #include "../lib/robin_hood.h"
// #include "../lib/tsl/robin_map.h"

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
// 
//
// MatrixAccumulator_ordered_map, MatrixAccumulator_tsl_hash
//  - Attempts at using ordered or unordered maps to add together values from the same matrix entry.
//    Uniformly slower, on the order of 

namespace BPCells {

template<typename T>
class MatrixAccumulator {
private:
    std::vector<uint32_t> row_data, col_data, row_buf, col_buf;
    std::vector<T> val_data, val_buf;

    uint32_t entries_stored = 0;
    uint32_t output_idx = UINT32_MAX;
    uint32_t load_capacity = 0;

    uint32_t sorted_count = 0;

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
        sorted_count += entries_stored;
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
        if ((out_idx + 1) * 2 > row_data.size()) {
            uint32_t new_size = 1 + row_data.size() + row_data.size()/2;
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
    ~MatrixAccumulator() {
        printf("Sorted_count=%d\n", sorted_count);
    }

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



template<typename T>
class MatrixAccumulator_ordered_map {
private:
    std::map<std::pair<uint32_t, uint32_t>, T> values; 

    uint32_t current_col;
    std::vector<uint32_t> row_data;
    std::vector<T> val_data;
public:

    void clear() {
        values.clear();
    }

    // Add an item to the accumulator, or no-op if val == 0
    inline void add_one(uint32_t col, uint32_t row, T val) {
        if (val == 0) return;
        
        auto res = values.emplace(std::pair<uint32_t, uint32_t>{col, row}, val);
        if (!res.second) {
            res.first->second += val;
        }
    }

    // Say whether the buffer is full enough to be ready for loading
    bool ready_for_loading() const {return true;}

    // Discard entries until `col`, and return whether or not the next available entry
    // is at least column `col` (or greater)
    bool discard_until(uint32_t min_col) {
        auto final_it = values.lower_bound({min_col, 0});
        values.erase(values.begin(), final_it);
        return values.size() > 0;
    }

    // Load accumulated entries up to max_col (inclusive) 
    bool load(uint32_t max_col, uint32_t max_entries) {
        if (values.size() == 0) return false;
        auto it = values.begin();
        current_col = values.begin()->first.first;
        if (current_col > max_col) return false;
        row_data.resize(0); val_data.resize(0);

        for (;it != values.end(); ++it) {
            if (it->first.first != current_col) break;
            if (row_data.size() >= max_entries) break;
            row_data.push_back(it->first.second);
            val_data.push_back(it->second);
        }
        values.erase(values.begin(), it);
        return row_data.size() > 0;
    }
    
    // Functions to access the loaded values
    uint32_t capacity() const {return row_data.size();}
    uint32_t currentCol() const {return current_col;}
    uint32_t* rowData() {return row_data.data();}
    T* valData() {return val_data.data();}
};




// template<typename T>
// class MatrixAccumulator_tsl_hash {
// private:
//     // Map from (col, row) -> value
//     robin_hood::unordered_flat_map<uint64_t, T> values; 
    
//     // Map from col -> vector of rows
//     robin_hood::unordered_flat_map<uint32_t, std::vector<uint32_t>> rows; 
    
//     // Heap of columns with some entries present
//     std::vector<uint32_t> cols;

//     std::vector<std::vector<uint32_t>> vec_pool;

//     uint32_t current_col;
//     std::vector<uint32_t> row_data;
//     std::vector<T> val_data;

//     uint32_t output_idx = UINT32_MAX;
//     uint32_t load_capacity;

//     static inline uint64_t make_key(uint32_t col, uint32_t row) {
//         return (((uint64_t) col) << 32) | row;
//     }
// public:

//     void clear() {
//         values.clear();
//     }

//     // Add an item to the accumulator, or no-op if val == 0
//     inline void add_one(uint32_t col, uint32_t row, T val) {
//         if (val == 0) return;
//         output_idx = UINT32_MAX;
//         // Add value
//         auto res = values.emplace(make_key(col, row), val);
//         if (!res.second) {
//             res.first->second += val;
//             return;
//         }

//         // Record the row is non-empty, taking a vector from the object pool to use
//         // as the value
//         if (vec_pool.size() == 0) vec_pool.resize(1);
//         auto res2 = rows.emplace(col, std::move(vec_pool.back()));        
//         if (res2.second) {
//             vec_pool.pop_back();
//             // Record the column as non-empty
//             cols.push_back(col); std::push_heap(cols.begin(), cols.end(), std::greater<>{});
//         }
//         res2.first->second.push_back(row);
//     }

//     // Say whether the buffer is full enough to be ready for loading
//     bool ready_for_loading() const {return true;}

//     // Discard entries until `col`, and return whether or not the next available entry
//     // is at least column `col` (or greater)
//     bool discard_until(uint32_t min_col) {
//         while (cols.front() < min_col) {
//             // Get vector of all entries in the column
//             for (auto row : rows.at(cols.front())) {
//                 values.erase(make_key(cols.front(), row));
//             }
//             rows.erase(cols.front());
//             std::pop_heap(cols.begin(), cols.end(), std::greater<>{}); cols.pop_back();
//         }
//         return values.size() > 0;
//     }

//     // Load accumulated entries up to max_col (inclusive) 
//     bool load(uint32_t max_col, uint32_t max_entries) {
//         // Handle if we have leftover entries from a completed column
//         if (output_idx != UINT32_MAX && max_col >= current_col) {
//             output_idx += load_capacity; 
//             if (output_idx < row_data.size()) {
//                 load_capacity = std::min((uint32_t) row_data.size() - output_idx, max_entries);
//                 return true;
//             }
//         }

//         // Load a new column of entries
//         if (values.size() == 0) return false;
//         current_col = cols.front();
//         if (current_col > max_col) return false;
//         std::pop_heap(cols.begin(), cols.end(), std::greater<>{});
//         cols.pop_back();
        
//         val_data.resize(0);

//         std::vector<uint32_t> &row_vec = rows.at(current_col);
//         row_data.resize(row_vec.size());
//         lsdRadixSortArrays(row_vec.size(), row_vec, row_data);
//         std::swap(row_data, row_vec);
//         for (auto row : row_data) {
//             auto it = values.find(make_key(current_col, row));
//             val_data.push_back(it->second);
//             values.erase(it);
//         }
//         row_vec.resize(0);
//         vec_pool.push_back(std::move(row_vec));
//         rows.erase(current_col);
//         load_capacity = std::min((uint32_t) row_data.size(), max_entries);
//         output_idx = 0;
//         return row_data.size() > 0;
//     }
    
//     // Functions to access the loaded values
//     uint32_t capacity() const {return load_capacity;}
//     uint32_t currentCol() const {return current_col;}
//     uint32_t* rowData() {return row_data.data() + output_idx;}
//     T* valData() {return val_data.data() + output_idx;}
// };



template<typename T>
class MatrixAccumulator_vec {
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
            if (counter.vals.size() <= row) counter.vals.resize(row+1);
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
        if (counter.vals.size() <= row) counter.vals.resize(row+1);
        counter.vals[row] += val;
        active_counters.push_back(std::move(counter));
        counter_pool.pop_back();
    }

    // Say whether the buffer is full enough to be ready for loading
    bool ready_for_loading() const {return true;}

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
            for (auto row : active_counters[i].rows)
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
                load_capacity = std::min((uint32_t) row_data.size() - output_idx, max_entries);
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
        for (auto row : row_data) {
            val_data.push_back(active_counters.back().vals[row]);
            active_counters.back().vals[row] = 0;
        }
        counter_pool.push_back(std::move(active_counters.back()));
        active_counters.pop_back();
        
        load_capacity = std::min((uint32_t) row_data.size(), max_entries);
        output_idx = 0;
        return row_data.size() > 0;
    }
    
    // Functions to access the loaded values
    uint32_t capacity() const {return load_capacity;}
    uint32_t currentCol() const {return current_col;}
    uint32_t* rowData() {return row_data.data() + output_idx;}
    T* valData() {return val_data.data() + output_idx;}
};

} // end namespace BPCells