#pragma once
#include <tuple>

// Provide a templated radix_sort method that works on 
namespace BPCells {


// Template helper to swap multiple vectors
inline void vec_swap_helper() {}

template<typename T, typename ...Rest>
inline void vec_swap_helper(std::vector<T> &a1, std::vector<Rest>&... a_tail,
                 std::vector<T> &b1, std::vector<Rest>&... b_tail) {
    std::swap(a1, b1);
    vec_swap_helper(a_tail..., b_tail...);
}

// Template helper to assign multiple vectors
inline void vec_assign_helper(uint32_t from, uint32_t to) {}

template<typename T, typename ...Rest>
inline void vec_assign_helper(uint32_t from, uint32_t to,
                 std::vector<T> &a1, std::vector<Rest>&... a_tail,
                 std::vector<T> &b1, std::vector<Rest>&... b_tail) {
    b1[to] = a1[from];
    vec_assign_helper(from, to, a_tail..., b_tail...);
}

// Do a stable radix sort on a "struct-of-arrays" oriented dataset using a uint32_t-typed key
// Must provide a buffers of equal size to the original dataset

// Example usage on fragments: sorting by end coordinate while making the corresponding
// swaps in a cell_id array.
// vector<uint32_t> ends, cell_id, end_buf, cell_buf;
// lsdRadixSortArrays<uint32_t>(ends.size(), ends, cell_id, end_buf, cell_buf);
template<typename ...Vals>
void lsdRadixSortArrays(uint32_t size, std::vector<uint32_t> &key, std::vector<Vals>&... vals, 
                                   std::vector<uint32_t> &key_scratch, std::vector<Vals>&... vals_scratch) {
    uint32_t radix_counts[4][256] = {{0}};
    bool skip_byte[4] = {false};

    // Count up how many we see of each byte combination in 1 pass of the data
    for (uint32_t i = 0; i < size; i++) {
        radix_counts[0][255 & ((uint32_t) key[i] >> (0*8))]++;
        radix_counts[1][255 & ((uint32_t) key[i] >> (1*8))]++;
        radix_counts[2][255 & ((uint32_t) key[i] >> (2*8))]++;
        radix_counts[3][255 & ((uint32_t) key[i] >> (3*8))]++;
    }

    // Tally up which output index each byte combination should start at
    // So radix_counts[i][j] will turn into the sum of all radix_counts[i][k] where k < j
    for (int i = 0; i < 4; i++) {
        uint32_t running_sum = 0;
        for (int j = 0; j < 256; j++) {
            if (radix_counts[i][j] == size) skip_byte[i] = true;
            running_sum += radix_counts[i][j];
            radix_counts[i][j] = running_sum - radix_counts[i][j];
        }
    }

    for (int i = 0; i < 4; i++) {
        // Skip if this byte is all identical
        if (skip_byte[i]) { continue; }
        for (uint32_t j = 0; j < size; j++) {
            uint32_t val = key[j];
            uint32_t bucket = 255 & ((uint32_t) val >> (i * 8));
            key_scratch[radix_counts[i][bucket]] = val;
            vec_assign_helper<Vals...>(j, radix_counts[i][bucket], vals..., vals_scratch...);
            radix_counts[i][bucket]++;
        }
        std::swap(key, key_scratch);
        vec_swap_helper<Vals...>(vals..., vals_scratch...);
    }
}

// Base case copy+paste from above, since I'm having some template trouble when no value arrays are supplied
inline void lsdRadixSortArrays(uint32_t size, std::vector<uint32_t> &key, std::vector<uint32_t> &key_scratch) {
    uint32_t radix_counts[4][256] = {{0}};
    bool skip_byte[4] = {false};

    // Count up how many we see of each byte combination in 1 pass of the data
    for (uint32_t i = 0; i < size; i++) {
        radix_counts[0][255 & ((uint32_t) key[i] >> (0*8))]++;
        radix_counts[1][255 & ((uint32_t) key[i] >> (1*8))]++;
        radix_counts[2][255 & ((uint32_t) key[i] >> (2*8))]++;
        radix_counts[3][255 & ((uint32_t) key[i] >> (3*8))]++;
    }

    // Tally up which output index each byte combination should start at
    // So radix_counts[i][j] will turn into the sum of all radix_counts[i][k] where k < j
    for (int i = 0; i < 4; i++) {
        uint32_t running_sum = 0;
        for (uint32_t j = 0; j < 256; j++) {
            if (radix_counts[i][j] == size) skip_byte[i] = true;
            running_sum += radix_counts[i][j];
            radix_counts[i][j] = running_sum - radix_counts[i][j];
        }
    }

    for (int i = 0; i < 4; i++) {
        // Skip if this byte is all identical
        if (skip_byte[i]) { continue; }
        for (uint32_t j = 0; j < size; j++) {
            uint32_t val = key[j];
            uint32_t bucket = 255 & ((uint32_t) val >> (i * 8));
            key_scratch[radix_counts[i][bucket]] = val;
            radix_counts[i][bucket]++;
        }
        std::swap(key, key_scratch);
    }
}


} // end namespace BPCells