#pragma once
#include "../bitpacking/simd_vec.h"
#include <algorithm>

// Helper class to keep together the data and scratch vectors
template <typename T> using VecRefPair = std::pair<std::vector<T> &, std::vector<T> &>;

// Template helper to swap multiple vectors
inline void vec_swap_helper() {}

template <typename T, typename... Rest>
inline void vec_swap_helper(VecRefPair<T> head, VecRefPair<Rest>... tail) {
    std::swap(head.first, head.second);
    vec_swap_helper(tail...);
}

// Template helper to assign multiple vectors
inline void vec_assign_helper(uint32_t from, uint32_t to) {}

template <typename T, typename... Rest>
inline void
vec_assign_helper(uint32_t from, uint32_t to, VecRefPair<T> head, VecRefPair<Rest>... tail) {
    head.second[to] = head.first[from];
    vec_assign_helper(from, to, tail...);
}

// Do a stable radix sort on a "struct-of-arrays" oriented dataset using a uint32_t-typed key
// Must provide a buffers of equal size to the original dataset

// Example usage on fragments: sorting by end coordinate while making the corresponding
// swaps in a cell_id array.
// vector<uint32_t> ends, cell_id, end_buf, cell_buf;
// lsdRadixSortArrays<uint32_t>(ends.size(), {ends, end_buf}, {cell_id, cell_buf});
template <typename... Vals>
void lsdRadixSortArrays(
    uint32_t size, uint32_t *key, uint32_t *key_scratch, VecRefPair<Vals>... vals
) {
    uint32_t radix_counts[4][256] = {{0}};
    bool skip_byte[4] = {false};

    // Count up how many we see of each byte combination in 1 pass of the data
    for (uint32_t i = 0; i < size; i++) {
        radix_counts[0][255 & ((uint32_t)key[i] >> (0 * 8))]++;
        radix_counts[1][255 & ((uint32_t)key[i] >> (1 * 8))]++;
        radix_counts[2][255 & ((uint32_t)key[i] >> (2 * 8))]++;
        radix_counts[3][255 & ((uint32_t)key[i] >> (3 * 8))]++;
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
        if (skip_byte[i]) {
            continue;
        }
        for (uint32_t j = 0; j < size; j++) {
            uint32_t val = key[j];
            uint32_t bucket = 255 & ((uint32_t)val >> (i * 8));
            key_scratch[radix_counts[i][bucket]] = val;
            vec_assign_helper<uint32_t, Vals...>(j, radix_counts[i][bucket], key, vals...);
            // vec_assign_helper<Vals...>(j, radix_counts[i][bucket], vals..., vals_scratch...);
            radix_counts[i][bucket]++;
        }
        // std::swap(key, key_scratch);
        vec_swap_helper<uint32_t, Vals...>(key, vals...);
        // vec_swap_helper<Vals...>(vals..., vals_scratch...);
    }
}

template <typename... Vals>
void lsdRadixSortArrays(uint32_t size, VecRefPair<uint32_t> key, VecRefPair<Vals>... vals) {}