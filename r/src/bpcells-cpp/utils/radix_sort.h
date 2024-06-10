// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once
#include <cstdint>
#include <utility>
#include <vector>

// Provide a templated radix_sort method that works on
namespace BPCells {


// Template helper to swap multiple vectors
template <typename... Ts>
inline void vec_swap_helper(
    std::pair<std::vector<Ts> &, std::vector<Ts> &> ...vecs
) {
    (std::swap(vecs.first, vecs.second), ...);
}

// Template helper to assign multiple vectors
template <typename... Ts>
inline void vec_assign_helper(
    uint32_t from,
    uint32_t to,
    std::pair<std::vector<Ts> &, std::vector<Ts> &> ...vecs
) {
    // Behold my C++ template-foo! This is called a fold expression
    ((vecs.second[to] = vecs.first[from]), ...);
}

// Helper template logic for converting to a uint32_t or uint64_t radix
template <typename T>
struct RadixType
{
};

template <>
struct RadixType<uint32_t>
{
    using type = uint32_t;
};

template <>
struct RadixType<int32_t>
{
    using type = uint32_t;
};

template <>
struct RadixType<float>
{
    using type = uint32_t;
};

template <>
struct RadixType<uint64_t>
{
    using type = uint64_t;
};

template <>
struct RadixType<int64_t>
{
    using type = uint64_t;
};

template <>
struct RadixType<double>
{
    using type = uint64_t;
};

template <typename T, typename Radix = typename RadixType<T>::type>
struct RadixTransform
{
    static void preprocessRadixKey(Radix *key, uint32_t size);
    static void postprocessRadixKey(Radix *key, uint32_t size);
};

template <>
struct RadixTransform<uint32_t, uint32_t>
{
    static void preprocessRadixKey(uint32_t *key, uint32_t size){};
    static void postprocessRadixKey(uint32_t *key, uint32_t size){};
};

// For ints, flip the sign bit
template <>
struct RadixTransform<int32_t, uint32_t>
{
    static void preprocessRadixKey(uint32_t *key, uint32_t size)
    {
        for (uint32_t i = 0; i < size; i++)
        {
            key[i] = key[i] ^ ((uint32_t)1 << 31);
        }
    };
    static void postprocessRadixKey(uint32_t *key, uint32_t size)
    {
        for (uint32_t i = 0; i < size; i++)
        {
            key[i] = key[i] ^ ((uint32_t)1 << 31);
        }
    };
};

template <>
struct RadixTransform<float, uint32_t>
{
    static void preprocessRadixKey(uint32_t *key, uint32_t size)
    {
        for (uint32_t i = 0; i < size; i++)
        {
            // If sign bit is 1, flip all bits. Otherwise just flip sign bit
            uint32_t mask = (0 - (key[i] >> 31)) | ((uint32_t)1 << 31);
            key[i] ^= mask;
        }
    };
    static void postprocessRadixKey(uint32_t *key, uint32_t size)
    {
        for (uint32_t i = 0; i < size; i++)
        {
            // If sign bit is 0, flip all bits. Otherwise just flip sign bit
            uint32_t mask = (UINT32_MAX + (key[i] >> 31)) | ((uint32_t)1 << 31);
            key[i] ^= mask;
        }
    };
};

template <>
struct RadixTransform<uint64_t, uint64_t>
{
    static void preprocessRadixKey(uint64_t *key, uint32_t size){};
    static void postprocessRadixKey(uint64_t *key, uint32_t size){};
};

// For ints, flip the sign bit
template <>
struct RadixTransform<int64_t, uint64_t>
{
    static void preprocessRadixKey(uint64_t *key, uint32_t size)
    {
        for (uint32_t i = 0; i < size; i++)
        {
            key[i] = key[i] ^ ((uint64_t)1 << 63);
        }
    };
    static void postprocessRadixKey(uint64_t *key, uint32_t size)
    {
        for (uint32_t i = 0; i < size; i++)
        {
            key[i] = key[i] ^ ((uint64_t)1 << 63);
        }
    };
};

template <>
struct RadixTransform<double, uint64_t>
{
    static void preprocessRadixKey(uint64_t *key, uint32_t size)
    {
        for (uint32_t i = 0; i < size; i++)
        {
            // If sign bit is 1, flip all bits. Otherwise just flip sign bit
            uint64_t mask = (0 - (key[i] >> 63)) | ((uint64_t)1 << 63);
            key[i] ^= mask;
        }
    };
    static void postprocessRadixKey(uint64_t *key, uint32_t size)
    {
        for (uint32_t i = 0; i < size; i++)
        {
            // If sign bit is 0, flip all bits. Otherwise just flip sign bit
            uint64_t mask = (UINT64_MAX + (key[i] >> 63)) | ((uint64_t)1 << 63);
            key[i] ^= mask;
        }
    };
};

// Do a stable radix sort on a "struct-of-arrays" oriented dataset using a uint32_t-typed key
// Must provide a buffers of equal size to the original dataset

// Example usage on fragments: sorting by end coordinate while making the corresponding
// swaps in a cell_id array.
// vector<uint32_t> ends, cell_id, end_buf, cell_buf;
// lsdRadixSortArrays<uint32_t>(ends.size(), {ends, end_buf}, {cell_id, cell_buf});
template <typename T, typename... Vals>
void lsdRadixSortArrays(
    uint32_t size,
    std::vector<T> &key,
    std::vector<Vals> &...vals,
    std::vector<T> &key_scratch,
    std::vector<Vals> &...vals_scratch
) 
{
    static_assert(sizeof(T) == 4 || sizeof(T) == 8);
    using Radix = typename RadixType<T>::type;
    const uint32_t radix_bytes = sizeof(T);

    uint32_t radix_counts[radix_bytes][256] = {{0}};
    bool skip_byte[radix_bytes] = {false};

    Radix *key_data = (Radix *)key.data();
    RadixTransform<T>::preprocessRadixKey(key_data, size);

    // Count up how many we see of each byte combination in 1 pass of the data
    if constexpr (radix_bytes == 4)
    {
        for (uint32_t i = 0; i < size; i++)
        {
            radix_counts[0][255 & (key_data[i] >> (0 * 8))]++;
            radix_counts[1][255 & (key_data[i] >> (1 * 8))]++;
            radix_counts[2][255 & (key_data[i] >> (2 * 8))]++;
            radix_counts[3][255 & (key_data[i] >> (3 * 8))]++;
        }
    }
    if constexpr (radix_bytes == 8)
    {
        for (uint32_t i = 0; i < size; i++)
        {
            radix_counts[0][255 & (key_data[i] >> (0 * 8))]++;
            radix_counts[1][255 & (key_data[i] >> (1 * 8))]++;
            radix_counts[2][255 & (key_data[i] >> (2 * 8))]++;
            radix_counts[3][255 & (key_data[i] >> (3 * 8))]++;
            radix_counts[4][255 & (key_data[i] >> (4 * 8))]++;
            radix_counts[5][255 & (key_data[i] >> (5 * 8))]++;
            radix_counts[6][255 & (key_data[i] >> (6 * 8))]++;
            radix_counts[7][255 & (key_data[i] >> (7 * 8))]++;
        }
    }

    // Tally up which output index each byte combination should start at
    // So radix_counts[i][j] will turn into the sum of all radix_counts[i][k] where k < j
    for (uint32_t i = 0; i < radix_bytes; i++)
    {
        uint32_t running_sum = 0;
        for (int j = 0; j < 256; j++)
        {
            if (radix_counts[i][j] == size)
                skip_byte[i] = true;
            running_sum += radix_counts[i][j];
            radix_counts[i][j] = running_sum - radix_counts[i][j];
        }
    }

    for (uint32_t i = 0; i < radix_bytes; i++)
    {
        // Skip if this byte is all identical
        if (skip_byte[i])
        {
            continue;
        }
        for (uint32_t j = 0; j < size; j++)
        {
            uint32_t bucket = 255 & (key_data[j] >> (i * 8));
            vec_assign_helper<T, Vals...>(j, radix_counts[i][bucket], {key, key_scratch}, {vals, vals_scratch}...);
            radix_counts[i][bucket]++;
        }
        vec_swap_helper<T, Vals...>({key, key_scratch}, {vals, vals_scratch}...);
        key_data = (Radix *)key.data();
    }
    RadixTransform<T>::postprocessRadixKey(key_data, size);
}

} // end namespace BPCells