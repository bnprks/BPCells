#pragma once
#include <assert.h>
#include <algorithm>
#include <vector>

#include "bp128.h"

namespace BPCells {

// Transparent vector
// Similar to std::vector<uint32_t>, but designed so I can
// point to external memory to avoid copies
typedef struct vec_uint32_t {
    uint32_t *data;
    uint32_t size;
    size_t capacity;

    inline uint32_t& operator[](uint32_t idx) { return data[idx]; }
    inline const uint32_t& operator[](uint32_t idx) const { return data[idx]; }
} vec_uint32;

typedef struct bp128_array_t {
    vec_uint32 data;
    vec_uint32 idx; //Start idx for each chunk of 128 (length chunks+1)
    uint32_t len;
} bp128_array;

typedef struct bp128_d1_array_t {
    vec_uint32 data;
    vec_uint32 idx; //Start idx for each chunk of 128 (length chunks+1)
    vec_uint32 starts;
    uint32_t len;
} bp128_d1_array;

typedef struct bp128_for_array_t {
    vec_uint32 data;
    vec_uint32 idx; //Start idx for each chunk of 128 (length chunks+1)
    vec_uint32 mins;
    uint32_t len;
} bp128_for_array;

typedef struct packed_frags_t {
    bp128_d1_array starts;
    bp128_array ends;
    // Original design considered using FOR packing for cell_ids.
    // Since we are currently just storing pos-sorted without
    // much cell chunking I'll disable the FOR encoding
    bp128_array cell_ids;
    
    vec_uint32 end_max; // Max end value across all preceding chunks of 128 including current idx

    inline uint32_t n_chunks() const {return end_max.size; }
    inline uint32_t n_frags() const {return starts.len; }
} packed_frags;


void unpack_128(const bp128_array &src, uint32_t *dst, const uint32_t idx);
void unpack_128_d1(const bp128_d1_array &src, uint32_t *dst, const uint32_t idx);
void unpack_128_for(const bp128_for_array &src, uint32_t *dst, const uint32_t idx);
void unpack_128_frags(const packed_frags &src, uint32_t *starts, uint32_t *ends, uint32_t* cell_ids, const uint32_t idx);

bp128_array pack(const uint32_t *src, const uint32_t len);
bp128_d1_array pack_d1(const uint32_t *src, const uint32_t len);
bp128_for_array pack_for(const uint32_t *src, const uint32_t len);
packed_frags pack_frags(const uint32_t *starts, 
        const uint32_t *ends, const uint32_t *cell_ids, const uint32_t len);

} // end namespace BPCells