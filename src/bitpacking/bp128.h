#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdio>

#include "simd_vec.h"
namespace BPCells {

// Find required number of bits per number to represent the 128 integers from in
// using standard bitpacking encoding
uint32_t simdmaxbits(const uint32_t *in);

// Find required number of bits per number to represent the 128 integers from in
// using d1 encoding with initial offset of "initvalue"
uint32_t simdmaxbitsd1(uint32_t initvalue, const uint32_t *in);

// Find required number of bits per number to represent the 128 integers from in
// using FOR bitpacking encoding.
// minvalue and bits are return values, corresponding to the number of bits
// required when the frame of reference is set to minvalue
void simdmaxbitsFOR(const uint32_t *in, uint32_t & bits, uint32_t & minvalue);

// Pack 128 integers from in to out. Packed integers use "bit" bits.
void simdpack(const uint32_t *in, uint32_t *out, const uint32_t bit);

/* reads 128 values from "in", writes  "bit" 128-bit vectors to "out".
 * The input values are assumed to be less than 1<<bit. */
// Temporarily removed since there's no need for unsafe operations when
// packing is already plenty fast enough. Unpacking speed is the main bottleneck.
// void simdpackwithoutmask(const uint32_t *in, uint32_t *out, const uint32_t bit) {
//     const vec *_in = (vec *)in;
//     vec *_out = (vec *)out;
//     switch (bit) {
//         BP128_SWITCH_CASE32(pack_nomask, _in, _out)
//     }
// }


// Unpack 128 integers from in to out. Packed integers use "bit" bits.
void simdunpack(const uint32_t *in, uint32_t *out, const uint32_t bit);

// Pack 128 integers from in to out. Packed integers use "bit" bits, and
// are d1 encoded (difference between adjacent values), with the first value
// given as an offest relative to initvalue
void simdpackd1(uint32_t initvalue, const uint32_t *in, uint32_t *out,
                   const uint32_t bit);

/* reads 128 values from "in", writes  "bit" 128-bit vectors to "out"
   integer values should be in sorted order (for best results).
   The difference values are assumed to be less than 1<<bit. */
// Temporarily removed since there's no need for unsafe operations when
// packing is already plenty fast enough. Unpacking speed is the main bottleneck.
// void simdpackwithoutmaskd1(uint32_t initvalue, const uint32_t *in, uint32_t *out,
//                            const uint32_t bit) {
//     vec init = splat(initvalue);
//     const vec *_in = (vec *)in;
//     vec *_out = (vec *)out;
//     switch (bit) {
//         BP128_SWITCH_CASE32(packd1_nomask, init, _in, _out)
//     }
// }

// Unpack 128 integers from in to out. Packed integers use "bit" bits, and
// are d1 encoded (difference between adjacent values), with the first value
// given as an offest relative to initvalue
void simdunpackd1(uint32_t initvalue, const uint32_t *in, uint32_t *out,
                     const uint32_t bit);

// Pack 128 integers from in to out. Packed integers use "bit" bits, and
// will have "initvalue" subtracted from them during packing
void simdpackFOR(uint32_t initvalue, const uint32_t *in, uint32_t *out,
                    const uint32_t bit);

// Temporarily removed since there's no need for unsafe operations when
// packing is already plenty fast enough. Unpacking speed is the main bottleneck.
// void simdpackwithoutmaskFOR(uint32_t initvalue, const uint32_t *in, uint32_t *out,
//                     const uint32_t bit) {
//     const vec *_in = (vec *)in;
//     vec *_out = (vec *)out;
//     switch (bit) {
//         BP128_SWITCH_CASE32(packFOR_nomask, initvalue, _in, _out)
//     }
// }

// Unpack 128 integers from in to out. Packed integers use "bit" bits, and
// will have "initvalue" added to them during unpacking
void simdunpackFOR(uint32_t initvalue, const uint32_t *in, uint32_t *out,
                    const uint32_t bit);

/************* TESTING FUNCTIONS FOR BITPACKING *************/

uint32_t random_uint32_t();
void fill_buf(uint32_t *buf, int bits);
bool equal_buf(uint32_t *buf1, uint32_t *buf2);

bool test_bitpacking();

} // end namespace BPCells