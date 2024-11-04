// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <cstdint>
#include <hwy/base.h>

namespace BPCells::simd::bp128 {


// ##############################################################################
// #                           REFERENCE RESOURCES                              #
// ##############################################################################
// - Original BP-128 paper by Daniel Lemire: https://arxiv.org/pdf/1209.2137
// - Check out pack128 and unpack128 from the tests in `test-bp128.cpp` for an easy-to-read implementation
//   of the vanilla BP-128 logic
// - Daniel Lemire repo with some of these functions implemented: https://github.com/lemire/simdcomp
//    - A bit tough to read, but was an important reference for the original implementation


/**
 * @brief Unpack 128 bitpacked 32-bit integers
 * 
 * @param in pointer to bitpacked integers
 * @param out pointer to output destination (must be non-overlapping with input!)
 * @param bit number of bits per integer
 */
void unpack(const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t bit);
/**
 * @brief Bitpack 128 32-bit integers
 * 
 * @param in pointer to input integers
 * @param out pointer to output destination (must be non-overlapping with input!)
 * @param bit number of bits per integer
 */
void pack(const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t bit);
/**
 * @brief Determine the number of bits needed to represent the largest of 128 32-bit integers
 * 
 * @param in pointer to input integers
 * @return uint32_t number of bits needed per integer for bitpacking
 */
uint32_t maxbits(const uint32_t *HWY_RESTRICT in);

/**
 * @brief Unpack 128 bitpacked 32-bit integers with d1 encoding
 * 
 * For d1 encoding, take the difference between adjacent inputs prior to packing
 * 
 * @param initvalue Value to be subtracted from the fist input
 * @param in pointer to input integers
 * @param out pointer to output destination (must be non-overlapping with input!)
 * @param bit number of bits per integer
 */
void unpack_d1(uint32_t initvalue, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t bit);
/**
 * @brief Bitpack 128 32-bit integers with d1 encoding
 * 
 * For d1 encoding, take the difference between adjacent inputs prior to packing
 * 
 * @param initvalue Value to be subtracted from the fist input
 * @param in pointer to input integers
 * @param out pointer to output destination (must be non-overlapping with input!)
 * @param bit number of bits per integer
 */
void pack_d1(uint32_t initvalue, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t bit);
/**
 * @brief Determine the number of bits needed to represent the largest of 128 32-bit integers with d1 encoding
 * 
 * For d1 encoding, take the difference between adjacent inputs prior to packing
 * 
 * @param initvalue Value to be subtracted from the fist input
 * @param in pointer to input integers
 * @return uint32_t number of bits needed per integer for bitpacking
 */
uint32_t maxbits_d1(uint32_t initvalue, const uint32_t *HWY_RESTRICT in);

/**
 * @brief Unpack 128 bitpacked 32-bit integers with d1z encoding
 * 
 * For d1z encoding, take the difference between adjacent inputs, then
 * zigzag encode (0,-1,1,-2,2...) -> (0,1,2,3,4,...)
 * 
 * @param initvalue Value to be subtracted from the fist input
 * @param in pointer to input integers
 * @param out pointer to output destination (must be non-overlapping with input!)
 * @param bit number of bits per integer
 */
void unpack_d1z(uint32_t initvalue, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t bit);
/**
 * @brief Bitpack 128 32-bit integers with d1z encoding
 * 
 * For d1z encoding, take the difference between adjacent inputs, then
 * zigzag encode (0,-1,1,-2,2...) -> (0,1,2,3,4,...)
 * 
 * @param initvalue Value to be subtracted from the fist input
 * @param in pointer to input integers
 * @param out pointer to output destination (must be non-overlapping with input!)
 * @param bit number of bits per integer
 */
void pack_d1z(uint32_t initvalue, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t bit);
/**
 * @brief Determine the number of bits needed to represent the largest of 128 32-bit integers with d1z encoding
 * 
 * For d1z encoding, take the difference between adjacent inputs, then
 * zigzag encode (0,-1,1,-2,2...) -> (0,1,2,3,4,...)
 * 
 * @param initvalue Value to be subtracted from the fist input
 * @param in pointer to input integers
 * @return uint32_t number of bits needed per integer for bitpacking
 */
uint32_t maxbits_d1z(uint32_t initvalue, const uint32_t *HWY_RESTRICT in);

/**
 * @brief Unpack 128 bitpacked 32-bit integers with frame of reference encoding
 * 
 * For frame of reference encoding, subtract a constant from inputs prior to packing
 * 
 * @param offset Value to apply as offset
 * @param in pointer to input integers
 * @param out pointer to output destination (must be non-overlapping with input!)
 * @param bit number of bits per integer
 */
void unpack_FOR(uint32_t offset, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t bit);
/**
 * @brief Bitpack 128 32-bit integers with frame of reference encoding
 * 
 * For frame of reference encoding, subtract a constant from inputs prior to packing
 * 
 * @param offset Value to apply as offset
 * @param in pointer to input integers
 * @param out pointer to output destination (must be non-overlapping with input!)
 * @param bit number of bits per integer
 */
void pack_FOR(uint32_t offset, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t bit);
/**
 * @brief Determine the number of bits needed to represent the largest of 128 32-bit integers with frame of reference encoding
 * 
 * For frame of reference encoding, subtract a constant from inputs prior to packing
 * 
 * @param offset Value to apply as offset
 * @param in pointer to input integers
 * @return uint32_t number of bits needed per integer for bitpacking
 */
uint32_t maxbits_FOR(uint32_t offset, const uint32_t *HWY_RESTRICT in);



/**
 * @brief Unpack 128 bitpacked 32-bit integers with difference-from-reference encoding
 * 
 * For difference-from-reference encoding, subtract a vector of 128 integers prior to encoding
 * 
 * @param ref pointer to 128 values to use as offset
 * @param in pointer to input integers
 * @param out pointer to output destination (must be non-overlapping with input!)
 * @param bit number of bits per integer
 */
void unpack_diff(const uint32_t *HWY_RESTRICT ref, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t bit);
/**
 * @brief Bitpack 128 32-bit integers with difference-from-reference encoding
 * 
 * For difference-from-reference encoding, subtract a vector of 128 integers prior to encoding
 * 
 * @param ref pointer to 128 values to use as offset
 * @param in pointer to input integers
 * @param out pointer to output destination (must be non-overlapping with input!)
 * @param bit number of bits per integer
 */
void pack_diff(const uint32_t *HWY_RESTRICT ref, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t bit);
/**
 * @brief Determine the number of bits needed to represent the largest of 128 32-bit integers with difference-from-reference encoding
 * 
 * For difference-from-reference encoding, subtract a vector of 128 integers prior to encoding
 * 
 * @param ref pointer to 128 values to use as offset
 * @param in pointer to input integers
 * @return uint32_t number of bits needed per integer for bitpacking
 */
uint32_t maxbits_diff(const uint32_t *HWY_RESTRICT ref, const uint32_t *HWY_RESTRICT in);


/**
 * @brief Unpack N*128 bitpacked 32-bit integers
 * 
 * @param n Number of 128-integer chunks
 * @param in (in, length sum(bit)*4)pointer to bitpacked integers
 * @param out (out, length N*128) pointer to output destination (must be non-overlapping with input!)
 * @param bit (in, length N) pointer to number of bits per integer
 */
void unpack_Nx128(uint32_t n, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t *HWY_RESTRICT bit);
/**
 * @brief Bitpack N*128 32-bit integers
 * 
 * Calculate optimal bitwidth for each of N 128-integer chunk (stored to bit), then write bitpacked data to out.
 * 
 * @param n Number of 128-integer chunks
 * @param in (in, length N*128) pointer to input integers
 * @param out (out, capacity N*128) pointer to output destination (must be non-overlapping with input!)
 * @param bit (out, length N) pointer to number of bits per integer
 */
void pack_Nx128(uint32_t n, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, uint32_t *HWY_RESTRICT bit);


/**
 * @brief Unpack N*128 bitpacked 32-bit integers with d1 encoding
 * 
 * For d1 encoding, take the difference between adjacent inputs prior to packing
 * 
 * @param n Number of 128-integer chunks
 * @param initvalue (in, length N) Value to be subtracted from the fist input
 * @param in (in, length sum(bit)*4)pointer to bitpacked integers
 * @param out (out, length N*128) pointer to output destination (must be non-overlapping with input!)
 * @param bit (in, length N) pointer to number of bits per integer
 */
void unpack_d1_Nx128(uint32_t n, const uint32_t *HWY_RESTRICT initvalue, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t *HWY_RESTRICT bit);
/**
 * @brief Bitpack 128 32-bit integers with d1 encoding
 * 
 * Calculate optimal bitwidth for each of N 128-integer chunk (stored to bit), then write bitpacked data to out.
 * 
 * For d1 encoding, take the difference between adjacent inputs prior to packing. Uses the first value from each
 * chunk of 128 as the initvalue.
 * 
 * @param n Number of 128-integer chunks
 * @param initvalue (out, length N) Value to be subtracted from the fist input during unpacking
 * @param in (in, length N*128) pointer to input integers
 * @param out (out, capacity N*128) pointer to output destination (must be non-overlapping with input!)
 * @param bit (out, length N) pointer to number of bits per integer
 */
void pack_d1_Nx128(uint32_t n, uint32_t *HWY_RESTRICT initvalue, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, uint32_t *HWY_RESTRICT bit);


/**
 * @brief Unpack N*128 bitpacked 32-bit integers with d1z encoding
 * 
 * For d1z encoding, take the difference between adjacent inputs, then
 * zigzag encode (0,-1,1,-2,2...) -> (0,1,2,3,4,...)
 * 
 * @param n Number of 128-integer chunks
 * @param initvalue (in, length N) Value to be subtracted from the fist input
 * @param in (in, length sum(bit)*4)pointer to bitpacked integers
 * @param out (out, length N*128) pointer to output destination (must be non-overlapping with input!)
 * @param bit (in, length N) pointer to number of bits per integer
 */
void unpack_d1z_Nx128(uint32_t n, const uint32_t *HWY_RESTRICT initvalue, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t *HWY_RESTRICT bit);

/**
 * @brief Bitpack N*128 32-bit integers with d1z encoding
 * 
 * Calculate optimal bitwidth for each of N 128-integer chunk (stored to bit), then write bitpacked data to out.
 * 
 * For d1z encoding, take the difference between adjacent inputs, then
 * zigzag encode (0,-1,1,-2,2...) -> (0,1,2,3,4,...)
 * Uses the first value from each chunk of 128 as the initvalue.
 * 
 * @param n Number of 128-integer chunks
 * @param initvalue (out, length N) Value to be subtracted from the fist input
 * @param in (in, length N*128) pointer to input integers
 * @param out (out, capacity N*128) pointer to output destination (must be non-overlapping with input!)
 * @param bit (out, length N) pointer to number of bits per integer
 */
void pack_d1z_Nx128(uint32_t n, uint32_t *HWY_RESTRICT initvalue, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, uint32_t *HWY_RESTRICT bit);


/**
 * @brief Unpack N*128 bitpacked 32-bit integers with frame of reference encoding
 * 
 * For frame of reference encoding, subtract a constant from inputs prior to packing
 * 
 * @param n Number of 128-integer chunks
 * @param offset Value to apply as offset
 * @param in (in, length sum(bit)*4)pointer to bitpacked integers
 * @param out (out, length N*128) pointer to output destination (must be non-overlapping with input!)
 * @param bit (in, length N) pointer to number of bits per integer
 */
void unpack_FOR_Nx128(uint32_t n, uint32_t offset, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t *HWY_RESTRICT bit);
/**
 * @brief Bitpack N*128 32-bit integers with frame of reference encoding
 * 
 * For frame of reference encoding, subtract a constant from inputs prior to packing
 * 
 * @param n Number of 128-integer chunks
 * @param offset Value to apply as offset
 * @param in (in, length N*128) pointer to input integers
 * @param out (out, capacity N*128) pointer to output destination (must be non-overlapping with input!)
 * @param bit (out, length N) pointer to number of bits per integer
 */
void pack_FOR_Nx128(uint32_t n, uint32_t offset, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, uint32_t *HWY_RESTRICT bit);


/**
 * @brief Unpack N*128 bitpacked 32-bit integers with difference-from-reference encoding
 * 
 * For difference-from-reference encoding, subtract a vector of 128 integers prior to encoding
 * 
 * @param n Number of 128-integer chunks
 * @param ref (in, length N*128) pointer to 128 values to use as offset
 * @param in (in, length sum(bit)*4)pointer to bitpacked integers
 * @param out (out, length N*128) pointer to output destination (must be non-overlapping with input!)
 * @param bit (in, length N) pointer to number of bits per integer
 */
void unpack_diff_Nx128(uint32_t n, const uint32_t *HWY_RESTRICT ref, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t *HWY_RESTRICT bit);
/**
 * @brief Bitpack N*128 32-bit integers with difference-from-reference encoding
 * 
 * For difference-from-reference encoding, subtract a vector of 128 integers prior to encoding
 * 
 * @param n Number of 128-integer chunks
 * @param ref (in, length N*128) pointer to 128 values to use as offset
 * @param in (in, length N*128) pointer to input integers
 * @param out (out, capacity N*128) pointer to output destination (must be non-overlapping with input!)
 * @param bit (out, length N) pointer to number of bits per integer
 */
void pack_diff_Nx128(uint32_t n, const uint32_t *HWY_RESTRICT ref, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, uint32_t *HWY_RESTRICT bit);


/**
 * @brief Return the current target (SIMD feature set) used for bp128 operations
 * 
 * @return const char* 
 */
const char *current_target();

} // end namespace BPCells