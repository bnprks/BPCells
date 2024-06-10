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
 * @brief Return the current target (SIMD feature set) used for bp128 operations
 * 
 * @return const char* 
 */
const char *current_target();

} // end namespace BPCells