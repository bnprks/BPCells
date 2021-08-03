#pragma once

#include <cstdint>

// #if defined(__SSSE4__) || defined(__AVX__)
#include <smmintrin.h>
// #endif

namespace BPCells {


// Write end - start to out for 128 integers
void simdsubtract(const uint32_t *end, const uint32_t *start, uint32_t *out);

// Write end + start to out for 128 integers
void simdadd(const uint32_t *start, const uint32_t *end, uint32_t *out);

// Get the maximum value of 128 integers
uint32_t simdmax(const uint32_t *in);

//From https://github.com/lemire/simdcomp/blob/dd9317ff7df8ad6350d492e6a54600a9a69e7919/src/simdcomputil.c#L16-L29
inline uint32_t bits(const uint32_t v) {
#ifdef _MSC_VER
  unsigned long answer;
  if (v == 0) {
    return 0;
  }
  _BitScanReverse(&answer, v);
  return answer + 1;
#else
  return v == 0 ? 0
                : 32 - __builtin_clz(
                           v); /* assume GCC-like compiler if not microsoft */
#endif
}


// #if defined(__SSSE4__) || defined(__AVX__)
using vec = __m128i;
// #endif

inline vec load(vec const* mem_addr) { return _mm_loadu_si128(mem_addr); }
inline vec load_aligned(vec const* mem_addr) { return _mm_load_si128(mem_addr); }

inline void store(vec * mem_addr, vec v) { _mm_storeu_si128(mem_addr, v); }
inline void store_aligned(vec * mem_addr, vec v) { _mm_store_si128(mem_addr, v); }

// Move a vector of 4 32-bit ints left/right by 1 or 2 elements
inline vec move_l_1(const vec &v) {return _mm_slli_si128(v, 4); }
inline vec move_l_2(const vec &v) {return _mm_slli_si128(v, 8); }
inline vec move_r_1(const vec &v) {return _mm_srli_si128(v, 4); }
inline vec move_r_2(const vec &v) {return _mm_srli_si128(v, 8); }

inline uint32_t extract_bottom(const vec &v) {return _mm_cvtsi128_si32(v); }

inline vec splat(int32_t i) { return _mm_set1_epi32(i); }

inline vec shift_l(const vec &v, uint32_t i) { return _mm_slli_epi32(v, i); }
inline vec shift_r(const vec &v, uint32_t i) { return _mm_srli_epi32(v, i); }

inline vec bitwise_or(const vec &v1, const vec &v2) { return _mm_or_si128(v1, v2); }
inline vec bitwise_and(const vec &v1, const vec &v2) { return _mm_and_si128(v1, v2); }

inline vec add(const vec &v1, const vec &v2) { return _mm_add_epi32(v1, v2); }
inline vec sub(const vec &v1, const vec &v2) { return _mm_sub_epi32(v1, v2); }

// Note: these min & max instructions are the only bits newer than SSE2
inline vec min(const vec &v1, const vec &v2) { return _mm_min_epu32(v1, v2); }
inline vec max(const vec &v1, const vec &v2) { return _mm_max_epu32(v1, v2); }

// Input: [a,b,c,d]; [?,?,?,h]; Output [a-h, b-a, c-b, d-c]
inline vec delta(const vec &v, const vec& initOffset) {
    return _mm_sub_epi32(v, _mm_alignr_epi8(v, initOffset, 12));
}

// Input: [a,b,c,d]; [?,?,?,h]; Output [a+h, a+b+h, a+b+c+h, a+b+c+d+h]
inline vec prefixSum(const vec&v, const vec& initOffset) {
    const vec _tmp1 = add(move_l_2(v), v);
    const vec _tmp2 = add(move_l_1(_tmp1), _tmp1);
    return add(_tmp2, _mm_shuffle_epi32(initOffset, 0xff));
}

// Return number of bits required to represent largest number in accumulator
inline uint32_t maxbitas32int(const vec accumulator) {
  const vec _tmp1 = bitwise_or(
      move_r_2(accumulator),
      accumulator); /* (A,B,C,D) or (0,0,A,B) = (A,B,C or A,D or B)*/
  const vec _tmp2 =
      bitwise_or(move_r_1(_tmp1),
                   _tmp1); /*  (A,B,C or A,D or B) or  (0,0,0,C or A)*/
  uint32_t ans = extract_bottom(_tmp2);
  return bits(ans);
}

//Return number of bits required to represent the difference between mins and maxs
//Used for setting up FOR encoding
//bit_count and minvalue are return variables
inline void maxdiffbitsas32int(const vec mins, const vec maxs, uint32_t & bit_count, uint32_t & minvalue) {
  //last of _tmp2 holds the min of mins
  const vec _tmp1 = min(move_r_2(mins), mins); 
  const vec _tmp2 = min(move_r_1(_tmp1), _tmp1);
  
  //last of _tmp4 holds the max of maxs
  const vec _tmp3 = max(move_r_2(maxs), maxs); 
  const vec _tmp4 = max(move_r_1(_tmp3), _tmp3);

  bit_count = bits(extract_bottom(sub(_tmp4, _tmp2)));
  minvalue = extract_bottom(_tmp2);
}


} // end namespace BPCells