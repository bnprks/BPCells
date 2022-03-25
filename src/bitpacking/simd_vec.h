#pragma once

// This provides a basic set of SIMD operations on vectors of 4 x uint32_t
// for ARM NEON, x86, and non-simd fallback implementations in standard C

#include <cstdint>

#define _SIMDBP128_X86_FULL 1
#define _SIMDBP128_X86 2
#define _SIMDBP128_ARM_NEON 3
#define _SIMDBP128_C_FALLBACK 0

#if defined(__SSE4_1__) || defined(__AVX__)
  // x86 with full SSE4.1 support
  #include <smmintrin.h>
  #define _SIMDBP128_MODE_ _SIMDBP128_X86_FULL
#elif defined(__SSE2__)
  // x86 with SSE2 but no SSE4.1
  #include <emmintrin.h>
  #define _SIMDBP128_MODE_ _SIMDBP128_X86
#elif defined(__ARM_NEON)
  // ARM NEON
  #include <arm_neon.h>
  #define _SIMDBP128_MODE_ _SIMDBP128_ARM_NEON
#else
  // Fallback portable C (non-simd)
  #define _SIMDBP128_MODE_ _SIMDBP128_C_FALLBACK
#endif

namespace BPCells {


// Write end - start to end for 128 integers
void simdsubtract(uint32_t *end, const uint32_t *start);

// Write end + start to end for 128 integers
void simdadd(uint32_t *end, const uint32_t *start);

// Get the maximum value of 128 integers
uint32_t simdmax(const uint32_t *in);

//From https://github.com/lemire/simdcomp/blob/dd9317ff7df8ad6350d492e6a54600a9a69e7919/src/simdcomputil.c#L16-L29
// Return number of bits needed to represent v
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


#if _SIMDBP128_MODE_ == _SIMDBP128_X86 || _SIMDBP128_MODE_ == _SIMDBP128_X86_FULL
using vec = __m128i;
#elif _SIMDBP128_MODE_ == _SIMDBP128_ARM_NEON
using vec = uint32x4_t;
#else
using vec = struct {
    uint32_t x0, x1, x2, x3;
};
#endif

#if _SIMDBP128_MODE_ == _SIMDBP128_X86 || _SIMDBP128_MODE_ == _SIMDBP128_X86_FULL
// Intrinsics Guide: https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html 
inline vec load(vec const* mem_addr) { return _mm_loadu_si128(mem_addr); }
inline vec load_aligned(vec const* mem_addr) { return _mm_load_si128(mem_addr); }

inline void store(vec * mem_addr, vec v) { _mm_storeu_si128(mem_addr, v); }
inline void store_aligned(vec * mem_addr, vec v) { _mm_store_si128(mem_addr, v); }

inline vec splat(uint32_t i) { return _mm_set1_epi32(i); }

// Move a vector of 4 32-bit ints left/right by 1 or 2 elements
inline vec move_r_1(const vec &v) {return _mm_slli_si128(v, 4); }
inline vec move_r_2(const vec &v) {return _mm_slli_si128(v, 8); }
inline vec move_r_3(const vec &v) {return _mm_slli_si128(v, 12); }
inline vec move_l_1(const vec &v) {return _mm_srli_si128(v, 4); }
inline vec move_l_2(const vec &v) {return _mm_srli_si128(v, 8); }
inline vec move_l_3(const vec &v) {return _mm_srli_si128(v, 12); }

inline uint32_t extract_left(const vec &v) {return _mm_cvtsi128_si32(v); }

inline vec shift_l(const vec &v, uint32_t i) { return _mm_slli_epi32(v, i); }
inline vec shift_r(const vec &v, uint32_t i) { return _mm_srli_epi32(v, i); }
inline vec shift_r_arith(const vec &v, uint32_t i) { return _mm_srai_epi32(v, i); }

inline vec bitwise_or(const vec &v1, const vec &v2) { return _mm_or_si128(v1, v2); }
inline vec bitwise_and(const vec &v1, const vec &v2) { return _mm_and_si128(v1, v2); }
inline vec bitwise_xor(const vec &v1, const vec &v2) { return _mm_xor_si128(v1, v2); }

inline vec add(const vec &v1, const vec &v2) { return _mm_add_epi32(v1, v2); }
inline vec sub(const vec &v1, const vec &v2) { return _mm_sub_epi32(v1, v2); }

// Signed comparisons -- if true all 1 bits, if false all 0 bits
inline vec cmp_lt_signed(const vec &v1, const vec &v2) {return _mm_cmplt_epi32(v1, v2); }
inline vec cmp_gt_signed(const vec &v1, const vec &v2) {return _mm_cmpgt_epi32(v1, v2); }

// Input: [a,b,c,d]; [?,?,?,h]; Output [a+h, a+b+h, a+b+c+h, a+b+c+d+h]
inline vec prefixSum(const vec&v, const vec& initOffset) {
    const vec _tmp1 = add(move_r_2(v), v);
    const vec _tmp2 = add(move_r_1(_tmp1), _tmp1);
    return add(_tmp2, _mm_shuffle_epi32(initOffset, 0xff));
}

#endif

#if _SIMDBP128_MODE_ == _SIMDBP128_X86_FULL
// Note: these min & max instructions are the only bits newer than SSE2
inline vec min(const vec &v1, const vec &v2) { return _mm_min_epu32(v1, v2); }
inline vec max(const vec &v1, const vec &v2) { return _mm_max_epu32(v1, v2); }

// Input: [a,b,c,d]; [?,?,?,h]; Output [a-h, b-a, c-b, d-c]
inline vec delta(const vec &v, const vec& initOffset) {
    return sub(v, _mm_alignr_epi8(v, initOffset, 12));
}
#elif _SIMDBP128_MODE_ == _SIMDBP128_X86
inline vec min(const vec &v1, const vec &v2) { 
  uint32_t a1[4], a2[4];
  store((vec *) a1, v1);
  store((vec *) a2, v2);
  a1[0] = a1[0] > a2[0] ? a2[0] : a1[0];
  a1[1] = a1[1] > a2[1] ? a2[1] : a1[1];
  a1[2] = a1[2] > a2[2] ? a2[2] : a1[2];
  a1[3] = a1[3] > a2[3] ? a2[3] : a1[3];
  return load((vec *) a1); 
}
inline vec max(const vec &v1, const vec &v2) { 
  uint32_t a1[4], a2[4];
  store((vec *) a1, v1);
  store((vec *) a2, v2);
  a1[0] = a1[0] < a2[0] ? a2[0] : a1[0];
  a1[1] = a1[1] < a2[1] ? a2[1] : a1[1];
  a1[2] = a1[2] < a2[2] ? a2[2] : a1[2];
  a1[3] = a1[3] < a2[3] ? a2[3] : a1[3];
  return load((vec *) a1); 
}

// Input: [a,b,c,d]; [?,?,?,h]; Output [a-h, b-a, c-b, d-c]
inline vec delta(const vec &v, const vec& initOffset) {
  return sub(v, bitwise_or(move_r_1(v), move_l_3(initOffset)));
}
#endif

#if _SIMDBP128_MODE_ == _SIMDBP128_ARM_NEON
// Intrinsics Guide: https://developer.arm.com/architectures/instruction-sets/intrinsics/
inline vec load(vec const* mem_addr) { return vld1q_u32((uint32_t *) mem_addr); }
inline vec load_aligned(vec const* mem_addr) { return vld1q_u32((uint32_t *) mem_addr); }

inline void store(vec * mem_addr, vec v) { vst1q_u32((uint32_t *) mem_addr, v); }
inline void store_aligned(vec * mem_addr, vec v) { vst1q_u32((uint32_t *) mem_addr, v); }

inline vec splat(uint32_t i) { return vdupq_n_u32(i); }

// Move a vector of 4 32-bit ints left/right by 1 or 2 elements
inline vec move_r_1(const vec &v) {return vextq_u32(splat(0), v, 3); }
inline vec move_r_2(const vec &v) {return vextq_u32(splat(0), v, 2); }
inline vec move_r_3(const vec &v) {return vextq_u32(splat(0), v, 1); }
inline vec move_l_1(const vec &v) {return vextq_u32(v, splat(0), 1); }
inline vec move_l_2(const vec &v) {return vextq_u32(v, splat(0), 2); }
inline vec move_l_3(const vec &v) {return vextq_u32(v, splat(0), 3); }

inline uint32_t extract_left(const vec &v) {return vgetq_lane_u32(v, 0); }


inline vec shift_l(const vec &v, const uint32_t i) { return vshlq_u32(v, vdupq_n_s32(i)); }
inline vec shift_r(const vec &v, const uint32_t i) { return vshlq_u32(v, vdupq_n_s32(-((int32_t) i))); }
// Arithmetic shift that shifts in the sign bit
inline vec shift_r_arith(const vec &v, const uint32_t i) { return (vec) vshlq_s32((int32x4_t) v, vdupq_n_s32(-((int32_t) i))); }

inline vec bitwise_or(const vec &v1, const vec &v2) { return vorrq_u32(v1, v2); }
inline vec bitwise_and(const vec &v1, const vec &v2) { return vandq_u32(v1, v2); }
inline vec bitwise_xor(const vec &v1, const vec &v2) { return veorq_u32(v1, v2); }

inline vec add(const vec &v1, const vec &v2) { return vaddq_u32(v1, v2); }
inline vec sub(const vec &v1, const vec &v2) { return vsubq_u32(v1, v2); }

// Signed comparisons -- if true all 1 bits, if false all 0 bits
inline vec cmp_lt_signed(const vec &v1, const vec &v2) {return vcltq_s32((int32x4_t) v1, (int32x4_t) v2); }
inline vec cmp_gt_signed(const vec &v1, const vec &v2) {return vcgtq_s32((int32x4_t) v1, (int32x4_t) v2); }

// Input: [a,b,c,d]; [?,?,?,h]; Output [a+h, a+b+h, a+b+c+h, a+b+c+d+h]
inline vec prefixSum(const vec&v, const vec& initOffset) {
    const vec _tmp1 = add(move_r_2(v), v);
    const vec _tmp2 = add(move_r_1(_tmp1), _tmp1);
    return add(_tmp2, vdupq_laneq_u32(initOffset, 3));
}

// Note: these min & max instructions are the only bits newer than SSE2
inline vec min(const vec &v1, const vec &v2) { return vminq_u32(v1, v2); }
inline vec max(const vec &v1, const vec &v2) { return vmaxq_u32(v1, v2); }

// Input: [a,b,c,d]; [?,?,?,h]; Output [a-h, b-a, c-b, d-c]
inline vec delta(const vec &v, const vec& initOffset) {
    return sub(v, vextq_u32(initOffset, v, 3));
}
#endif


#if _SIMDBP128_MODE_ == _SIMDBP128_C_FALLBACK

inline vec load(vec const* mem_addr) { uint32_t* m = (uint32_t*) mem_addr; return {m[0], m[1], m[2], m[3]}; }
inline vec load_aligned(vec const* mem_addr) { return load(mem_addr); }

inline void store(vec * mem_addr, vec v) { uint32_t* m = (uint32_t*) mem_addr; m[0] = v.x0; m[1] = v.x1; m[2] = v.x2; m[3] = v.x3;}
inline void store_aligned(vec * mem_addr, vec v) { store(mem_addr, v); }

inline vec splat(uint32_t i) { return {i, i, i, i}; }

// Move a vector of 4 32-bit ints left/right by 1 or 2 elements
inline vec move_r_1(const vec &v) {return {0, v.x0, v.x1, v.x2}; }
inline vec move_r_2(const vec &v) {return {0, 0, v.x0, v.x1}; }
inline vec move_r_3(const vec &v) {return {0, 0, 0, v.x0}; }
inline vec move_l_1(const vec &v) {return {v.x1, v.x2, v.x3, 0}; }
inline vec move_l_2(const vec &v) {return {v.x2, v.x3, 0, 0}; }
inline vec move_l_3(const vec &v) {return {v.x3, 0, 0 , 0}; }

inline uint32_t extract_left(const vec &v) {return v.x0; }


inline vec shift_l(const vec &v, uint32_t i) { return {v.x0 << i, v.x1 << i, v.x2 << i, v.x3 << i}; }
inline vec shift_r(const vec &v, uint32_t i) { return {v.x0 >> i, v.x1 >> i, v.x2 >> i, v.x3 >> i}; }
// Arithmetic shift that shifts in the sign bit
inline vec shift_r_arith(const vec &v, uint32_t i) { 
  // GCC, clang, and MSVC all implement signed shifts as arithmetic, 
  // although technically this is implementation-defined in the C spec
  return {
    (uint32_t) ((int32_t) v.x0 >> i), 
    (uint32_t) ((int32_t) v.x1 >> i), 
    (uint32_t) ((int32_t) v.x2 >> i), 
    (uint32_t) ((int32_t) v.x3 >> i)};
}

inline vec bitwise_or(const vec &v1, const vec &v2) { return {v1.x0 | v2.x0, v1.x1 | v2.x1, v1.x2 | v2.x2, v1.x3 | v2.x3}; }
inline vec bitwise_and(const vec &v1, const vec &v2) { return {v1.x0 & v2.x0, v1.x1 & v2.x1, v1.x2 & v2.x2, v1.x3 & v2.x3}; }
inline vec bitwise_xor(const vec &v1, const vec &v2) { return {v1.x0 ^ v2.x0, v1.x1 ^ v2.x1, v1.x2 ^ v2.x2, v1.x3 ^ v2.x3}; }

inline vec add(const vec &v1, const vec &v2) { return {v1.x0 + v2.x0, v1.x1 + v2.x1, v1.x2 + v2.x2, v1.x3 + v2.x3}; }
inline vec sub(const vec &v1, const vec &v2) { return {v1.x0 - v2.x0, v1.x1 - v2.x1, v1.x2 - v2.x2, v1.x3 - v2.x3}; }

// Signed comparisons -- if true all 1 bits, if false all 0 bits
inline vec cmp_lt_signed(const vec &v1, const vec &v2) {
  return {
    (int32_t) v1.x0 < (int32_t) v2.x0 ? 0xFFFFFFFF : 0,
    (int32_t) v1.x1 < (int32_t) v2.x1 ? 0xFFFFFFFF : 0,
    (int32_t) v1.x2 < (int32_t) v2.x2 ? 0xFFFFFFFF : 0,
    (int32_t) v1.x3 < (int32_t) v2.x3 ? 0xFFFFFFFF : 0
  }; 
}
inline vec cmp_gt_signed(const vec &v1, const vec &v2) {
  return {
    (int32_t) v1.x0 > (int32_t) v2.x0 ? 0xFFFFFFFF : 0,
    (int32_t) v1.x1 > (int32_t) v2.x1 ? 0xFFFFFFFF : 0,
    (int32_t) v1.x2 > (int32_t) v2.x2 ? 0xFFFFFFFF : 0,
    (int32_t) v1.x3 > (int32_t) v2.x3 ? 0xFFFFFFFF : 0
  }; 
}

// Input: [a,b,c,d]; [?,?,?,h]; Output [a+h, a+b+h, a+b+c+h, a+b+c+d+h]
inline vec prefixSum(const vec&v, const vec& initOffset) {
    return {
      v.x0 + initOffset.x3, 
      v.x0 + v.x1 + initOffset.x3,
      v.x0 + v.x1 + v.x2 + initOffset.x3, 
      v.x0 + v.x1 + v.x2 + v.x3 + initOffset.x3
    };
}

// Note: these min & max instructions are the only bits newer than SSE2
inline vec min(const vec &v1, const vec &v2) { 
  return {
      v1.x0 < v2.x0 ? v1.x0 : v2.x0, 
      v1.x1 < v2.x1 ? v1.x1 : v2.x1, 
      v1.x2 < v2.x2 ? v1.x2 : v2.x2, 
      v1.x3 < v2.x3 ? v1.x3 : v2.x3
    };
 }
inline vec max(const vec &v1, const vec &v2) { 
  return {
      v1.x0 > v2.x0 ? v1.x0 : v2.x0, 
      v1.x1 > v2.x1 ? v1.x1 : v2.x1, 
      v1.x2 > v2.x2 ? v1.x2 : v2.x2, 
      v1.x3 > v2.x3 ? v1.x3 : v2.x3
    };
 }

// Input: [a,b,c,d]; [?,?,?,h]; Output [a-h, b-a, c-b, d-c]
inline vec delta(const vec &v, const vec& initOffset) {
    return {
      v.x0 - initOffset.x3, 
      v.x1 - v.x0,
      v.x2 - v.x1,
      v.x3 - v.x2
    };
}
#endif



// Return number of bits required to represent largest number in accumulator
inline uint32_t maxbitas32int(const vec accumulator) {
  const vec _tmp1 = bitwise_or(
      move_l_2(accumulator),
      accumulator); /* (A,B,C,D) or (0,0,A,B) = (A,B,C or A,D or B)*/
  const vec _tmp2 =
      bitwise_or(move_l_1(_tmp1),
                   _tmp1); /*  (A,B,C or A,D or B) or  (0,0,0,C or A)*/
  uint32_t ans = extract_left(_tmp2);
  return bits(ans);
}

//Return number of bits required to represent the difference between mins and maxs
//Used for setting up FOR encoding
//bit_count and minvalue are return variables
inline void maxdiffbitsas32int(const vec mins, const vec maxs, uint32_t & bit_count, uint32_t & minvalue) {
  //last of _tmp2 holds the min of mins
  const vec _tmp1 = min(move_l_2(mins), mins); 
  const vec _tmp2 = min(move_l_1(_tmp1), _tmp1);
  
  //last of _tmp4 holds the max of maxs
  const vec _tmp3 = max(move_l_2(maxs), maxs); 
  const vec _tmp4 = max(move_l_1(_tmp3), _tmp3);

  bit_count = bits(extract_left(sub(_tmp4, _tmp2)));
  minvalue = extract_left(_tmp2);
}


} // end namespace BPCells