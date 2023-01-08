#pragma once

// Wrapper to handle providing sleef math functions in a uniform way to downstream users
#include <cmath>
#include <cstdint>


#define _BPCELLS_SLEEF_FALLBACK 0
#define _BPCELLS_SLEEF_SSE2 1
#define _BPCELLS_SLEEF_AVX 2
#define _BPCELLS_SLEEF_AVX2 3
#define _BPCELLS_SLEEF_NEON 4

#if defined(__GNU_C__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wpedantic"
#elif defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#pragma clang diagnostic ignored "-Wc99-extensions"
#endif

// Decide what architecture to use
#ifdef BPCELLS_SLEEF_FALLBACK
// Allow a forced fallback to C for testing purposes
#define BPCELLS_SLEEF_MODE _BPCELLS_SLEEF_FALLBACK
#elif defined(__AVX2__)

#include <cstring>
#include <x86intrin.h>

#include "sleef/sleefinline_avx2.h"
#define BPCELLS_SLEEF_MODE _BPCELLS_SLEEF_AVX2

#elif defined(__AVX__)

#include <cstring>
#include <x86intrin.h>

#include "sleef/sleefinline_avx.h"
#define BPCELLS_SLEEF_MODE _BPCELLS_SLEEF_AVX

#elif defined(__SSE2__)

#include <cstring>
#include <x86intrin.h>

#include "sleef/sleefinline_sse2.h"
#define BPCELLS_SLEEF_MODE _BPCELLS_SLEEF_SSE2

#elif defined(__ARM_NEON)

#include <arm_neon.h>
#include <cstring>

#include "sleef/sleefinline_advsimd.h"
#define BPCELLS_SLEEF_MODE _BPCELLS_SLEEF_NEON

#else

#define BPCELLS_SLEEF_MODE _BPCELLS_SLEEF_FALLBACK

#endif

#ifdef __clang__
#pragma clang diagnostic pop
#elif defined(__GNU_C__)
#pragma GCC diagnostic pop
#endif

namespace BPCells {

// Note: may need to improve
#if BPCELLS_SLEEF_MODE == _BPCELLS_SLEEF_AVX2

#define BPCELLS_VEC_FLOAT_SIZE 8
#define BPCELLS_VEC_DOUBLE_SIZE 4

using vec_float = __m256;
using vec_double = __m256d;

inline vec_float load_float(float const *mem_addr) { return _mm256_loadu_ps(mem_addr); }
inline void store_float(float *mem_addr, vec_float v) { _mm256_storeu_ps(mem_addr, v); }
inline vec_float splat_float(float f) { return _mm256_set1_ps(f); }

inline vec_double load_double(double const *mem_addr) { return _mm256_loadu_pd(mem_addr); }
inline void store_double(double *mem_addr, vec_double v) { _mm256_storeu_pd(mem_addr, v); }
inline vec_double splat_double(double f) { return _mm256_set1_pd(f); }

inline vec_float log1p_f(const vec_float &v) { return Sleef_log1pf8_u10avx2(v); }
inline vec_double log1p_d(const vec_double &v) { return Sleef_log1pd4_u10avx2(v); }

inline vec_float expm1_f(const vec_float &v) { return Sleef_expm1f8_u10avx2(v); }
inline vec_double expm1_d(const vec_double &v) { return Sleef_expm1d4_u10avx2(v); }

inline vec_float pow_f(const vec_float &v, const vec_float &exp) {
    return Sleef_powf8_u10avx2(v, exp);
}
inline vec_double pow_d(const vec_double &v, const vec_double &exp) {
    return Sleef_powd4_u10avx2(v, exp);
}

#elif BPCELLS_SLEEF_MODE == _BPCELLS_SLEEF_AVX

#define BPCELLS_VEC_FLOAT_SIZE 8
#define BPCELLS_VEC_DOUBLE_SIZE 4

using vec_float = __m256;
using vec_double = __m256d;

inline vec_float load_float(float const *mem_addr) { return _mm256_loadu_ps(mem_addr); }
inline void store_float(float *mem_addr, vec_float v) { _mm256_storeu_ps(mem_addr, v); }
inline vec_float splat_float(float f) { return _mm256_set1_ps(f); }

inline vec_double load_double(double const *mem_addr) { return _mm256_loadu_pd(mem_addr); }
inline void store_double(double *mem_addr, vec_double v) { _mm256_storeu_pd(mem_addr, v); }
inline vec_double splat_double(double f) { return _mm256_set1_pd(f); }

inline vec_float log1p_f(const vec_float &v) { return Sleef_log1pf8_u10avx(v); }
inline vec_double log1p_d(const vec_double &v) { return Sleef_log1pd4_u10avx(v); }

inline vec_float expm1_f(const vec_float &v) { return Sleef_expm1f8_u10avx(v); }
inline vec_double expm1_d(const vec_double &v) { return Sleef_expm1d4_u10avx(v); }

inline vec_float pow_f(const vec_float &v, const vec_float &exp) {
    return Sleef_powf8_u10avx(v, exp);
}
inline vec_double pow_d(const vec_double &v, const vec_double &exp) {
    return Sleef_powd4_u10avx(v, exp);
}

#elif BPCELLS_SLEEF_MODE == _BPCELLS_SLEEF_SSE2

#define BPCELLS_VEC_FLOAT_SIZE 4
#define BPCELLS_VEC_DOUBLE_SIZE 2

using vec_float = __m128;
using vec_double = __m128d;

inline vec_float load_float(float const *mem_addr) { return _mm_loadu_ps(mem_addr); }
inline void store_float(float *mem_addr, vec_float v) { _mm_storeu_ps(mem_addr, v); }
inline vec_float splat_float(float f) { return _mm_set1_ps(f); }

inline vec_double load_double(double const *mem_addr) { return _mm_loadu_pd(mem_addr); }
inline void store_double(double *mem_addr, vec_double v) { _mm_storeu_pd(mem_addr, v); }
inline vec_double splat_double(double f) { return _mm_set1_pd(f); }

inline vec_float log1p_f(const vec_float &v) { return Sleef_log1pf4_u10sse2(v); }
inline vec_double log1p_d(const vec_double &v) { return Sleef_log1pd2_u10sse2(v); }

inline vec_float expm1_f(const vec_float &v) { return Sleef_expm1f4_u10sse2(v); }
inline vec_double expm1_d(const vec_double &v) { return Sleef_expm1d2_u10sse2(v); }

inline vec_float pow_f(const vec_float &v, const vec_float &exp) {
    return Sleef_powf4_u10sse2(v, exp);
}
inline vec_double pow_d(const vec_double &v, const vec_double &exp) {
    return Sleef_powd2_u10sse2(v, exp);
}

#elif BPCELLS_SLEEF_MODE == _BPCELLS_SLEEF_NEON

#define BPCELLS_VEC_FLOAT_SIZE 4
#define BPCELLS_VEC_DOUBLE_SIZE 2

using vec_float = float32x4_t;
using vec_double = float64x2_t;

inline vec_float load_float(float const *mem_addr) { return vld1q_f32(mem_addr); }
inline void store_float(float *mem_addr, vec_float v) { vst1q_f32(mem_addr, v); }
inline vec_float splat_float(float f) { return vdupq_n_f32(f); }

inline vec_double load_double(double const *mem_addr) { return vld1q_f64(mem_addr); }
inline void store_double(double *mem_addr, vec_double v) { vst1q_f64(mem_addr, v); }
inline vec_double splat_double(double f) { return vdupq_n_f64(f); }

inline vec_float log1p_f(const vec_float &v) { return Sleef_log1pf4_u10advsimd(v); }
inline vec_double log1p_d(const vec_double &v) { return Sleef_log1pd2_u10advsimd(v); }

inline vec_float expm1_f(const vec_float &v) { return Sleef_expm1f4_u10advsimd(v); }
inline vec_double expm1_d(const vec_double &v) { return Sleef_expm1d2_u10advsimd(v); }

inline vec_float pow_f(const vec_float &v, const vec_float &exp) {
    return Sleef_powf4_u10advsimd(v, exp);
}
inline vec_double pow_d(const vec_double &v, const vec_double &exp) {
    return Sleef_powd2_u10advsimd(v, exp);
}

#else

#define BPCELLS_VEC_FLOAT_SIZE 4
#define BPCELLS_VEC_DOUBLE_SIZE 4

using vec_float = struct {
    float x0, x1, x2, x3;
};
using vec_double = struct {
    double x0, x1, x2, x3;
};

inline vec_float load_float(float const *mem_addr) {
    return {mem_addr[0], mem_addr[1], mem_addr[2], mem_addr[3]};
}
inline void store_float(float *mem_addr, vec_float v) {
    mem_addr[0] = v.x0;
    mem_addr[1] = v.x1;
    mem_addr[2] = v.x2;
    mem_addr[3] = v.x3;
}
inline vec_float splat_float(float f) { return {f, f, f, f}; }

inline vec_double load_double(double const *mem_addr) {
    return {mem_addr[0], mem_addr[1], mem_addr[2], mem_addr[3]};
}
inline void store_double(double *mem_addr, vec_double v) {
    mem_addr[0] = v.x0;
    mem_addr[1] = v.x1;
    mem_addr[2] = v.x2;
    mem_addr[3] = v.x3;
}
inline vec_double splat_double(double f) { return {f, f, f, f}; }

inline vec_float log1p_f(const vec_float &v) {
    return {log1pf(v.x0), log1pf(v.x1), log1pf(v.x2), log1pf(v.x3)};
}
inline vec_double log1p_d(const vec_double &v) {
    return {log1p(v.x0), log1p(v.x1), log1p(v.x2), log1p(v.x3)};
}

inline vec_float expm1_f(const vec_float &v) {
    return {expm1f(v.x0), expm1f(v.x1), expm1f(v.x2), expm1f(v.x3)};
}
inline vec_double expm1_d(const vec_double &v) {
    return {expm1(v.x0), expm1(v.x1), expm1(v.x2), expm1(v.x3)};
}

inline vec_float pow_f(const vec_float &v, const vec_float &exp) {
    return {powf(v.x0, exp.x0), powf(v.x1, exp.x1), powf(v.x2, exp.x2), powf(v.x3, exp.x3)};
}
inline vec_double pow_d(const vec_double &v, const vec_double &exp) {
    return {pow(v.x0, exp.x0), pow(v.x1, exp.x1), pow(v.x2, exp.x2), pow(v.x3, exp.x3)};
}

#endif



} // end namespace BPCells