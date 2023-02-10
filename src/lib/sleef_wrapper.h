#pragma once

// Wrapper to handle providing sleef math functions in a uniform way to downstream users
#include <cmath>
#include <cstdint>

#define _BPCELLS_SLEEF_FALLBACK 0
#define _BPCELLS_SLEEF_SSE2 1
#define _BPCELLS_SLEEF_AVX 2
#define _BPCELLS_SLEEF_AVX2 3
#define _BPCELLS_SLEEF_NEON 4

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#pragma clang diagnostic ignored "-Wpedantic"
#pragma clang diagnostic ignored "-Wc99-extensions"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
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

#include <algorithm>
#define BPCELLS_SLEEF_MODE _BPCELLS_SLEEF_FALLBACK

#endif

#ifdef __clang__
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

namespace BPCells {

#if BPCELLS_SLEEF_MODE == _BPCELLS_SLEEF_AVX2

#define BPCELLS_VEC_FLOAT_SIZE 8
#define BPCELLS_VEC_DOUBLE_SIZE 4

using vec_float = __m256;
using vec_double = __m256d;

inline vec_float load_float(float const *mem_addr) { return _mm256_loadu_ps(mem_addr); }
inline void store_float(float *mem_addr, const vec_float &v) { _mm256_storeu_ps(mem_addr, v); }
inline vec_float splat_float(float f) { return _mm256_set1_ps(f); }

inline vec_double load_double(double const *mem_addr) { return _mm256_loadu_pd(mem_addr); }
inline void store_double(double *mem_addr, const vec_double &v) { _mm256_storeu_pd(mem_addr, v); }
inline vec_double splat_double(double f) { return _mm256_set1_pd(f); }

inline vec_float load_double_to_float(double *mem_addr) {
    // See: https://stackoverflow.com/a/11117004
    __m128 v1 = _mm256_cvtpd_ps(load_double(mem_addr));
    __m128 v2 = _mm256_cvtpd_ps(load_double(mem_addr + 4));
    return _mm256_insertf128_ps(_mm256_castps128_ps256(v1), v2, 1);
}

inline void store_float_to_double(double *mem_addr, const vec_float &v) {
    // Write low floats
    store_double(mem_addr, _mm256_cvtps_pd(_mm256_castps256_ps128(v)));
    // Write high floats
    store_double(mem_addr + 4, _mm256_cvtps_pd(_mm256_extractf128_ps(v, 1)));
}

inline void store_double_to_float(float *mem_addr, const vec_double &v) {
    _mm_storeu_ps(mem_addr, _mm256_cvtpd_ps(v));
}

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

// Use dedicated FMA instruction. Returns (a*b) + c
inline vec_float fma_f(const vec_float &a, const vec_float &b, const vec_float &c) {
    return Sleef_fmaf8_avx2(a, b, c);
}
inline vec_double fma_d(const vec_double &a, const vec_double &b, const vec_double &c) {
    return Sleef_fmad4_avx2(a, b, c);
}

inline vec_float add_f(const vec_float &a, const vec_float &b) { return _mm256_add_ps(a, b); }
inline vec_double add_d(const vec_double &a, const vec_double &b) { return _mm256_add_pd(a, b); }

inline vec_float sub_f(const vec_float &a, const vec_float &b) { return _mm256_sub_ps(a, b); }
inline vec_double sub_d(const vec_double &a, const vec_double &b) { return _mm256_sub_pd(a, b); }

inline vec_float mul_f(const vec_float &a, const vec_float &b) { return _mm256_mul_ps(a, b); }
inline vec_double mul_d(const vec_double &a, const vec_double &b) { return _mm256_mul_pd(a, b); }

inline vec_float min_f(const vec_float &a, const vec_float &b) { return _mm256_min_ps(a, b); }
inline vec_double min_d(const vec_double &a, const vec_double &b) { return _mm256_min_pd(a, b); }

inline vec_float max_f(const vec_float &a, const vec_float &b) { return _mm256_max_ps(a, b); }
inline vec_double max_d(const vec_double &a, const vec_double &b) { return _mm256_max_pd(a, b); }

inline vec_float neg_f(const vec_float &a) {
    // Constant from Eigen "Core/arch/AVX/PacketMath.h"
    vec_float mask = _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000));
    return _mm256_xor_ps(a, mask);
}
inline vec_double neg_d(const vec_double &a) {
    // Constant from Eigen "Core/arch/AVX/PacketMath.h"
    vec_double mask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x8000000000000000ULL));
    return _mm256_xor_pd(a, mask);
}

template <int It> inline vec_float rsqrt_newton_f(const vec_float &a, const vec_float &approx);
template <int It> inline vec_double rsqrt_newton_d(const vec_double &a, const vec_double &approx);

inline vec_float rsqrt_f(const vec_float &a) {
    // Newton iter count from "Core/arch/AVX/MathFunctions.h"
    return rsqrt_newton_f<1>(a, _mm256_rsqrt_ps(a));
}
inline vec_double rsqrt_d(const vec_double &a) {
    // I don't think the newton iteration is worth it given the lack of fast approximation in double
    // precision
    return _mm256_div_pd(splat_double(1.0), _mm256_sqrt_pd(a));
}

#elif BPCELLS_SLEEF_MODE == _BPCELLS_SLEEF_AVX

#define BPCELLS_VEC_FLOAT_SIZE 8
#define BPCELLS_VEC_DOUBLE_SIZE 4

using vec_float = __m256;
using vec_double = __m256d;

inline vec_float load_float(float const *mem_addr) { return _mm256_loadu_ps(mem_addr); }
inline void store_float(float *mem_addr, const vec_float &v) { _mm256_storeu_ps(mem_addr, v); }
inline vec_float splat_float(float f) { return _mm256_set1_ps(f); }

inline vec_double load_double(double const *mem_addr) { return _mm256_loadu_pd(mem_addr); }
inline void store_double(double *mem_addr, const vec_double &v) { _mm256_storeu_pd(mem_addr, v); }
inline vec_double splat_double(double f) { return _mm256_set1_pd(f); }

inline vec_float load_double_to_float(double *mem_addr) {
    // See: https://stackoverflow.com/a/11117004
    __m128 v1 = _mm256_cvtpd_ps(load_double(mem_addr));
    __m128 v2 = _mm256_cvtpd_ps(load_double(mem_addr + 4));
    return _mm256_insertf128_ps(_mm256_castps128_ps256(v1), v2, 1);
}

inline void store_float_to_double(double *mem_addr, const vec_float &v) {
    // See: https://stackoverflow.com/a/66537016
    // Write low floats
    store_double(mem_addr, _mm256_cvtps_pd(_mm256_castps256_ps128(v)));
    // Write high floats
    store_double(mem_addr + 4, _mm256_cvtps_pd(_mm256_extractf128_ps(v, 1)));
}

inline void store_double_to_float(float *mem_addr, const vec_double &v) {
    _mm_storeu_ps(mem_addr, _mm256_cvtpd_ps(v));
}

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

// Use dedicated FMA instruction. Returns (a*b) + c
inline vec_float fma_f(const vec_float &a, const vec_float &b, const vec_float &c) {
    return Sleef_fmaf8_avx(a, b, c);
}
inline vec_double fma_d(const vec_double &a, const vec_double &b, const vec_double &c) {
    return Sleef_fmad4_avx(a, b, c);
}

inline vec_float add_f(const vec_float &a, const vec_float &b) { return _mm256_add_ps(a, b); }
inline vec_double add_d(const vec_double &a, const vec_double &b) { return _mm256_add_pd(a, b); }

inline vec_float sub_f(const vec_float &a, const vec_float &b) { return _mm256_sub_ps(a, b); }
inline vec_double sub_d(const vec_double &a, const vec_double &b) { return _mm256_sub_pd(a, b); }

inline vec_float mul_f(const vec_float &a, const vec_float &b) { return _mm256_mul_ps(a, b); }
inline vec_double mul_d(const vec_double &a, const vec_double &b) { return _mm256_mul_pd(a, b); }

inline vec_float min_f(const vec_float &a, const vec_float &b) { return _mm256_min_ps(a, b); }
inline vec_double min_d(const vec_double &a, const vec_double &b) { return _mm256_min_pd(a, b); }

inline vec_float max_f(const vec_float &a, const vec_float &b) { return _mm256_max_ps(a, b); }
inline vec_double max_d(const vec_double &a, const vec_double &b) { return _mm256_max_pd(a, b); }

inline vec_float neg_f(const vec_float &a) {
    // Constant from Eigen "Core/arch/AVX/PacketMath.h"
    vec_float mask = _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000));
    return _mm256_xor_ps(a, mask);
}
inline vec_double neg_d(const vec_double &a) {
    // Constant from Eigen "Core/arch/AVX/PacketMath.h"
    vec_double mask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x8000000000000000ULL));
    return _mm256_xor_pd(a, mask);
}

template <int It> inline vec_float rsqrt_newton_f(const vec_float &a, const vec_float &approx);
template <int It> inline vec_double rsqrt_newton_d(const vec_double &a, const vec_double &approx);

inline vec_float rsqrt_f(const vec_float &a) {
    // Newton iter count from "Core/arch/AVX/MathFunctions.h"
    return rsqrt_newton_f<1>(a, _mm256_rsqrt_ps(a));
}
inline vec_double rsqrt_d(const vec_double &a) {
    // I don't think the newton iteration is worth it given the lack of fast approximation in double
    // precision
    return _mm256_div_pd(splat_double(1.0), _mm256_sqrt_pd(a));
}

#elif BPCELLS_SLEEF_MODE == _BPCELLS_SLEEF_SSE2

#define BPCELLS_VEC_FLOAT_SIZE 4
#define BPCELLS_VEC_DOUBLE_SIZE 2

using vec_float = __m128;
using vec_double = __m128d;

inline vec_float load_float(float const *mem_addr) { return _mm_loadu_ps(mem_addr); }
inline void store_float(float *mem_addr, const vec_float &v) { _mm_storeu_ps(mem_addr, v); }
inline vec_float splat_float(float f) { return _mm_set1_ps(f); }

inline vec_double load_double(double const *mem_addr) { return _mm_loadu_pd(mem_addr); }
inline void store_double(double *mem_addr, const vec_double &v) { _mm_storeu_pd(mem_addr, v); }
inline vec_double splat_double(double f) { return _mm_set1_pd(f); }

inline vec_float load_double_to_float(double *mem_addr) {
    __m128 v1 = _mm_cvtpd_ps(load_double(mem_addr));
    __m128 v2 = _mm_cvtpd_ps(load_double(mem_addr + 2));
    return _mm_shuffle_ps(v1, v2, _MM_SHUFFLE(1, 0, 1, 0));
}

inline void store_float_to_double(double *mem_addr, const vec_float &v) {
    // Write low floats
    store_double(mem_addr, _mm_cvtps_pd(v));

    // Write high floats
    store_double(mem_addr + 2, _mm_cvtps_pd(_mm_shuffle_ps(v, v, _MM_SHUFFLE(3, 2, 3, 2))));
}

inline void store_double_to_float(float *mem_addr, const vec_double &v) {
    _mm_storeu_si64(mem_addr, _mm_castps_si128(_mm_cvtpd_ps(v)));
}

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

// Use dedicated FMA instruction. Returns (a*b) + c
inline vec_float fma_f(const vec_float &a, const vec_float &b, const vec_float &c) {
    return Sleef_fmaf4_sse2(a, b, c);
}
inline vec_double fma_d(const vec_double &a, const vec_double &b, const vec_double &c) {
    return Sleef_fmad2_sse2(a, b, c);
}

inline vec_float add_f(const vec_float &a, const vec_float &b) { return _mm_add_ps(a, b); }
inline vec_double add_d(const vec_double &a, const vec_double &b) { return _mm_add_pd(a, b); }

inline vec_float sub_f(const vec_float &a, const vec_float &b) { return _mm_sub_ps(a, b); }
inline vec_double sub_d(const vec_double &a, const vec_double &b) { return _mm_sub_pd(a, b); }

inline vec_float mul_f(const vec_float &a, const vec_float &b) { return _mm_mul_ps(a, b); }
inline vec_double mul_d(const vec_double &a, const vec_double &b) { return _mm_mul_pd(a, b); }

inline vec_float min_f(const vec_float &a, const vec_float &b) { return _mm_min_ps(a, b); }
inline vec_double min_d(const vec_double &a, const vec_double &b) { return _mm_min_pd(a, b); }

inline vec_float max_f(const vec_float &a, const vec_float &b) { return _mm_max_ps(a, b); }
inline vec_double max_d(const vec_double &a, const vec_double &b) { return _mm_max_pd(a, b); }

inline vec_float neg_f(const vec_float &a) {
    // Constant from Eigen "Core/arch/SSE/PacketMath.h"
    vec_float mask =
        _mm_castsi128_ps(_mm_setr_epi32(0x80000000, 0x80000000, 0x80000000, 0x80000000));
    return _mm_xor_ps(a, mask);
}
inline vec_double neg_d(const vec_double &a) {
    // Constant from Eigen "Core/arch/SSE/PacketMath.h"
    vec_double mask = _mm_castsi128_pd(_mm_setr_epi32(0x0, 0x80000000, 0x0, 0x80000000));
    return _mm_xor_pd(a, mask);
}

template <int It> inline vec_float rsqrt_newton_f(const vec_float &a, const vec_float &approx);
template <int It> inline vec_double rsqrt_newton_d(const vec_double &a, const vec_double &approx);

inline vec_float rsqrt_f(const vec_float &a) {
    // Newton iter count from "Core/arch/SSE/MathFunctions.h"
    return rsqrt_newton_f<1>(a, _mm_rsqrt_ps(a));
}
inline vec_double rsqrt_d(const vec_double &a) {
    // I don't think the newton iteration is worth it given the lack of fast approximation in double
    // precision
    return _mm_div_pd(splat_double(1.0), _mm_sqrt_pd(a));
}

#elif BPCELLS_SLEEF_MODE == _BPCELLS_SLEEF_NEON

#define BPCELLS_VEC_FLOAT_SIZE 4
#define BPCELLS_VEC_DOUBLE_SIZE 2

using vec_float = float32x4_t;
using vec_double = float64x2_t;

inline vec_float load_float(float const *mem_addr) { return vld1q_f32(mem_addr); }
inline void store_float(float *mem_addr, const vec_float &v) { vst1q_f32(mem_addr, v); }
inline vec_float splat_float(float f) { return vdupq_n_f32(f); }

inline vec_double load_double(double const *mem_addr) { return vld1q_f64(mem_addr); }
inline void store_double(double *mem_addr, const vec_double &v) { vst1q_f64(mem_addr, v); }
inline vec_double splat_double(double f) { return vdupq_n_f64(f); }

inline vec_float load_double_to_float(double *mem_addr) {
    float32x2_t v1 = vcvt_f32_f64(load_double(mem_addr));
    float32x2_t v2 = vcvt_f32_f64(load_double(mem_addr + 2));
    return vcombine_f32(v1, v2);
}

inline void store_float_to_double(double *mem_addr, const vec_float &v) {
    // Write low floats
    store_double(mem_addr, vcvt_high_f64_f32(vextq_f32(v, v, 2)));

    // Write high floats
    store_double(mem_addr + 2, vcvt_high_f64_f32(v));
}

inline void store_double_to_float(float *mem_addr, const vec_double &v) {
    vst1_f32(mem_addr, vcvt_f32_f64(v));
}

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

// Use dedicated FMA instruction. Returns (a*b) + c
inline vec_float fma_f(const vec_float &a, const vec_float &b, const vec_float &c) {
    return Sleef_fmaf4_advsimd(a, b, c);
}
inline vec_double fma_d(const vec_double &a, const vec_double &b, const vec_double &c) {
    return Sleef_fmad2_advsimd(a, b, c);
}

inline vec_float add_f(const vec_float &a, const vec_float &b) { return vaddq_f32(a, b); }
inline vec_double add_d(const vec_double &a, const vec_double &b) { return vaddq_f64(a, b); }

inline vec_float sub_f(const vec_float &a, const vec_float &b) { return vsubq_f32(a, b); }
inline vec_double sub_d(const vec_double &a, const vec_double &b) { return vsubq_f64(a, b); }

inline vec_float mul_f(const vec_float &a, const vec_float &b) { return vmulq_f32(a, b); }
inline vec_double mul_d(const vec_double &a, const vec_double &b) { return vmulq_f64(a, b); }

inline vec_float min_f(const vec_float &a, const vec_float &b) { return vminq_f32(a, b); }
inline vec_double min_d(const vec_double &a, const vec_double &b) { return vminq_f64(a, b); }

inline vec_float max_f(const vec_float &a, const vec_float &b) { return vmaxq_f32(a, b); }
inline vec_double max_d(const vec_double &a, const vec_double &b) { return vmaxq_f64(a, b); }

inline vec_float neg_f(const vec_float &a) { return vnegq_f32(a); }
inline vec_double neg_d(const vec_double &a) { return vnegq_f64(a); }

template <int It> inline vec_float rsqrt_newton_f(const vec_float &a, const vec_float &approx);
template <int It> inline vec_double rsqrt_newton_d(const vec_double &a, const vec_double &approx);

inline vec_float rsqrt_f(const vec_float &a) {
    // Newton iter count from "Core/arch/NEON/MathFunctions.h"
    return rsqrt_newton_f<2>(a, vrsqrteq_f32(a));
}
inline vec_double rsqrt_d(const vec_double &a) {
    // Newton iter count from "Core/arch/NEON/MathFunctions.h"
    return rsqrt_newton_d<3>(a, vrsqrteq_f64(a));
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
inline void store_float(float *mem_addr, const vec_float &v) {
    mem_addr[0] = v.x0;
    mem_addr[1] = v.x1;
    mem_addr[2] = v.x2;
    mem_addr[3] = v.x3;
}
inline vec_float splat_float(float f) { return {f, f, f, f}; }

inline vec_double load_double(double const *mem_addr) {
    return {mem_addr[0], mem_addr[1], mem_addr[2], mem_addr[3]};
}
inline void store_double(double *mem_addr, const vec_double &v) {
    mem_addr[0] = v.x0;
    mem_addr[1] = v.x1;
    mem_addr[2] = v.x2;
    mem_addr[3] = v.x3;
}
inline vec_double splat_double(double f) { return {f, f, f, f}; }

inline vec_float load_double_to_float(double *mem_addr) {
    return {
        static_cast<float>(mem_addr[0]),
        static_cast<float>(mem_addr[1]),
        static_cast<float>(mem_addr[2]),
        static_cast<float>(mem_addr[3])};
}

inline void store_float_to_double(double *mem_addr, const vec_float &v) {
    mem_addr[0] = v.x0;
    mem_addr[1] = v.x1;
    mem_addr[2] = v.x2;
    mem_addr[3] = v.x3;
}

inline void store_double_to_float(float *mem_addr, const vec_double &v) {
    mem_addr[0] = v.x0;
    mem_addr[1] = v.x1;
    mem_addr[2] = v.x2;
    mem_addr[3] = v.x3;
}

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

// Use dedicated FMA instruction. Returns (a*b) + c
inline vec_float fma_f(const vec_float &a, const vec_float &b, const vec_float &c) {
    return {a.x0 * b.x0 + c.x0, a.x1 * b.x1 + c.x1, a.x2 * b.x2 + c.x2, a.x3 * b.x3 + c.x3};
}
inline vec_double fma_d(const vec_double &a, const vec_double &b, const vec_double &c) {
    return {a.x0 * b.x0 + c.x0, a.x1 * b.x1 + c.x1, a.x2 * b.x2 + c.x2, a.x3 * b.x3 + c.x3};
}

inline vec_float add_f(const vec_float &a, const vec_float &b) {
    return {a.x0 + b.x0, a.x1 + b.x1, a.x2 + b.x2, a.x3 + b.x3};
}
inline vec_double add_d(const vec_double &a, const vec_double &b) {
    return {a.x0 + b.x0, a.x1 + b.x1, a.x2 + b.x2, a.x3 + b.x3};
}

inline vec_float sub_f(const vec_float &a, const vec_float &b) {
    return {a.x0 - b.x0, a.x1 - b.x1, a.x2 - b.x2, a.x3 - b.x3};
}
inline vec_double sub_d(const vec_double &a, const vec_double &b) {
    return {a.x0 - b.x0, a.x1 - b.x1, a.x2 - b.x2, a.x3 - b.x3};
}

inline vec_float mul_f(const vec_float &a, const vec_float &b) {
    return {a.x0 * b.x0, a.x1 * b.x1, a.x2 * b.x2, a.x3 * b.x3};
}
inline vec_double mul_d(const vec_double &a, const vec_double &b) {
    return {a.x0 * b.x0, a.x1 * b.x1, a.x2 * b.x2, a.x3 * b.x3};
}

inline vec_float min_f(const vec_float &a, const vec_float &b) {
    return {std::min(a.x0, b.x0), std::min(a.x1, b.x1), std::min(a.x2, b.x2), std::min(a.x3, b.x3)};
}
inline vec_double min_d(const vec_double &a, const vec_double &b) {
    return {std::min(a.x0, b.x0), std::min(a.x1, b.x1), std::min(a.x2, b.x2), std::min(a.x3, b.x3)};
}

inline vec_float max_f(const vec_float &a, const vec_float &b) {
    return {std::max(a.x0, b.x0), std::max(a.x1, b.x1), std::max(a.x2, b.x2), std::max(a.x3, b.x3)};
}
inline vec_double max_d(const vec_double &a, const vec_double &b) {
    return {std::max(a.x0, b.x0), std::max(a.x1, b.x1), std::max(a.x2, b.x2), std::max(a.x3, b.x3)};
}

inline vec_float neg_f(const vec_float &a) { return {-a.x0, -a.x1, -a.x2, -a.x3}; }
inline vec_double neg_d(const vec_double &a) { return {-a.x0, -a.x1, -a.x2, -a.x3}; }

inline vec_float rsqrt_f(const vec_float &a) {
    return {1 / sqrtf(a.x0), 1 / sqrtf(a.x1), 1 / sqrtf(a.x2), 1 / sqrtf(a.x3)};
}
inline vec_double rsqrt_d(const vec_double &a) {
    return {1 / sqrt(a.x0), 1 / sqrt(a.x1), 1 / sqrt(a.x2), 1 / sqrt(a.x3)};
}

#endif

// Out-of-line definitions of rsqrt Netwon-Raphson iterations

// Improve a sqrt estimate via Newton-Raphson iterations
template <int It> inline vec_float rsqrt_newton_f(const vec_float &a, const vec_float &approx) {

    // Comment from Eigen "Core/MathFunctionsImpl.h":
    // Refine the approximation using one Newton-Raphson step:
    //   x_{n+1} = x_n * (1.5 + (-0.5 * x_n) * (a * x_n)).
    // The approximation is expressed this way to avoid over/under-flows.
    const vec_float one_point_five = splat_float(1.5);
    const vec_float minus_half = splat_float(-0.5);
    vec_float newton =
        mul_f(approx, fma_f(mul_f(minus_half, approx), mul_f(a, approx), one_point_five));

    for (int i = 1; i < It; i++) {
        newton = mul_f(newton, fma_f(mul_f(minus_half, newton), mul_f(a, newton), one_point_five));
    }
    return newton;
}

// Improve a sqrt estimate via Newton-Raphson iterations
template <int It> inline vec_double rsqrt_newton_d(const vec_double &a, const vec_double &approx) {
    // Comment from Eigen "Core/MathFunctionsImpl.h":
    // Refine the approximation using one Newton-Raphson step:
    //   x_{n+1} = x_n * (1.5 + (-0.5 * x_n) * (a * x_n)).
    // The approximation is expressed this way to avoid over/under-flows.
    const vec_double one_point_five = splat_double(1.5);
    const vec_double minus_half = splat_double(-0.5);
    vec_double newton =
        mul_d(approx, fma_d(mul_d(minus_half, approx), mul_d(a, approx), one_point_five));

    for (int i = 1; i < It; i++) {
        newton = mul_d(newton, fma_d(mul_d(minus_half, newton), mul_d(a, newton), one_point_five));
    }
    return newton;
}

} // end namespace BPCells