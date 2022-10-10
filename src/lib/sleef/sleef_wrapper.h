#pragma once

// Wrapper to handle providing sleef math functions in a uniform way to downstream users
#include <cstdint>

// Note: may need to improve
#if defined(__AVX2__)

    #include <x86intrin.h>
    #include <cstring>

    // This nonsense is to avoid getting -Wunused-function warnings from inside sleef
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-function"
    #pragma GCC diagnostic ignored "-Wpedantic"
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wunused-function"
    #include "sleefinline_avx2.h"
    #pragma clang diagnostic pop
    #pragma GCC diagnostic pop

    #define BPCELLS_F32_VEC_SIZE 8
    #define BPCELLS_F64_VEC_SIZE 4


    inline void bpcells_log1pf_vec(const float* in, float *out) {
        auto x = _mm256_loadu_ps(in);
        _mm256_storeu_ps(out, Sleef_log1pf8_u10avx2(x));
    }

    inline void bpcells_log1pd_vec(const double* in, double *out) {
        auto x = _mm256_loadu_pd(in);
        _mm256_storeu_pd(out, Sleef_log1pd4_u10avx2(x));
    }

#elif defined(__AVX__)

    #include <x86intrin.h>
    #include <cstring>

    // This nonsense is to avoid getting -Wunused-function warnings from inside sleef
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-function"
    #pragma GCC diagnostic ignored "-Wpedantic"
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wunused-function"
    #include "sleefinline_avx.h"
    #pragma clang diagnostic pop
    #pragma GCC diagnostic pop
    

    #define BPCELLS_F32_VEC_SIZE 8
    #define BPCELLS_F64_VEC_SIZE 4

    inline void bpcells_log1pf_vec(const float* in, float *out) {
        auto x = _mm256_loadu_ps(in);
        _mm256_storeu_ps(out, Sleef_log1pf8_u10avx2(x));
    }

    inline void bpcells_log1pd_vec(const double* in, double *out) {
        auto x = _mm256_loadu_pd(in);
        _mm256_storeu_pd(out, Sleef_log1pd4_u10avx2(x));
    }

#elif defined(__SSE2__)

    #include <x86intrin.h>
    #include <cstring>

    // This nonsense is to avoid getting -Wunused-function warnings from inside sleef
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-function"
    #pragma GCC diagnostic ignored "-Wpedantic"
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wunused-function"
    #include "sleefinline_sse2.h"
    #pragma clang diagnostic pop
    #pragma GCC diagnostic pop
    

    #define BPCELLS_F32_VEC_SIZE 4
    #define BPCELLS_F64_VEC_SIZE 2

    inline void bpcells_log1pf_vec(const float* in, float *out) {
        auto x = _mm_loadu_ps(in);
        _mm_storeu_ps(out, Sleef_log1pf4_u10sse2(x));
    }

    inline void bpcells_log1pd_vec(const double* in, double *out) {
        auto x = _mm_loadu_pd(in);
        _mm_storeu_pd(out, Sleef_log1pd2_u10sse2(x));
    }


#elif defined(__ARM_NEON)

    #include <arm_neon.h>
    #include <cstring>

    // This nonsense is to avoid getting -Wunused-function warnings from inside sleef
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-function"
    #pragma GCC diagnostic ignored "-Wpedantic"
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wunused-function"
    #include "sleefinline_advsimd.h"
    #pragma clang diagnostic pop
    #pragma GCC diagnostic pop
    
    
    #define BPCELLS_F32_VEC_SIZE 4
    #define BPCELLS_F64_VEC_SIZE 2

    inline void bpcells_log1pf_vec(const float* in, float *out) {
        auto x = vld1q_f32(in);
        vst1q_f32(out, Sleef_log1pf4_u10advsimd(x));
    }

    inline void bpcells_log1pd_vec(const double* in, double *out) {
        auto x = vld1q_f64(in);
        vst1q_f64(out, Sleef_log1pd2_u10advsimd(x));
    }


#else
    
    #include <cmath>
    #define BPCELLS_F32_VEC_SIZE 4
    #define BPCELLS_F64_VEC_SIZE 2

    inline void bpcells_log1pf_vec(const float* in, float *out) {
        out[0] = log1pf(in[0]);
        out[1] = log1pf(in[1]);
        out[2] = log1pf(in[2]);
        out[3] = log1pf(in[3]);
    }

    inline void bpcells_log1pd_vec(const double* in, double *out) {
        out[0] = log1p(in[0]);
        out[1] = log1p(in[1]);
    }

#endif


