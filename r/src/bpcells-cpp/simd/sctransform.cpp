// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "simd/sctransform.cpp"
#include <hwy/foreach_target.h>

#include <hwy/contrib/math/math-inl.h>
#include <hwy/highway.h>

#include "transform-inl.h"
#include "sctransform.h"

HWY_BEFORE_NAMESPACE();

namespace BPCells::simd::HWY_NAMESPACE {

namespace hn = hwy::HWY_NAMESPACE;

namespace {

// Helper to handle Newton-Raphson iterations
template <int It, class V> inline V ReciprocalSqrtNewtonIter(const V a, const V approx) {
    hn::DFromV<V> d;
    V minus_half = hn::Set(d, -0.5);
    V one_point_five = hn::Set(d, 1.5);

    // Comment from Eigen "Core/MathFunctionsImpl.h":
    // Refine the approximation using one Newton-Raphson step:
    //   x_{n+1} = x_n * (1.5 + (-0.5 * x_n) * (a * x_n)).
    // The approximation is expressed this way to avoid over/under-flows.
    V newton = hn::Mul(
        approx, hn::MulAdd(hn::Mul(minus_half, approx), hn::Mul(a, approx), one_point_five)
    );
    if constexpr (It > 1) {
        newton = ReciprocalSqrtNewtonIter<It-1>(a, newton);
    }
    return newton;
}

} // namespace

// Reasonably accurate reciprocal sqrt implementation using Newton-Raphson iterations to improve the
// fast approximations. Makes a guess on how many iterations to use based on the architecture
template <class V> inline V ReciprocalSqrt(V a) {
    V approx = hn::ApproximateReciprocalSqrt(a);

    if constexpr (HWY_ARCH_ARM) {
        if constexpr (std::is_same_v<hn::TFromV<V>, float>) {
            // 2 iter for ARM float
            return ReciprocalSqrtNewtonIter<2>(a, approx);
        } else {
            // 3 iter for ARM double
            return ReciprocalSqrtNewtonIter<3>(a, approx);
        }
    } else {
        // 1 iter for Intel and other architectures (WASM might not actually need an iteration)
        return ReciprocalSqrtNewtonIter<1>(a, approx);
    }
}

/**
 * @brief Vectorized SCTransform(0) Pearson residuals calculation
 *
 * Normalized values are calculated as `(X - mu) / sqrt(mu + mu^2/theta)` where `mu` is calculated
 * as `cell_factor * gene_beta`. Additionally, the denominator (variance) can be clipped at a minimu
 * threshold, and the final value can be clipped at both min and max thresholds. (Only min is needed
 * for this function, since we assume the clip max is positive and mu is negative)
 *
 * @tparam V
 * @param cell_factor Cell size factor, often number of reads
 * @param gene_beta Gene abundance factor, often fraction of total reads
 * @param theta_inv 1/theta
 * @param sd_inv_max max threshold for `1/sqrt(mu + mu^2/theta)`
 * @param clip_min min threshold for final output
 * @return V
 */
template <class V>
inline V SCTransformZero(V cell_factor, V gene_beta, V theta_inv, V sd_inv_max, V clip_min) {
    // mu = cell_factor * gene_beta
    V mu = hn::Mul(cell_factor, gene_beta);
    // 1/sd = 1/sqrt(mu + mu*mu*theta_inv)
    V sd_inv = ReciprocalSqrt(hn::MulAdd(hn::Mul(mu, theta_inv), mu, mu));
    // Apply clipping to standard deviation
    sd_inv = hn::Min(sd_inv_max, sd_inv);

    // -mu / sd
    V val = hn::Mul(hn::Neg(mu), sd_inv);
    // Apply final clipping
    return hn::Max(clip_min, val);
}

/**
 * @brief Vectorized SCTransform(X) - SCTransform(0) Pearson residuals calculation
 *
 * Normalized values are calculated as `(X - mu) / sqrt(mu + mu^2/theta)` where `mu` is calculated
 * as `cell_factor * gene_beta`. Additionally, the denominator (variance) can be clipped at a minimu
 * threshold, and the final value can be clipped at both min and max thresholds.
 *
 * @tparam V
 * @param x raw counts value X
 * @param cell_factor Cell size factor, often number of reads
 * @param gene_beta Gene abundance factor, often fraction of total reads
 * @param theta_inv 1/theta
 * @param sd_inv_max max threshold for `1/sqrt(mu + mu^2/theta)`
 * @param clip_min min threshold for final output
 * @param clip_max max threshold for final output
 * @return V
 */
template <class V>
inline V SCTransformZeroSubtracted(
    V x, V cell_factor, V gene_beta, V theta_inv, V sd_inv_max, V clip_min, V clip_max
) {
    
    // mu = cell_factor * gene_beta
    V mu = hn::Mul(cell_factor, gene_beta);
    // 1/sd = 1/sqrt(mu + mu*mu*theta_inv)
    V sd_inv = ReciprocalSqrt(hn::MulAdd(hn::Mul(mu, theta_inv), mu, mu));
    // Apply clipping to standard deviation
    sd_inv = hn::Min(sd_inv_max, sd_inv);

    // Calculate the zero value
    V zero_val = hn::Max(hn::Mul(hn::Neg(mu), sd_inv), clip_min);

    // (X - mu) / sd
    V val = hn::Mul(hn::Sub(x, mu), sd_inv);

    // Apply clamping
    val = hn::Clamp(val, clip_min, clip_max);

    // Subtract out zero value
    return hn::Sub(val, zero_val);
}

// Load sctransform(inout) - sctransform(0) into inout
// Convert to single precision internally for faster computations
void sctransform_zero_subtracted(
    double *HWY_RESTRICT inout,
    const float cell_factor_scalar,
    const uint32_t *HWY_RESTRICT row,
    const float *HWY_RESTRICT gene_beta,
    const float *HWY_RESTRICT theta_inv,
    const SCTransformClipParam &clip,
    size_t n
) {
    return transform3_sparse_downcast(
        inout,
        row,
        n,
        gene_beta,
        theta_inv,
        [cell_factor_scalar, &clip](auto d, auto x, auto gene_beta, auto theta_inv) {
            auto cell_factor = hn::Set(d, cell_factor_scalar);
            auto sd_inv_max = hn::Set(d, clip.sd_inv_max);
            auto clip_min = hn::Set(d, clip.clip_min);
            auto clip_max = hn::Set(d, clip.clip_max);
            return SCTransformZeroSubtracted(
                x, cell_factor, gene_beta, theta_inv, sd_inv_max, clip_min, clip_max
            );
        }
    );
}

// Same as sctransform_zero_subtracted, but transposed so columns are genes
// and rows are cells
void sctransform_zero_subtracted_transpose(
    double *HWY_RESTRICT inout,
    const float *HWY_RESTRICT cell_factor,
    const uint32_t *HWY_RESTRICT row,
    const float gene_beta_scalar,
    const float theta_inv_scalar,
    const SCTransformClipParam &clip,
    size_t n
) {
    return transform2_sparse_downcast(
        inout,
        row,
        n,
        cell_factor,
        [gene_beta_scalar, theta_inv_scalar, &clip](auto d, auto x, auto cell_factor) {
            auto gene_beta = hn::Set(d, gene_beta_scalar);
            auto theta_inv = hn::Set(d, theta_inv_scalar);
            auto sd_inv_max = hn::Set(d, clip.sd_inv_max);
            auto clip_min = hn::Set(d, clip.clip_min);
            auto clip_max = hn::Set(d, clip.clip_max);
            return SCTransformZeroSubtracted(
                x, cell_factor, gene_beta, theta_inv, sd_inv_max, clip_min, clip_max
            );
        }
    );
}

void sctransform_load_zero(
    double *HWY_RESTRICT out,
    const float cell_factor_scalar,
    const float *HWY_RESTRICT gene_beta,
    const float *HWY_RESTRICT theta_inv,
    const SCTransformClipParam &clip,
    size_t n
) {
    transform3_downcast(
        out,
        n,
        gene_beta,
        theta_inv,
        [cell_factor_scalar, &clip](auto d, auto x, auto gene_beta, auto theta_inv) {
            auto cell_factor = hn::Set(d, cell_factor_scalar);
            auto sd_inv_max = hn::Set(d, clip.sd_inv_max);
            auto clip_min = hn::Set(d, clip.clip_min);
            return SCTransformZero(cell_factor, gene_beta, theta_inv, sd_inv_max, clip_min);
        }
    );
}

void sctransform_load_zero_transpose(
    double *HWY_RESTRICT out,
    const float *HWY_RESTRICT cell_factor,
    const float gene_beta_scalar,
    const float theta_inv_scalar,
    const SCTransformClipParam &clip,
    size_t n
) {
    transform2_downcast(
        out,
        n,
        cell_factor,
        [gene_beta_scalar, theta_inv_scalar, &clip](auto d, auto x, auto cell_factor) {
            auto gene_beta = hn::Set(d, gene_beta_scalar);
            auto theta_inv = hn::Set(d, theta_inv_scalar);
            auto sd_inv_max = hn::Set(d, clip.sd_inv_max);
            auto clip_min = hn::Set(d, clip.clip_min);
            return SCTransformZero(cell_factor, gene_beta, theta_inv, sd_inv_max, clip_min);
        }
    );
}


// Multiply a vector entry by a column of the zero-transformed matrix, and accumulate into out
void sctransform_multiply_right_zero(
    float *HWY_RESTRICT out,
    const float vec_entry,
    const float cell_factor_scalar,
    const float *HWY_RESTRICT gene_beta,
    const float *HWY_RESTRICT theta_inv,
    const SCTransformClipParam &clip,
    size_t rows
) {
    transform3(
        out,
        rows,
        gene_beta,
        theta_inv,
        [vec_entry, cell_factor_scalar, &clip](auto d, auto x, auto gene_beta, auto theta_inv) {
            auto vec = hn::Set(d, vec_entry);
            auto cell_factor = hn::Set(d, cell_factor_scalar);
            auto sd_inv_max = hn::Set(d, clip.sd_inv_max);
            auto clip_min = hn::Set(d, clip.clip_min);
            return hn::MulAdd(vec, SCTransformZero(cell_factor, gene_beta, theta_inv, sd_inv_max, clip_min), x);
        }
    );
}

// Multiply a vector entry by a row of the zero-transformed matrix, and accumulate into out
void sctransform_multiply_left_zero(
    float *HWY_RESTRICT out,
    const float vec_entry,
    const float *HWY_RESTRICT cell_factor,
    const float gene_beta_scalar,
    const float theta_inv_scalar,
    const SCTransformClipParam &clip,
    size_t cols
) {
    transform2(
        out,
        cols,
        cell_factor,
        [vec_entry, gene_beta_scalar, theta_inv_scalar, &clip](auto d, auto x, auto cell_factor) {
            auto vec = hn::Set(d, vec_entry);
            auto gene_beta = hn::Set(d, gene_beta_scalar);
            auto theta_inv = hn::Set(d, theta_inv_scalar);
            auto sd_inv_max = hn::Set(d, clip.sd_inv_max);
            auto clip_min = hn::Set(d, clip.clip_min);
            return hn::MulAdd(vec, SCTransformZero(cell_factor, gene_beta, theta_inv, sd_inv_max, clip_min), x);
        }
    );
}

} // namespace BPCells::simd::HWY_NAMESPACE
HWY_AFTER_NAMESPACE();

#if HWY_ONCE

namespace BPCells::simd {

HWY_EXPORT(sctransform_zero_subtracted);
HWY_EXPORT(sctransform_load_zero);
HWY_EXPORT(sctransform_zero_subtracted_transpose);
HWY_EXPORT(sctransform_load_zero_transpose);
HWY_EXPORT(sctransform_multiply_right_zero);
HWY_EXPORT(sctransform_multiply_left_zero);

void sctransform_zero_subtracted(
    double *HWY_RESTRICT inout,
    const float cell_factor,
    const uint32_t *HWY_RESTRICT row,
    const float *HWY_RESTRICT gene_beta,
    const float *HWY_RESTRICT theta_inv,
    const SCTransformClipParam &clip,
    size_t n
) {
    return HWY_DYNAMIC_DISPATCH(sctransform_zero_subtracted)(inout, cell_factor, row, gene_beta, theta_inv, clip, n);
}

void sctransform_load_zero(
    double *HWY_RESTRICT out,
    const float cell_factor,
    const float *HWY_RESTRICT gene_beta,
    const float *HWY_RESTRICT theta_inv,
    const SCTransformClipParam &clip,
    size_t n
) {
    return HWY_DYNAMIC_DISPATCH(sctransform_load_zero)(out, cell_factor, gene_beta, theta_inv, clip, n);
}

void sctransform_zero_subtracted_transpose(
    double *HWY_RESTRICT inout,
    const float *HWY_RESTRICT cell_factor,
    const uint32_t *HWY_RESTRICT row,
    const float gene_beta,
    const float theta_inv,
    const SCTransformClipParam &clip,
    size_t n
) {
    return HWY_DYNAMIC_DISPATCH(sctransform_zero_subtracted_transpose)(inout, cell_factor, row, gene_beta, theta_inv, clip, n);
}

void sctransform_load_zero_transpose(
    double *HWY_RESTRICT out,
    const float *HWY_RESTRICT cell_factor,
    const float gene_beta,
    const float theta_inv,
    const SCTransformClipParam &clip,
    size_t n
) {
    return HWY_DYNAMIC_DISPATCH(sctransform_load_zero_transpose)(out, cell_factor, gene_beta, theta_inv, clip, n);
}

void sctransform_multiply_right_zero(
    float *HWY_RESTRICT out,
    const float vec_entry,
    const float cell_factor,
    const float *HWY_RESTRICT gene_beta,
    const float *HWY_RESTRICT theta_inv,
    const SCTransformClipParam &clip,
    size_t rows
) {
    return HWY_DYNAMIC_DISPATCH(sctransform_multiply_right_zero)(out, vec_entry, cell_factor, gene_beta, theta_inv, clip, rows);
}

void sctransform_multiply_left_zero(
    float *HWY_RESTRICT out,
    const float vec_entry,
    const float *HWY_RESTRICT cell_factor,
    const float gene_beta,
    const float theta_inv,
    const SCTransformClipParam &clip,
    size_t cols
) {
    return HWY_DYNAMIC_DISPATCH(sctransform_multiply_left_zero)(out, vec_entry, cell_factor, gene_beta, theta_inv, clip, cols);
}

} // namespace BPCells::simd
#endif