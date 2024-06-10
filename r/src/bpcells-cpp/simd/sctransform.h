// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <hwy/base.h>

namespace BPCells::simd {

struct SCTransformClipParam {
    float sd_inv_max, clip_min, clip_max;
};

// Load sctransform(inout) - sctransform(0) into inout
// Convert to single precision internally for faster computations
void sctransform_zero_subtracted(
    double *HWY_RESTRICT inout,
    const float cell_factor,
    const uint32_t *HWY_RESTRICT row,
    const float *HWY_RESTRICT gene_beta,
    const float *HWY_RESTRICT theta_inv,
    const SCTransformClipParam &clip,
    size_t n
);

// Same as sctransform_zero_subtracted, but transposed so columns are genes
// and rows are cells
void sctransform_zero_subtracted_transpose(
    double *HWY_RESTRICT inout,
    const float *HWY_RESTRICT cell_factor,
    const uint32_t *HWY_RESTRICT row,
    const float gene_beta,
    const float theta_inv,
    const SCTransformClipParam &clip,
    size_t n
);

// Load sctransform(0) into out (data for an individual cell)
void sctransform_load_zero(
    double *HWY_RESTRICT out,
    const float cell_factor,
    const float *HWY_RESTRICT gene_beta,
    const float *HWY_RESTRICT theta_inv,
    const SCTransformClipParam &clip,
    size_t n
);

// Load sctransform(0) into out (data for an individual gene)
void sctransform_load_zero_transpose(
    double *HWY_RESTRICT out,
    const float *HWY_RESTRICT cell_factor,
    const float gene_beta,
    const float theta_inv,
    const SCTransformClipParam &clip,
    size_t n
);

// Multiply a vector entry by a column of the zero-transformed matrix, and add into out
void sctransform_multiply_right_zero(
    float *HWY_RESTRICT out,
    const float vec_entry,
    const float cell_factor,
    const float *HWY_RESTRICT gene_beta,
    const float *HWY_RESTRICT theta_inv,
    const SCTransformClipParam &clip,
    size_t rows
);

// Multiply a vector entry by a row of the zero-transformed matrix, and add into out
void sctransform_multiply_left_zero(
    float *HWY_RESTRICT out,
    const float vec_entry,
    const float *HWY_RESTRICT cell_factor,
    const float gene_beta,
    const float theta_inv,
    const SCTransformClipParam &clip,
    size_t cols
);

} // namespace BPCells::simd