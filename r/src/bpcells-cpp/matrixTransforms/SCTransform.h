// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <atomic>
#include <algorithm>

#include "MatrixTransform.h"

namespace BPCells {

// Pearson residuals of an SCTransform model, with SIMD optimizations
// There is support for SCT models with no continuous covariates (i.e.
// there are only a few unique covariate vectors repeated across cells)
//
// Rows are genes/features
// First row of row params is 1/theta. Second row is per-gene beta parameters (linear scale).
//
// Cols are cells
// Cell params are the linear UMI counts per cell
//
// Global params are: max inverse standard deviation, clip min, clip max
// For cell c, gene g:
// mu = row_params[1, g] * col_params[0,c]
// theta_inv = row_params[0,g]
// sd = sqrt(mu + mu^2/theta_inv)
// Transform = (X - mu) / min(max_sd, sd)
// Then transform is also clipped according to clip-min and clip-max
class SCTransformPearsonSIMD : public MatrixTransformDense {
  private:
    Eigen::ArrayXf theta_inv;
    // Linear-scale counts per cell
    Eigen::ArrayXf cell_read_counts; 
    // Linear-scale gene betas
    Eigen::ArrayXf gene_beta;  

    float sd_inv_max, clip_min, clip_max;
  public:
    SCTransformPearsonSIMD(std::unique_ptr<MatrixLoader<double>> &&loader, TransformFit fit);

    bool loadZeroSubtracted(MatrixLoader<double> &loader) override;
    void loadZero(double *values, uint32_t count, uint32_t start_row, uint32_t col) override;

    void vecMultiplyRightZero(
        Eigen::VectorXd &out,
        const Eigen::Map<Eigen::VectorXd> v,
        std::atomic<bool> *user_interrupt = NULL
    ) override;

    void vecMultiplyLeftZero(
        Eigen::VectorXd &out,
        const Eigen::Map<Eigen::VectorXd> v,
        std::atomic<bool> *user_interrupt = NULL
    ) override;
};

// Same as before, but now assuming rows are cells and cols are genes
class SCTransformPearsonTransposeSIMD : public MatrixTransformDense {
  private:
    Eigen::ArrayXf theta_inv;
    // Linear-scale counts per cell
    Eigen::ArrayXf cell_read_counts; 
    // Linear-scale gene betas
    Eigen::ArrayXf gene_beta;  

    float sd_inv_max, clip_min, clip_max;
  public:
    SCTransformPearsonTransposeSIMD(std::unique_ptr<MatrixLoader<double>> &&loader, TransformFit fit);

    bool loadZeroSubtracted(MatrixLoader<double> &loader) override;
    void loadZero(double *values, uint32_t count, uint32_t start_row, uint32_t col) override;

    void vecMultiplyRightZero(
        Eigen::VectorXd &out,
        const Eigen::Map<Eigen::VectorXd> v,
        std::atomic<bool> *user_interrupt = NULL
    ) override;

    void vecMultiplyLeftZero(
        Eigen::VectorXd &out,
        const Eigen::Map<Eigen::VectorXd> v,
        std::atomic<bool> *user_interrupt = NULL
    ) override;
};

// Non-SIMD, full precision variants
class SCTransformPearson : public MatrixTransformDense {
  public:
    using MatrixTransformDense::MatrixTransformDense;

    bool loadZeroSubtracted(MatrixLoader<double> &loader) override;
    void loadZero(double *values, uint32_t count, uint32_t start_row, uint32_t col) override;
};

class SCTransformPearsonTranspose : public MatrixTransformDense {
  public:
    using MatrixTransformDense::MatrixTransformDense;

    bool loadZeroSubtracted(MatrixLoader<double> &loader) override;
    void loadZero(double *values, uint32_t count, uint32_t start_row, uint32_t col) override;
};

} // end namespace BPCells