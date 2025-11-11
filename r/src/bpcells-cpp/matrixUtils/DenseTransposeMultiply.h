// Copyright 2025 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <atomic>
#include <memory>

#ifndef RCPP_EIGEN
#include <Eigen/Dense>
#else
#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <RcppEigen.h>
#endif

#include "../matrixIterators/MatrixIterator.h"


/**
 * Calculate mat * t(mat), returning a dense matrix output
 * @param mat Unique pointer to the matrix
 * @param buffer_bytes Number of bytes to use for buffering matrix columns. Must be at least
 * mat.rows()*48. Recommended to be just below the size of the largest CPU cache (several MB)
 * @param threads Number of threads to use
 * @param user_interrupt Signal for user-requested early termination of computations
 * 
 * This uses a one-pass algorithm to load contiguous column chunks of mat, perform sparse-sparse
 * matrix multiplies in memory and update a dense output matrix.
 */
Eigen::MatrixXd dense_transpose_multiply(
    std::unique_ptr<BPCells::MatrixLoader<double>> &&mat,
    size_t buffer_bytes,
    size_t threads,
    std::atomic<bool> *user_interrupt
);