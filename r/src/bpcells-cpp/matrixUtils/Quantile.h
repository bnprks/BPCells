// Copyright 2024 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#ifndef RCPP_EIGEN
#include <Eigen/SparseCore>
#else
#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <RcppEigen.h>
#endif

#include "../matrixIterators/ColwiseRank.h"
namespace BPCells {


// Find the quantile for each row in an IterableMatrix.
// Args:
// - mat: matrix to compute quantile from
// - quantile: quantile to compute from each row, between [0,1]
// Note: This function will probably crash out for large matrices
template <typename T>
std::vector<T> matrix_quantile_per_row(std::unique_ptr<MatrixLoader<T>>&& mat, 
                                       double quantile, 
                                       std::atomic<bool> *user_interrupt);


// Find the `quantile`th value for each column in an IterableMatrix.
// Args:
// - mat: matrix to compute quantiles from
// - quantile: quantile to compute from each column, between [0,1]
template <typename T>
std::vector<T> matrix_quantile_per_col(std::unique_ptr<MatrixLoader<T>>&& mat, 
                                       double quantile, 
                                       std::atomic<bool> *user_interrupt);
} // namespace BPCells