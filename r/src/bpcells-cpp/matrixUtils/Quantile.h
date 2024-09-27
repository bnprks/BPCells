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

// Find the `quantile`th value for each column in an IterableMatrix.
// Uses type 7 quantile calculation, which is the default in R.
// Args:
// - mat: matrix to compute quantiles from
// - quantile: quantile to compute from each column, between [0,1]
// Returns:
// - A vector of quantile values, one for each column.
template <typename T>
std::vector<T> matrix_quantile_per_col(std::unique_ptr<MatrixLoader<T>>&& mat, 
                                       double quantile, 
                                       std::atomic<bool> *user_interrupt);
} // namespace BPCells