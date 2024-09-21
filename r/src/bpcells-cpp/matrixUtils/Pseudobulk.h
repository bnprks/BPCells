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

struct PseudobulkStats {
    Eigen::ArrayXXd non_zeros;
    Eigen::ArrayXXd sum;
    Eigen::ArrayXXd mean;
    Eigen::ArrayXXd var;
};

// Enum to choose which statistic to compute for pseudobulk matrix creation
// enum class Stats { None = 0, NonZeroCount = 1, Mean = 2, Variance = 3 };

// Take a matrix and do a pseudobulk aggregation on the matrix.  Av
// Args:
// - mat: Matrix to compute variance from.
// - cell_groups: Mapping of columns to groups.
// - method: Method to use.  One of "non-zeros", "sum", "mean", "variance",
// - transpose: Whether the matrix is transposed.
// Returns:
// - A matrix of variances in the shape of (num_features, num_groups).
template <typename T>
PseudobulkStats pseudobulk_matrix(std::unique_ptr<MatrixLoader<T>> &&mat,
                                  const std::vector<uint32_t>& cell_groups,
                                  const std::string& method,
                                  bool transpose,
                                  std::atomic<bool> *user_interrupt);
} // namespace BPCells
