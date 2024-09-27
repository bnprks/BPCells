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
// Enum to choose which statistic to compute for pseudobulk matrix creation
enum class PseudobulkStatsMethod {
    NonZeros = 0,
    Sum = 1,
    Mean = 2,
    Variance = 3
};
struct PseudobulkStats {
    Eigen::ArrayXXd non_zeros;
    Eigen::ArrayXXd sum;
    Eigen::ArrayXXd mean;
    Eigen::ArrayXXd var;
};


// Take a matrix and do a pseudobulk aggregation on the matrix.
// Note that when calculating a more complex statistic, the simpler ones are also calculated.  The order of complexity is NonZeros < Sum < Mean < Variance.
// Args:
// - mat: Matrix to compute variance from.
// - cell_groups: Mapping of columns to groups.
// - Method to use for calculating pseudobulk stats.
// - transpose: Whether the matrix is transposed. Whether the matrix is transposed. If false, groups columns and calculates stats per-row. If true, groups rows and calculates stats per-column.
// Returns:
// - PseudobulkStats struct containing the requested pseudobulk aggregation, in the shape of (num_features, num_groups).
template <typename T>
PseudobulkStats pseudobulk_matrix(std::unique_ptr<MatrixLoader<T>> &&mat,
                                  const std::vector<uint32_t>& cell_groups,
                                  PseudobulkStatsMethod method,
                                  bool transpose,
                                  std::atomic<bool> *user_interrupt);
} // namespace BPCells
