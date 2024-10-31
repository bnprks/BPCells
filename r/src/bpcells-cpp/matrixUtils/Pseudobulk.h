// Copyright 2024 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once
#include "../matrixIterators/MatrixIterator.h"
namespace BPCells {

// Enum to choose which statistic to compute for pseudobulk matrix creation, with bitwise operators.
enum class PseudobulkStatsMethod { 
    NonZeros = 1 << 0,
    Sum = 1 << 1,
    Mean = 1 << 2,
    Variance = 1 << 3
};

inline PseudobulkStatsMethod operator&(PseudobulkStatsMethod a, PseudobulkStatsMethod b) {
    return static_cast<PseudobulkStatsMethod>(static_cast<int>(a) & static_cast<int>(b));
}
inline PseudobulkStatsMethod operator|(PseudobulkStatsMethod a, PseudobulkStatsMethod b) {
    return static_cast<PseudobulkStatsMethod>(static_cast<int>(a) | static_cast<int>(b));
}

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
// - method: method to use for calculating pseudobulk stats, constructed using bitflags from PseudobulkStatsMethod.
// - transpose: Whether the matrix is transposed. If false, groups columns and calculates stats per-row. If true, groups rows and calculates stats per-column.
// Returns:
// - PseudobulkStats struct containing each of the requested pseudobulk aggregations, in the shape of (num_features, num_groups).
template <typename T>
PseudobulkStats pseudobulk_matrix(std::unique_ptr<MatrixLoader<T>> &&mat,
                                  const std::vector<uint32_t>& cell_groups,
                                  PseudobulkStatsMethod method,
                                  bool transpose,
                                  std::atomic<bool> *user_interrupt);
} // namespace BPCells
