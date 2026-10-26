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

// Caller-owned output arrays of shape (num_features, num_groups), initialized to zero.
// Arrays of statistics that are not calculated may be empty.
struct PseudobulkStats {
    Eigen::Ref<Eigen::ArrayXXd> non_zeros;
    Eigen::Ref<Eigen::ArrayXXd> sum;
    Eigen::Ref<Eigen::ArrayXXd> mean;
    Eigen::Ref<Eigen::ArrayXXd> var;
};

// Take a matrix and do a pseudobulk aggregation on the matrix.
// Note that when calculating a more complex statistic, the simpler ones are also calculated.  The order of complexity is NonZeros < Sum < Mean < Variance.
// Args:
// - mat: Matrix to compute variance from.
// - cell_groups: Mapping of columns to groups.
// - method: method to use for calculating pseudobulk stats, constructed using bitflags from PseudobulkStatsMethod.
// - transpose: Whether the matrix is transposed. If false, groups columns and calculates stats per-row. If true, groups rows and calculates stats per-column.
// - feature_offset: Output row of the first feature in mat. If transposed, disjoint column ranges of a matrix
//   can thus be aggregated in parallel into the same output.
// - res: Output arrays, which hold running means and variances until pseudobulk_finalize() is called.
template <typename T>
void pseudobulk_matrix(std::unique_ptr<MatrixLoader<T>> &&mat,
                       const std::vector<uint32_t>& cell_groups,
                       PseudobulkStatsMethod method,
                       bool transpose,
                       uint32_t feature_offset,
                       PseudobulkStats &res,
                       std::atomic<bool> *user_interrupt);

// Correct the running means and variances in res for the zeros of each group, once all features are aggregated.
void pseudobulk_finalize(PseudobulkStats &res, const std::vector<uint32_t>& cell_groups, PseudobulkStatsMethod method);
} // namespace BPCells
