// Copyright 2024 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "Pseudobulk.h"
namespace BPCells {

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
                                  std::atomic<bool> *user_interrupt) {
    MatrixIterator<T> it(std::move(mat));
    if (transpose && (it.rows() != cell_groups.size())) {
        throw std::invalid_argument("Pseudobulk_matrix(): Cell groups must match the number of columns in the matrix");
    } else if (!transpose && (it.cols() != cell_groups.size())) {
        throw std::invalid_argument("Pseudobulk_matrix(): Cell groups must match the number of columns in the matrix");
    }
    std::vector<double> group_count_vec;
    struct PseudobulkStats res;
    for (auto& cell_group : cell_groups) {
        if (cell_group >= group_count_vec.size()) group_count_vec.resize(cell_group + 1);
        group_count_vec[cell_group]++;
    }
    Eigen::Map<Eigen::Array<double, 1, Eigen::Dynamic>> group_count(group_count_vec.data(), group_count_vec.size());
    // If transposed, also keep matrix transposed with groups as rows for cache-friendliness
    uint32_t num_rows = transpose ? group_count_vec.size() : it.rows();
    uint32_t num_cols = transpose ? it.cols() : group_count_vec.size();
    if ((method & PseudobulkStatsMethod::NonZeros) == PseudobulkStatsMethod::NonZeros) res.non_zeros = Eigen::ArrayXXd::Zero(num_rows, num_cols);
    if ((method & PseudobulkStatsMethod::Sum) == PseudobulkStatsMethod::Sum) res.sum = Eigen::ArrayXXd::Zero(num_rows, num_cols);
    if ((method & PseudobulkStatsMethod::Mean) == PseudobulkStatsMethod::Mean) res.mean = Eigen::ArrayXXd::Zero(num_rows, num_cols);
    if ((method & PseudobulkStatsMethod::Variance)== PseudobulkStatsMethod::Variance) res.var = Eigen::ArrayXXd::Zero(num_rows, num_cols);
    // Use a one pass strategy that can do non-zero, mean, and variance at the same time
    while (it.nextCol()) {
        const uint32_t col = it.currentCol();
        if (user_interrupt != NULL && *user_interrupt) return res;
        while (it.load()) {
            const uint32_t *row_data = it.rowData();
            const T *val_data = it.valData();
            const uint32_t count = it.capacity();
            // Do transpose logic to calculate row_idx, col_idx out of loop, to avoid repeating conditional checks
            if (transpose) {
                // Doing variance calculation also does mean and non-zero count, so specific order is required to make sure we don't double count
                if ((method & PseudobulkStatsMethod::Variance) == PseudobulkStatsMethod::Variance) {
                    for (uint32_t i = 0; i < count; i++) update_mean_and_variance(val_data[i], res.non_zeros(cell_groups[row_data[i]], col), res.mean(cell_groups[row_data[i]], col), res.var(cell_groups[row_data[i]], col));
                } else if ((method & PseudobulkStatsMethod::Mean) == PseudobulkStatsMethod::Mean) {
                    for (uint32_t i = 0; i < count; i++) update_mean(val_data[i], res.non_zeros(cell_groups[row_data[i]], col), res.mean(cell_groups[row_data[i]], col));
                } else if ((method & PseudobulkStatsMethod::NonZeros) == PseudobulkStatsMethod::NonZeros) {
                    for (uint32_t i = 0; i < count; i++) res.non_zeros(cell_groups[row_data[i]], col)++;
                }
                if ((method & PseudobulkStatsMethod::Sum) == PseudobulkStatsMethod::Sum) {
                    for (uint32_t i = 0; i < count; i++) res.sum(cell_groups[row_data[i]], col) += val_data[i];
                }
            } else {
                if ((method & PseudobulkStatsMethod::Variance) == PseudobulkStatsMethod::Variance) {
                    for (uint32_t i = 0; i < count; i++) update_mean_and_variance(val_data[i], res.non_zeros(row_data[i], cell_groups[col]), res.mean(row_data[i], cell_groups[col]), res.var(row_data[i], cell_groups[col]));
                } else if ((method & PseudobulkStatsMethod::Mean) == PseudobulkStatsMethod::Mean) {
                    for (uint32_t i = 0; i < count; i++) update_mean(val_data[i], res.non_zeros(row_data[i], cell_groups[col]), res.mean(row_data[i], cell_groups[col]));
                } else if ((method & PseudobulkStatsMethod::NonZeros) == PseudobulkStatsMethod::NonZeros) {
                    for (uint32_t i = 0; i < count; i++) res.non_zeros(row_data[i], cell_groups[col])++;
                }
                if ((method & PseudobulkStatsMethod::Sum) == PseudobulkStatsMethod::Sum) {
                    for (uint32_t i = 0; i < count; i++) res.sum(row_data[i], cell_groups[col]) += val_data[i];
                }
            }
        }
    }
    if (transpose) {
        //transpose all res to match format of non-transposed res
        if ((method & PseudobulkStatsMethod::NonZeros) == PseudobulkStatsMethod::NonZeros) res.non_zeros.transposeInPlace();
        if ((method & PseudobulkStatsMethod::Sum) == PseudobulkStatsMethod::Sum) res.sum.transposeInPlace();
        if ((method & PseudobulkStatsMethod::Mean) == PseudobulkStatsMethod::Mean) res.mean.transposeInPlace();
        if ((method & PseudobulkStatsMethod::Variance) == PseudobulkStatsMethod::Variance) res.var.transposeInPlace();
        num_rows = num_cols;
    }
    // Fix variance calculation to consider non-zero count
    if ((method & PseudobulkStatsMethod::Variance) == PseudobulkStatsMethod::Variance) {
        res.var = (res.var + (res.mean.square() * res.non_zeros * (((-res.non_zeros).rowwise() + group_count).rowwise() / group_count))).rowwise() / (group_count - 1);
    }
    // Fix mean calculation to consider non-zero count
    if ((method & PseudobulkStatsMethod::Mean) == PseudobulkStatsMethod::Mean) {
        res.mean = (res.mean * res.non_zeros).rowwise() / group_count;
    }
    return res;
};
template PseudobulkStats pseudobulk_matrix<uint32_t>(std::unique_ptr<MatrixLoader<uint32_t>>&& mat,
                                                     const std::vector<uint32_t>& cell_groups,
                                                     PseudobulkStatsMethod method,
                                                     bool transpose,
                                                     std::atomic<bool> *user_interrupt);
template PseudobulkStats pseudobulk_matrix<float>(std::unique_ptr<MatrixLoader<float>>&& mat,
                                                  const std::vector<uint32_t>& cell_groups,
                                                  PseudobulkStatsMethod method,
                                                  bool transpose,
                                                  std::atomic<bool> *user_interrupt);
template PseudobulkStats pseudobulk_matrix<double>(std::unique_ptr<MatrixLoader<double>>&& mat,
                                                   const std::vector<uint32_t>& cell_groups,
                                                   PseudobulkStatsMethod method,
                                                   bool transpose,
                                                   std::atomic<bool> *user_interrupt);
} // namespace BPCells