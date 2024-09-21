// Copyright 2024 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "Pseudobulk.h"
namespace BPCells {


// Calculate variance of the values in a matrix, for each pseudobulk.  Requires the means of each pseudobulk to be passed in.
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
                                  std::atomic<bool> *user_interrupt) {
    MatrixIterator<T> it(std::move(mat));
    // K = group name, V = count of group
    std::unordered_map<std::string, double> group_count;
    // K = group name, V = group index
    std::unordered_map<std::string, uint32_t> group_map;
    struct PseudobulkStats res;
    uint32_t group_num = 0;
    std::vector<std::string> group_names;
    for (auto& cell_group : cell_groups) {
        if (group_map.find(std::to_string(cell_group)) == group_map.end()) {
            group_map[std::to_string(cell_group)] = group_num;
            group_num++;
            group_names.push_back(std::to_string(cell_group));
            group_count[std::to_string(cell_group)] = 1;
        } else {
            group_count[std::to_string(cell_group)]++;
        }
    }
    uint32_t num_rows = transpose ? it.cols() : it.rows();
    res.non_zeros = Eigen::ArrayXXd::Zero(num_rows, group_num);
    if (method != "non-zeros") res.sum = Eigen::ArrayXXd::Zero(num_rows, group_num);
    if (method == "var" || method == "mean") res.mean = Eigen::ArrayXXd::Zero(num_rows, group_num);
    if (method == "var") res.var = Eigen::ArrayXXd::Zero(num_rows, group_num);
    // Use a one pass strategy that can do non-zero, mean, and variance at the same time
    while (it.nextCol()) {
        const uint32_t col = it.currentCol();
        if (user_interrupt != NULL && *user_interrupt) return res;
        while (it.load()) {
            const uint32_t *row_data = it.rowData();
            const T *val_data = it.valData();
            const uint32_t count = it.capacity();
            for (uint32_t i = 0; i < count; i++) {
                uint32_t row_idx, col_idx;
                if (transpose) {
                    row_idx = col;
                    col_idx = group_map.at(std::to_string(cell_groups[row_data[i]]));
                } else {
                    row_idx = row_data[i];
                    col_idx = group_map.at(std::to_string(cell_groups[col]));
                }
                if (method == "sum") {
                    res.sum(row_idx, col_idx) += val_data[i];
                }
                // Each method below requires the previous summary statistic to be calculated
                // in the order of non-zeroes, mean, variance
                res.non_zeros(row_idx, col_idx)++;
                // Calculate using a running mean using formula `running_mean = running_mean + (val - running_mean) / n`
                if (method == "mean") {
                    res.mean(row_idx, col_idx) += (val_data[i] - res.mean(row_idx, col_idx)) / res.non_zeros(row_idx, col_idx);
                } else if (method == "var") {
                    double delta = val_data[i] - res.mean(row_idx, col_idx);
                    res.mean(row_idx, col_idx) += delta / res.non_zeros(row_idx, col_idx);
                    double delta2 = val_data[i] - res.mean(row_idx, col_idx);
                    res.var(row_idx, col_idx) += delta * delta2;
                }
            }
        }
    }
    if (method == "mean" || method == "var") {
        // fix up mean for non zeros using mean * num_non_zeros / num_total
        for (int group_idx = 0; group_idx < res.mean.cols(); group_idx++) {
            for (int row_idx = 0; row_idx < res.mean.rows(); row_idx++) {
                // any way of element wise vectorizing this without instantiating a matrix for group_totals?
                if (method == "var") {
                    res.var(row_idx, group_idx) = 
                        (res.var(row_idx, group_idx) + 
                        std::pow(res.mean(row_idx, group_idx), 2) *
                        res.non_zeros(row_idx, group_idx) * 
                        ((group_count.at(group_names[group_idx]) - res.non_zeros(row_idx, group_idx)) / group_count.at(group_names[group_idx]))) / 
                        (group_count.at(group_names[group_idx]) - 1);
                }
                res.mean(row_idx, group_idx) = res.mean(row_idx, group_idx) * res.non_zeros(row_idx, group_idx) / group_count.at(group_names[group_idx]);
            }
        }
    }
    return res;
};
template PseudobulkStats pseudobulk_matrix<uint32_t>(std::unique_ptr<MatrixLoader<uint32_t>>&& mat,
                                                     const std::vector<uint32_t>& cell_groups,
                                                     const std::string& method,
                                                     bool transpose,
                                                     std::atomic<bool> *user_interrupt);
template PseudobulkStats pseudobulk_matrix<float>(std::unique_ptr<MatrixLoader<float>>&& mat,
                                                  const std::vector<uint32_t>& cell_groups,
                                                  const std::string& method,
                                                  bool transpose,
                                                  std::atomic<bool> *user_interrupt);
template PseudobulkStats pseudobulk_matrix<double>(std::unique_ptr<MatrixLoader<double>>&& mat,
                                                   const std::vector<uint32_t>& cell_groups,
                                                   const std::string& method,
                                                   bool transpose,
                                                   std::atomic<bool> *user_interrupt);
} // namespace BPCells