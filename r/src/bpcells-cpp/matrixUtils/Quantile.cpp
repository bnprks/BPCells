// Copyright 2024 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "Quantile.h"
namespace BPCells {


// Find the quantile for each row in an IterableMatrix.
// Args:
// - mat: matrix to compute quantile from
// - quantile: quantile to compute from each row, between [0,1]
// Note: This function will probably crash out for large matrices
template <typename T>
std::vector<T> matrix_quantile_per_row(std::unique_ptr<MatrixLoader<T>>&& mat, 
                                       double quantile, 
                                       std::atomic<bool> *user_interrupt) {
    MatrixIterator<T> it(std::move(mat));
    std::vector<T> res;
    quantile = std::min(1.0, std::max(0.0, quantile));
    uint32_t quantile_idx = std::max(0.0, (std::ceil(quantile * it.cols()) - 1));
    std::vector<std::vector<T>> matrix(it.rows(), std::vector<T>(it.cols(), 0));
    // need to hold all the rows in memory until all columns have been iterated through
    // since matrix itself is still in CSR format
    std::vector<uint32_t> row_count(it.rows(), 0);
    while (it.nextCol()) {
        if (user_interrupt != NULL && *user_interrupt) return res;
        while (it.nextValue()) {
            matrix[it.row()][row_count[it.row()]] = it.val();
            row_count[it.row()]++;
        }
    }
    for (auto& row : matrix) {
        std::nth_element(row.begin(), row.begin() + quantile_idx, row.end());
        res.push_back(row[quantile_idx]);
    }
    return res;
};

// Find the `quantile`th value for each column in an IterableMatrix.
// Args:
// - mat: matrix to compute quantiles from
// - quantile: quantile to compute from each column, between [0,1]
template <typename T>
std::vector<T> matrix_quantile_per_col(std::unique_ptr<MatrixLoader<T>>&& mat, 
                                            double quantile, 
                                            std::atomic<bool> *user_interrupt) {
    MatrixIterator<T> it(std::move(mat));
    uint32_t curr_num;
    std::vector<T> res, curr;
    quantile = std::min(1.0, std::max(0.0, quantile));
    uint32_t quantile_idx = std::max(0.0, (std::ceil(quantile * it.rows()) - 1));
    while (it.nextCol()) {
        if (user_interrupt != NULL && *user_interrupt) return res;
        curr = std::vector<T>(it.rows(), 0);
        curr_num = 0;
        while (it.nextValue()) {
            curr[curr_num] = it.val();
            curr_num++;
        }
        std::nth_element(curr.begin(), curr.begin() + quantile_idx, curr.end());
        res.push_back(curr[quantile_idx]);
    }
    return res;
};
template std::vector<double> matrix_quantile_per_row<double>(std::unique_ptr<MatrixLoader<double>>&& mat, 
                                                             double quantile, 
                                                             std::atomic<bool> *user_interrupt);
template std::vector<float> matrix_quantile_per_row<float>(std::unique_ptr<MatrixLoader<float>>&& mat, 
                                                           double quantile, 
                                                           std::atomic<bool> *user_interrupt);
template std::vector<uint32_t> matrix_quantile_per_row<uint32_t>(std::unique_ptr<MatrixLoader<uint32_t>>&& mat,
                                                                 double quantile, 
                                                                 std::atomic<bool> *user_interrupt);
template std::vector<float> matrix_quantile_per_col<float>(std::unique_ptr<MatrixLoader<float>>&& mat, 
                                                           double quantile, 
                                                           std::atomic<bool> *user_interrupt);
template std::vector<double> matrix_quantile_per_col<double>(std::unique_ptr<MatrixLoader<double>>&& mat, 
                                                             double quantile, 
                                                             std::atomic<bool> *user_interrupt);
template std::vector<uint32_t> matrix_quantile_per_col<uint32_t>(std::unique_ptr<MatrixLoader<uint32_t>>&& mat, 
                                                             double quantile, 
                                                             std::atomic<bool> *user_interrupt);
} // namespace BPCells