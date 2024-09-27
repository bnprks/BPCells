// Copyright 2024 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "Quantile.h"
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
                                       std::atomic<bool> *user_interrupt) {
    MatrixIterator<T> it(std::move(mat));
    uint32_t curr_num;
    std::vector<T> res, curr;
    quantile = std::min(1.0, std::max(0.0, quantile));
    // uint32_t quantile_idx = std::max(0.0, (std::ceil(quantile * it.rows()) - 1));
    uint32_t num_neg = 0;
    while (it.nextCol()) {
        if (user_interrupt != NULL && *user_interrupt) return res;
        curr.clear();
        curr_num = 0;
        while (it.nextValue()) {
            if (it.val() < 0) num_neg++;
            curr.push_back(it.val());
            curr_num++;
        }
        uint32_t num_zeros = it.rows() - curr_num;
        // calculate using type 7 quantile
        double pos = quantile * (it.rows() - 1);
        uint32_t index = static_cast<uint32_t>(std::floor(pos));
        double h = pos - index;
        std::sort(curr.begin(), curr.end());
        double col_num;
        // type 7 quantiles requires a little bit of extra logic to handle edge cases when there are zeros
        if (index < num_neg) {
            if (index == num_neg - 1) {
                if (num_zeros > 0) {
                    // at the boundary between non-zeros and zeros
                    col_num = curr[index] * (1-h);
                } else if (index < it.rows() - 1) {
                    // there are no zeroes, and index is the last negative value
                    col_num = curr[index] * (1-h) + curr[index + 1] * h;
                } else {
                    //there are no zeroes, and index is the last value
                    col_num = curr[index];
                }
            } else {
                // regular case for calculating quantiles
                col_num = curr[index] * (1-h) + curr[index + 1] * h;
            }
        } else if (index < num_neg + num_zeros) {
            if (index == num_neg + num_zeros - 1) {
                if (index < it.rows() - 1) {
                    // at the boundary between zeros and positive numbers
                    col_num = curr[index - num_zeros + 1] * h;
                } else {
                    // at the last value
                    col_num = 0;
                }
            } else {
                // in the middle of the zeros
                col_num = 0;
            }
        } else {
            if (index == it.rows() - 1) {
                // at the last positive value
                col_num = curr[-1];
            } else {
                // regular case for calculating quantiles
                col_num = curr[index - num_zeros] * (1-h) + curr[index - num_zeros  + 1] * h;
            }
        }
        res.push_back(col_num);
    }
    return res;
};

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