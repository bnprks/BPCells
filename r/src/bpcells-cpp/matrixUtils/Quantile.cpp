// Copyright 2024 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "Quantile.h"
#include "../utils/radix_sort.h"
namespace BPCells {


// Find the nth order statistic given the rank, and the number of values of each sign.
// Used for quantile calculation where we do not look at zero values within a vector.
// Args:
// - rank: the rank of the order statistic to find
// - sorted_nonzero_values: Sorted non-zero values in the vector.
// - num_neg: number of negative values in the vector
// - num_zero: number of zero values in the vector
// Returns:
// - The nth order statistic in the vector.
template <typename T>
T order_statistic(const std::vector<T>& sorted_nonzero_values, uint32_t rank, uint32_t num_neg, uint32_t num_zero, uint32_t num_pos) {
    if (rank < num_neg) {
        return sorted_nonzero_values[rank];
    } else if (rank < num_neg + num_zero) {
        return 0;
    } else {
        return sorted_nonzero_values[rank - num_zero];
    }
}

// Find the `quantile`th value for each column in an IterableMatrix.
// Please refer to the `Statistics.quantile` function in `julialang` for more information on how quantiles are calculated.
// Args:
// - mat: matrix to compute quantiles from
// - quantile: quantile to compute from each column, between [0,1]
// - alpha: parameter for the quantile calculation
// - beta: parameter for the quantile calculation
// Returns:
// - A vector of quantile values, one for each column.
template <typename T>
std::vector<T> matrix_quantile_per_col(std::unique_ptr<MatrixLoader<T>>&& mat, 
                                       double quantile,
                                       double alpha,
                                       double beta,
                                       std::atomic<bool> *user_interrupt) {
    MatrixIterator<T> it(std::move(mat));
    // keep track of the current number of values in a column
    uint32_t curr_num;
    std::vector<T> res, curr;
    // clamp quantile to [0,1]
    quantile = std::min(1.0, std::max(0.0, quantile));
    
    // clamp alpha and beta to [0,1]
    alpha = std::min(1.0, std::max(0.0, alpha));
    beta = std::min(1.0, std::max(0.0, beta));
    // stats used for quantile calculation
    double m = alpha + quantile * (1.0-alpha-beta);
    double aleph = (quantile * it.rows() +m);
    // index represents one higher than the true index
    uint32_t index = static_cast<uint32_t>(std::floor(std::max(std::min(it.rows()-1.0, aleph), 1.0)));
    double gamma = std::max(std::min(aleph - index, 1.0), 0.0);

    // do not want to consider zero values in sorting/quantile, so we keep track of the number of negative values
    uint32_t num_neg = 0;
    std::vector<T> buffer;
    double col_num;
    while (it.nextCol()) {
        num_neg = 0;
        if (user_interrupt != NULL && *user_interrupt) return res;
        curr.clear();
        curr_num = 0;
        while (it.nextValue()) {
            if (it.val() < 0) num_neg++;
            curr.push_back(it.val());
            curr_num++;
        }
        uint32_t num_zeros = it.rows() - curr_num;
        buffer.resize(curr.size());
        BPCells::lsdRadixSortArrays(curr.size(), curr, buffer);
        if (index >= it.rows()) {
            col_num = curr[curr.size() - 1];
        } else {
            col_num = 
                order_statistic(curr, index - 1, num_neg, num_zeros, it.rows() - num_neg - num_zeros)* (1-gamma) +
                order_statistic(curr, index, num_neg, num_zeros, it.rows() - num_neg - num_zeros) * gamma;
        }
        res.push_back(col_num);
    }
    return res;
};

template std::vector<float> matrix_quantile_per_col<float>(std::unique_ptr<MatrixLoader<float>>&& mat,
                                                           double quantile,
                                                           double alpha,
                                                           double beta,
                                                           std::atomic<bool> *user_interrupt);
template std::vector<double> matrix_quantile_per_col<double>(std::unique_ptr<MatrixLoader<double>>&& mat, 
                                                             double quantile,
                                                             double alpha,
                                                             double beta,
                                                             std::atomic<bool> *user_interrupt);
template std::vector<uint32_t> matrix_quantile_per_col<uint32_t>(std::unique_ptr<MatrixLoader<uint32_t>>&& mat, 
                                                                 double quantile, 
                                                                 double alpha,
                                                                 double beta,
                                                                 std::atomic<bool> *user_interrupt);
} // namespace BPCells