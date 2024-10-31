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


// Find the `n`th order statistic given the rank, and the number of values of each sign in a sorted vector.
// Used for quantile calculation where we do not include zero values within the sorted vector.
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

// Find the `quantile`th value(s) for each column in an IterableMatrix.
// Please refer to the `Statistics.quantile` function in `julialang` for more information on how quantiles are calculated.
// Args:
// - mat: matrix to compute quantiles from
// - quantile: quantile to compute from each column, between [0,1]
// - alpha: parameter for the quantile calculation
// - beta: parameter for the quantile calculation
// Returns:
// - A matrix of quantiles values, with each row corresponding to a quantile and each column corresponding to a column in the matrix.
template <typename T>
Eigen::ArrayXXd matrix_quantile_per_col(std::unique_ptr<MatrixLoader<T>>&& mat, 
                                        std::vector<double> quantile,
                                        double alpha,
                                        double beta,
                                        std::atomic<bool> *user_interrupt) {
    MatrixIterator<T> it(std::move(mat));
    Eigen::ArrayXXd res(it.cols(), quantile.size());
    std::vector<T> curr;
    // clamp quantile, alpha, beta to [0,1]
    for (auto& q: quantile) q = std::min(1.0, std::max(0.0, q));
    alpha = std::min(1.0, std::max(0.0, alpha));
    beta = std::min(1.0, std::max(0.0, beta));
    // stats used for quantile index calculation
    // store the index and gamma for each quantile
    std::vector<uint32_t> indexes;
    std::vector<double> gammas;
    for (uint32_t q_idx = 0; q_idx < quantile.size(); q_idx++) {
        double m = alpha + quantile[q_idx] * (1.0-alpha-beta);
        double aleph = (quantile[q_idx] * it.rows() + m);
        // index represents one higher than the true index
        indexes.push_back(static_cast<uint32_t>(std::floor(std::max(std::min(it.rows()-1.0, aleph), 1.0))));
        gammas.push_back(std::max(std::min(aleph - indexes[q_idx], 1.0), 0.0));
    }

    // do not want to consider zero values in sorting, so we keep track of the number of negative values
    uint32_t num_neg = 0;
    std::vector<T> buffer;
    while (it.nextCol()) {
        num_neg = 0;
        if (user_interrupt != NULL && *user_interrupt) return res;
        curr.clear();
        while (it.nextValue()) {
            if (it.val() < 0) num_neg++;
            curr.push_back(it.val());
        }
        uint32_t num_zeros = it.rows() - curr.size();
        buffer.resize(curr.size());
        BPCells::lsdRadixSortArrays(curr.size(), curr, buffer);
        for (uint32_t q_idx = 0; q_idx < quantile.size(); q_idx++) {
            double quantile_num = order_statistic(curr, indexes[q_idx] - 1, num_neg, num_zeros, it.rows() - num_neg - num_zeros)* (1-gammas[q_idx]) +
                order_statistic(curr, indexes[q_idx], num_neg, num_zeros, it.rows() - num_neg - num_zeros) * gammas[q_idx];
            res(it.currentCol(), q_idx) = quantile_num;
        }
    }
    return res;
};

template Eigen::ArrayXXd matrix_quantile_per_col<float>(std::unique_ptr<MatrixLoader<float>>&& mat,
                                                        std::vector<double> quantile,
                                                        double alpha,
                                                        double beta,
                                                        std::atomic<bool> *user_interrupt);
template Eigen::ArrayXXd matrix_quantile_per_col<double>(std::unique_ptr<MatrixLoader<double>>&& mat, 
                                                         std::vector<double> quantile,
                                                         double alpha,
                                                         double beta,
                                                         std::atomic<bool> *user_interrupt);
template Eigen::ArrayXXd matrix_quantile_per_col<uint32_t>(std::unique_ptr<MatrixLoader<uint32_t>>&& mat, 
                                                           std::vector<double> quantile, 
                                                           double alpha,
                                                           double beta,
                                                           std::atomic<bool> *user_interrupt);
} // namespace BPCells