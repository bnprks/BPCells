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
T order_statistic(const std::vector<T>& sorted_nonzero_values, 
                  uint32_t rank, 
                  uint32_t num_neg, 
                  uint32_t num_zero, 
                  uint32_t num_pos);


// Find the `quantile`th value(s) for each column in an IterableMatrix.
// Please refer to the `Statistics.quantile` function in `julialang` for more information on how quantiles are calculated. https://docs.julialang.org/en/v1/stdlib/Statistics/#Statistics.quantile
// Args:
// - mat: matrix to compute quantiles from
// - quantile: quantile to compute from each column, between [0,1]
// - alpha: parameter for the quantile calculation
// - beta: parameter for the quantile calculation
// Returns:
// - A matrix of quantiles values, with each column corresponding to a quantile and each row corresponding to a column in the input matrix.
template <typename T>
Eigen::ArrayXXd matrix_quantile_per_col(std::unique_ptr<MatrixLoader<T>>&& mat, 
                                        std::vector<double> quantile,
                                        double alpha,
                                        double beta,
                                        std::atomic<bool> *user_interrupt);
} // namespace BPCells