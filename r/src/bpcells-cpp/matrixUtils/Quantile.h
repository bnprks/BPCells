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
T order_statistic(const std::vector<T>& sorted_nonzero_values, 
                  uint32_t rank, 
                  uint32_t num_neg, 
                  uint32_t num_zero, 
                  uint32_t num_pos);


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
                                       std::atomic<bool> *user_interrupt);
} // namespace BPCells