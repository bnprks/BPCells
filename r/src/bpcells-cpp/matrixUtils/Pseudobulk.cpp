// Copyright 2024 BPCells contributors
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "Pseudobulk.h"
namespace BPCells {

template <typename T>
void pseudobulk_matrix(std::unique_ptr<MatrixLoader<T>> &&mat,
                       const std::vector<uint32_t>& cell_groups,
                       PseudobulkStatsMethod method,
                       bool transpose,
                       uint32_t feature_offset,
                       PseudobulkStats &res,
                       std::atomic<bool> *user_interrupt) {
    MatrixIterator<T> it(std::move(mat));
    if (transpose && (it.rows() != cell_groups.size())) {
        throw std::invalid_argument("Pseudobulk_matrix(): Cell groups must match the number of columns in the matrix");
    } else if (!transpose && (it.cols() != cell_groups.size())) {
        throw std::invalid_argument("Pseudobulk_matrix(): Cell groups must match the number of columns in the matrix");
    }
    const bool non_zeros = (method & PseudobulkStatsMethod::NonZeros) == PseudobulkStatsMethod::NonZeros;
    const bool sum = (method & PseudobulkStatsMethod::Sum) == PseudobulkStatsMethod::Sum;
    const bool mean = (method & PseudobulkStatsMethod::Mean) == PseudobulkStatsMethod::Mean;
    const bool variance = (method & PseudobulkStatsMethod::Variance) == PseudobulkStatsMethod::Variance;

    if (!transpose) {
        // Use a one pass strategy that can do non-zero, mean, and variance at the same time
        while (it.nextCol()) {
            const uint32_t group = cell_groups[it.currentCol()];
            if (user_interrupt != NULL && *user_interrupt) return;
            while (it.load()) {
                const uint32_t *row_data = it.rowData();
                const T *val_data = it.valData();
                const uint32_t count = it.capacity();
                // Doing variance calculation also does mean and non-zero count, so specific order is required to make sure we don't double count
                if (variance) {
                    for (uint32_t i = 0; i < count; i++) update_mean_and_variance(val_data[i], res.non_zeros(feature_offset + row_data[i], group), res.mean(feature_offset + row_data[i], group), res.var(feature_offset + row_data[i], group));
                } else if (mean) {
                    for (uint32_t i = 0; i < count; i++) update_mean(val_data[i], res.non_zeros(feature_offset + row_data[i], group), res.mean(feature_offset + row_data[i], group));
                } else if (non_zeros) {
                    for (uint32_t i = 0; i < count; i++) res.non_zeros(feature_offset + row_data[i], group)++;
                }
                if (sum) {
                    for (uint32_t i = 0; i < count; i++) res.sum(feature_offset + row_data[i], group) += val_data[i];
                }
            }
        }
        return;
    }

    // Transposed: each column is a feature. Accumulate a block of features with groups as rows for
    // cache-friendliness, then copy the block into the output rows of its features.
    const uint32_t block_size = 64;
    const Eigen::Index num_groups = std::max({res.non_zeros.cols(), res.sum.cols(), res.mean.cols(), res.var.cols()});
    Eigen::ArrayXXd block_non_zeros, block_sum, block_mean, block_var;
    if (non_zeros) block_non_zeros = Eigen::ArrayXXd::Zero(num_groups, block_size);
    if (sum) block_sum = Eigen::ArrayXXd::Zero(num_groups, block_size);
    if (mean) block_mean = Eigen::ArrayXXd::Zero(num_groups, block_size);
    if (variance) block_var = Eigen::ArrayXXd::Zero(num_groups, block_size);
    uint32_t block_start = 0;
    auto copy_block = [&](Eigen::ArrayXXd &block, Eigen::Ref<Eigen::ArrayXXd> &out, uint32_t features) {
        if (block.size() == 0) return;
        out.middleRows(feature_offset + block_start, features) = block.leftCols(features).transpose();
        block.setZero();
    };
    auto flush_block = [&](uint32_t features) {
        copy_block(block_non_zeros, res.non_zeros, features);
        copy_block(block_sum, res.sum, features);
        copy_block(block_mean, res.mean, features);
        copy_block(block_var, res.var, features);
        block_start += features;
    };

    while (it.nextCol()) {
        while (it.currentCol() - block_start >= block_size) flush_block(block_size);
        const uint32_t col = it.currentCol() - block_start;
        if (user_interrupt != NULL && *user_interrupt) return;
        while (it.load()) {
            const uint32_t *row_data = it.rowData();
            const T *val_data = it.valData();
            const uint32_t count = it.capacity();
            // Doing variance calculation also does mean and non-zero count, so specific order is required to make sure we don't double count
            if (variance) {
                for (uint32_t i = 0; i < count; i++) update_mean_and_variance(val_data[i], block_non_zeros(cell_groups[row_data[i]], col), block_mean(cell_groups[row_data[i]], col), block_var(cell_groups[row_data[i]], col));
            } else if (mean) {
                for (uint32_t i = 0; i < count; i++) update_mean(val_data[i], block_non_zeros(cell_groups[row_data[i]], col), block_mean(cell_groups[row_data[i]], col));
            } else if (non_zeros) {
                for (uint32_t i = 0; i < count; i++) block_non_zeros(cell_groups[row_data[i]], col)++;
            }
            if (sum) {
                for (uint32_t i = 0; i < count; i++) block_sum(cell_groups[row_data[i]], col) += val_data[i];
            }
        }
    }
    while (block_start < it.cols()) flush_block(std::min(block_size, it.cols() - block_start));
}

void pseudobulk_finalize(PseudobulkStats &res, const std::vector<uint32_t>& cell_groups, PseudobulkStatsMethod method) {
    std::vector<double> group_count_vec;
    for (auto& cell_group : cell_groups) {
        if (cell_group >= group_count_vec.size()) group_count_vec.resize(cell_group + 1);
        group_count_vec[cell_group]++;
    }
    Eigen::Map<Eigen::Array<double, 1, Eigen::Dynamic>> group_count(group_count_vec.data(), group_count_vec.size());
    // Fix variance calculation to consider non-zero count
    if ((method & PseudobulkStatsMethod::Variance) == PseudobulkStatsMethod::Variance) {
        res.var = (res.var + (res.mean.square() * res.non_zeros * (((-res.non_zeros).rowwise() + group_count).rowwise() / group_count))).rowwise() / (group_count - 1);
    }
    // Fix mean calculation to consider non-zero count
    if ((method & PseudobulkStatsMethod::Mean) == PseudobulkStatsMethod::Mean) {
        res.mean = (res.mean * res.non_zeros).rowwise() / group_count;
    }
}

template void pseudobulk_matrix<uint32_t>(std::unique_ptr<MatrixLoader<uint32_t>>&& mat,
                                          const std::vector<uint32_t>& cell_groups,
                                          PseudobulkStatsMethod method,
                                          bool transpose,
                                          uint32_t feature_offset,
                                          PseudobulkStats &res,
                                          std::atomic<bool> *user_interrupt);
template void pseudobulk_matrix<float>(std::unique_ptr<MatrixLoader<float>>&& mat,
                                       const std::vector<uint32_t>& cell_groups,
                                       PseudobulkStatsMethod method,
                                       bool transpose,
                                       uint32_t feature_offset,
                                       PseudobulkStats &res,
                                       std::atomic<bool> *user_interrupt);
template void pseudobulk_matrix<double>(std::unique_ptr<MatrixLoader<double>>&& mat,
                                        const std::vector<uint32_t>& cell_groups,
                                        PseudobulkStatsMethod method,
                                        bool transpose,
                                        uint32_t feature_offset,
                                        PseudobulkStats &res,
                                        std::atomic<bool> *user_interrupt);
} // namespace BPCells
