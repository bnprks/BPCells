// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include <cmath>
#include <stdexcept>

#include "../matrixIterators/ColwiseRank.h"
#include "WilcoxonRankSum.h"
namespace BPCells {

// Perform Wilcoxon rank sum test (aka Mann-Whitney U test) on each column of a matrix.
// Given the groupings specified by `groups`, perform a 1-vs-rest test, and return a matrix
// of p-values with dimensions (# groups) x (# columns)
template <typename T>
Eigen::MatrixXd
wilcoxon_rank_sum(std::unique_ptr<MatrixLoader<T>> &&mat, const std::vector<uint32_t> &groups, std::atomic<bool> *user_interrupt) {
    if (groups.size() != mat->rows()) {
        throw std::runtime_error("Error in wilcoxon_rank_sum: groups length != mat.rows()");
    }
    ColwiseRank<T> ranks(std::move(mat));

    uint32_t n_groups = 0;
    for (auto x : groups) {
        n_groups = std::max(n_groups, x + 1);
    }

    Eigen::VectorXd group_size(n_groups);
    Eigen::VectorXd rank_sum(n_groups);
    Eigen::MatrixXd pval(n_groups, ranks.cols());

    const double N = groups.size();
    

    group_size.setZero();
    for (auto x : groups) {
        group_size[x]++;
    }

    while (ranks.nextCol()) {
        if (user_interrupt != NULL && *user_interrupt) break;
        uint32_t col = ranks.currentCol();
        rank_sum.setZero();
        double total_rank = 0;
        while (ranks.load()) {
            uint32_t cap = ranks.capacity();
            double *rank_data = ranks.valData();
            uint32_t *row_data = ranks.rowData();
            for (uint32_t i = 0; i < cap; i++) {
                rank_sum[groups[row_data[i]]] += rank_data[i];
                total_rank += rank_data[i];
            }
        }
        for (size_t i = 0; i < n_groups; i++) {
            // Variable definisions as per wikipedia page formulas
            const double n1 = group_size[i];
            const double n2 = N - group_size[i];

            // Test statistic
            const double rank_offset = (N * (N + 1) / 2 - total_rank) / N;
            const double U1 = rank_sum[i] + rank_offset * n1 - n1 * (n1 + 1) / 2;
            const double U2 = total_rank - rank_sum[i] + rank_offset * n2 - n2 * (n2 + 1) / 2;
            const double U = std::max(U1, U2);
            // Calculate significance value from normal distribution (see wikipedia article)
            // To calculate U_mean accounting for the offset
            const double U_mean = n1 * n2 / 2;

            const double U_std_dev =
                std::sqrt(n1 * n2 / 12 * (N + 1 - ranks.tieStatistic() / (N * (N - 1))));

            // The -0.5 in here is for the continuity correction, as mentioned in R's
            // wilcox.test documentation. I believe it is also discussed in Lehmann's book
            // Nonparametrics: statistical methods based on ranks, though I have been
            // unable to view a clean copy of it myself
            const double continuity_correction = (U > U_mean) * 0.5;
            const double z_score = (U - continuity_correction - U_mean)  / U_std_dev;
            // If all values are tied (likely 0), then U_std_dev is 0
            const double p = U_std_dev != 0 ? std::erfc(z_score / std::sqrt(2)) : 1;
            pval(i, col) = p;
        }
    }
    return pval;
}

template Eigen::MatrixXd wilcoxon_rank_sum<uint32_t>(std::unique_ptr<MatrixLoader<uint32_t>> &&mat, const std::vector<uint32_t> &groups, std::atomic<bool> *user_interrupt);
template Eigen::MatrixXd wilcoxon_rank_sum<float>(std::unique_ptr<MatrixLoader<float>> &&mat, const std::vector<uint32_t> &groups, std::atomic<bool> *user_interrupt);
template Eigen::MatrixXd wilcoxon_rank_sum<uint64_t>(std::unique_ptr<MatrixLoader<uint64_t>> &&mat, const std::vector<uint32_t> &groups, std::atomic<bool> *user_interrupt);
template Eigen::MatrixXd wilcoxon_rank_sum<double>(std::unique_ptr<MatrixLoader<double>> &&mat, const std::vector<uint32_t> &groups, std::atomic<bool> *user_interrupt);
} // namespace BPCells
