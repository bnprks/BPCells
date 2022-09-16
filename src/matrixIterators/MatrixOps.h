#pragma once

#include <stdexcept>

#include "MatrixIterator.h"
#include "MatrixStats.h"

namespace BPCells {
// Calculate matrix-matrix product A*B where A (this) is sparse and B is a dense matrix.
template <typename T>
Eigen::MatrixXd MatrixLoader<T>::denseMultiplyRight(
    const Eigen::Map<Eigen::MatrixXd> B, void (*checkInterrupt)(void)
) {
    // Use transposed output so that results write contiguously in memory.
    // Note that left multiply will have better performance properties due
    // to better cache locality on the input matrix.
    // So if the performance is key and the size of the input is small,
    // it may be worth paying the price of additional transpose operations
    if (cols() != B.rows()) throw std::runtime_error("Incompatible dimensions for matrix multiply");
    Eigen::MatrixXd res(B.cols(), rows());
    res.setZero();
    restart();
    while (nextCol()) {
        const uint32_t col = currentCol();
        if (checkInterrupt != NULL && col % 128 == 0) checkInterrupt();
        while (load()) {
            const T *val_data = valData();
            const uint32_t *row_data = rowData();
            const uint32_t count = capacity();
            for (uint32_t i = 0; i < count; i++) {
                res.col(row_data[i]) += ((double)val_data[i]) * B.row(col);
            }
        }
    }
    return res.transpose();
}
template <typename T>
Eigen::MatrixXd MatrixLoader<T>::denseMultiplyLeft(
    const Eigen::Map<Eigen::MatrixXd> B, void (*checkInterrupt)(void)
) {
    if (rows() != B.cols()) throw std::runtime_error("Incompatible dimensions for matrix multiply");
    Eigen::MatrixXd res(B.rows(), cols());
    res.setZero();
    restart();
    while (nextCol()) {
        const uint32_t col = currentCol();
        if (checkInterrupt != NULL && col % 128 == 0) checkInterrupt();
        while (load()) {
            const T *val_data = valData();
            const uint32_t *row_data = rowData();
            const uint32_t count = capacity();
            for (uint32_t i = 0; i < count; i++) {
                res.col(col) += ((double)val_data[i]) * B.col(row_data[i]);
            }
        }
    }
    return res;
}
// Calculate matrix-vector product A*v where A (this) is sparse and B is a dense matrix.
template <typename T>
Eigen::VectorXd MatrixLoader<T>::vecMultiplyRight(
    const Eigen::Map<Eigen::VectorXd> v, void (*checkInterrupt)(void)
) {
    if (cols() != v.rows()) throw std::runtime_error("Incompatible dimensions for vector multiply");
    Eigen::VectorXd res(rows());
    res.setZero();
    restart();
    while (nextCol()) {
        const uint32_t col = currentCol();
        double v_col = v(col);
        if (checkInterrupt != NULL && col % 128 == 0) checkInterrupt();
        while (load()) {
            const T *val_data = valData();
            const uint32_t *row_data = rowData();
            const uint32_t count = capacity();
            // I considered manual loop unrolling here too, but it was <10% performance gain for the
            // function so I discarded that idea
            for (uint32_t i = 0; i < count; i++) {
                res(row_data[i]) += ((double)val_data[i]) * v_col;
            }
        }
    }
    return res;
}
template <typename T>
Eigen::VectorXd MatrixLoader<T>::vecMultiplyLeft(
    const Eigen::Map<Eigen::VectorXd> v, void (*checkInterrupt)(void)
) {
    if (rows() != v.rows()) throw std::runtime_error("Incompatible dimensions for vector multiply");
    Eigen::VectorXd res(cols());
    res.setZero();
    restart();
    while (nextCol()) {
        const uint32_t col = currentCol();
        if (checkInterrupt != NULL && col % 128 == 0) checkInterrupt();
        while (load()) {
            const T *val_data = valData();
            const uint32_t *row_data = rowData();
            const uint32_t count = capacity();
            // Use unrolling trick from Eigen:
            // https://gitlab.com/libeigen/eigen/-/blob/master/Eigen/src/SparseCore/SparseDenseProduct.h#L66-80
            // I tried 4x unrolling and the benefit was marginal so here we are
            uint32_t i;
            double tmp1 = 0;
            double tmp2 = 0;
            for (i = 0; i + 2 <= count; i += 2) {
                tmp1 += ((double)val_data[i]) * v(row_data[i]);
                tmp2 += ((double)val_data[i + 1]) * v(row_data[i + 1]);
            }
            for (; i < count; i++) {
                tmp1 += ((double)val_data[i]) * v(row_data[i]);
            }
            res(col) += tmp1 + tmp2;
        }
    }
    return res;
}

// Calculate row/column sums of the matrix
template <typename T> std::vector<T> MatrixLoader<T>::colSums(void (*checkInterrupt)(void)) {
    std::vector<T> sums(cols());
    restart();
    while (nextCol()) {
        const uint32_t col = currentCol();
        if (checkInterrupt != NULL && col % 128 == 0) checkInterrupt();
        while (load()) {
            const T *val_data = valData();
            const uint32_t count = capacity();
            for (uint32_t i = 0; i < count; i++) {
                sums[col] += val_data[i];
            }
        }
    }
    return sums;
}
template <typename T> std::vector<T> MatrixLoader<T>::rowSums(void (*checkInterrupt)(void)) {
    std::vector<T> sums(rows());
    restart();
    while (nextCol()) {
        const uint32_t col = currentCol();
        if (checkInterrupt != NULL && col % 128 == 0) checkInterrupt();
        while (load()) {
            const uint32_t *row_data = rowData();
            const T *val_data = valData();
            const uint32_t count = capacity();
            for (uint32_t i = 0; i < count; i++) {
                sums[row_data[i]] += val_data[i];
            }
        }
    }
    return sums;
}

// Helpers for calculating running mean + running variance
inline void update_mean(double val, double &count, double &running_mean) {
    count += 1;
    running_mean += (val - running_mean) / count;
}

inline void update_mean_and_variance(
    double val, double &count, double &running_mean, double &running_variance
) {
    count += 1;
    double delta = val - running_mean;
    running_mean += delta / count;
    double delta2 = val - running_mean;
    running_variance += delta * delta2;
}

// Calculate stats on the rows or columns of a matrix in a single pass.
// For each of rows and columns, the user can choose to from the following
// Available statistics:
// 1. None
// 2. NonZeroCount
// 3. Mean
// 4. Variance
// Each of the later statistics is always calculated simultaneously to the earlier statistics
// Outputs results to matrices row_output and col_output, which are col-major matrices,
// with one column per # rows or # columns as appropriate, and one row per output statistic
template <typename T>
StatsResult MatrixLoader<T>::computeMatrixStats(
    Stats row_stats, Stats col_stats, void (*checkInterrupt)(void)
) {
    restart();
    StatsResult res{
        Eigen::ArrayXXd((int)row_stats, rows()), Eigen::ArrayXXd((int)col_stats, cols())};
    res.row_stats.setZero();
    res.col_stats.setZero();

    while (nextCol()) {
        const uint32_t col = currentCol();
        if (checkInterrupt != NULL && col % 128 == 0) checkInterrupt();
        while (load()) {
            const uint32_t *row_data = rowData();
            const T *val_data = valData();
            const uint32_t count = capacity();

            // Update row stats if needed
            if (row_stats == Stats::NonZeroCount) {
                for (size_t i = 0; i < count; i++) {
                    res.row_stats(0, row_data[i]) += 1;
                }
            } else if (row_stats == Stats::Mean) {
                for (size_t i = 0; i < count; i++) {
                    update_mean(
                        val_data[i], res.row_stats(0, row_data[i]), res.row_stats(1, row_data[i])
                    );
                }
            } else if (row_stats == Stats::Variance) {
                for (size_t i = 0; i < count; i++) {
                    update_mean_and_variance(
                        val_data[i],
                        res.row_stats(0, row_data[i]),
                        res.row_stats(1, row_data[i]),
                        res.row_stats(2, row_data[i])
                    );
                }
            }
            // Update col stats if needed
            if (col_stats == Stats::NonZeroCount) {
                for (size_t i = 0; i < count; i++) {
                    res.col_stats(0, col) += 1;
                }
            } else if (col_stats == Stats::Mean) {
                for (size_t i = 0; i < count; i++) {
                    update_mean(val_data[i], res.col_stats(0, col), res.col_stats(1, col));
                }
            } else if (col_stats == Stats::Variance) {
                for (size_t i = 0; i < count; i++) {
                    update_mean_and_variance(
                        val_data[i],
                        res.col_stats(0, col),
                        res.col_stats(1, col),
                        res.col_stats(2, col)
                    );
                }
            }
        }
    }
    uint32_t n_rows = rows();
    uint32_t n_cols = cols();
    // fix up mean + variance accounting for the zero entries
    if (col_stats == Stats::Mean) {
        auto nnz = res.col_stats.row(0);
        auto mean = res.col_stats.row(1);
        res.col_stats.row(1) = mean * nnz / n_rows;
    } else if (col_stats == Stats::Variance) {
        // m2 += mean^2 * nnz*(nrow - nnz)/nrow;
        auto nnz = res.col_stats.row(0);
        auto mean = res.col_stats.row(1);
        auto m2 = res.col_stats.row(2);
        res.col_stats.row(2) =
            (m2 + mean.square() * nnz * ((double)n_rows - nnz) / n_rows) / (n_rows - 1);
        res.col_stats.row(1) = mean * nnz / n_rows;
    }

    if (row_stats == Stats::Mean) {
        auto nnz = res.row_stats.row(0);
        auto mean = res.row_stats.row(1);
        res.row_stats.row(1) = mean * nnz / n_cols;
    } else if (row_stats == Stats::Variance) {
        // m2 += mean^2 * nnz*(ncol - nnz)/ncol;
        auto nnz = res.row_stats.row(0);
        auto mean = res.row_stats.row(1);
        auto m2 = res.row_stats.row(2);
        res.row_stats.row(2) =
            (m2 + mean.square() * nnz * ((double)n_cols - nnz) / n_cols) / (n_cols - 1);
        res.row_stats.row(1) = mean * nnz / n_cols;
    }

    return res;
}

} // end namespace BPCells