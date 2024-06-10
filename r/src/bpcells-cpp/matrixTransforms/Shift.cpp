// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include <atomic>

#include "Shift.h"

namespace BPCells {

// Shift rows of a matrix
// out[i,j] = in[i,j] + row_params[0,i]
bool ShiftRows::loadZeroSubtracted(MatrixLoader<double> &loader) {
    // shift(data) - shift(0) = data by definition
    return loader.load();
}

void ShiftRows::loadZero(double *values, uint32_t count, uint32_t start_row, uint32_t col) {
    for (uint32_t i = 0; i < count; i++) {
        values[i] = fit.row_params(0, start_row + i);
    }
}

// Math tip: if A=untransformed matrix, and s = shift params as a column vector, ones = ones in a
// row vector then transform = A + s * (ones).
Eigen::MatrixXd ShiftRows::denseMultiplyRight(
    const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt
) {
    Eigen::MatrixXd res = loader->denseMultiplyRight(B, user_interrupt);
    res += fit.row_params.row(0).transpose().matrix() * B.colwise().sum();
    return res;
}
Eigen::MatrixXd ShiftRows::denseMultiplyLeft(
    const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt
) {
    Eigen::MatrixXd res = loader->denseMultiplyLeft(B, user_interrupt);
    res.colwise() += B * fit.row_params.row(0).transpose().matrix();
    return res;
}
// Calculate matrix-vector product A*v where A=this and B is a dense matrix.
Eigen::VectorXd ShiftRows::vecMultiplyRight(
    const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt
) {
    Eigen::VectorXd res = loader->vecMultiplyRight(v, user_interrupt);
    res += fit.row_params.row(0).transpose().matrix() * v.sum();
    return res;
}
Eigen::VectorXd
ShiftRows::vecMultiplyLeft(const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt) {
    Eigen::VectorXd res = loader->vecMultiplyLeft(v, user_interrupt);
    res.rowwise() += fit.row_params.row(0).matrix() * v;
    return res;
}

StatsResult
ShiftRows::computeMatrixStats(Stats row_stats, Stats col_stats, std::atomic<bool> *user_interrupt) {
    // If we shifted the rows and we haven't asked for any col_stats, then we don't have to iterate
    // through the whole dense matrix
    if (col_stats == Stats::None) {
        StatsResult res = loader->computeMatrixStats(row_stats, col_stats, user_interrupt);
        if (res.row_stats.rows() > 0) {
            res.row_stats.row(0) = cols(); // All entries count as non-zero
        }
        if (res.row_stats.rows() > 1) {
            res.row_stats.row(1) += fit.row_params.row(0); // Means offset by the row shifts
        }
        // Don't need to adjust variance at all
        return res;
    } else {
        return MatrixTransformDense::computeMatrixStats(row_stats, col_stats, user_interrupt);
    }
}

// Shift cols of a matrix
// out[i,j] = in[i,j] + col_params[0,j]
bool ShiftCols::loadZeroSubtracted(MatrixLoader<double> &loader) {
    // shift(data) - shift(0) = data by definition
    return loader.load();
}

void ShiftCols::loadZero(double *values, uint32_t count, uint32_t start_row, uint32_t col) {
    double val = fit.col_params(0, col);
    for (uint32_t i = 0; i < count; i++) {
        values[i] = val;
    }
}

// Math tip: if A=untransformed matrix, and s = shift params as a row vector, ones = ones in a col
// vector then transform = A + ones * s.
Eigen::MatrixXd ShiftCols::denseMultiplyRight(
    const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt
) {
    Eigen::MatrixXd res = loader->denseMultiplyRight(B, user_interrupt);
    res.rowwise() += fit.col_params.row(0).matrix() * B;
    return res;
}
Eigen::MatrixXd ShiftCols::denseMultiplyLeft(
    const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt
) {
    Eigen::MatrixXd res = loader->denseMultiplyLeft(B, user_interrupt);
    res += B.rowwise().sum() * fit.col_params.row(0).matrix();
    return res;
}
// Calculate matrix-vector product A*v where A=this and B is a dense matrix.
Eigen::VectorXd ShiftCols::vecMultiplyRight(
    const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt
) {
    Eigen::VectorXd res = loader->vecMultiplyRight(v, user_interrupt);
    res.rowwise() += fit.col_params.row(0).matrix() * v;
    return res;
}
Eigen::VectorXd
ShiftCols::vecMultiplyLeft(const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt) {
    Eigen::VectorXd res = loader->vecMultiplyLeft(v, user_interrupt);
    res += fit.col_params.row(0).transpose().matrix() * v.sum();
    return res;
}

StatsResult
ShiftCols::computeMatrixStats(Stats row_stats, Stats col_stats, std::atomic<bool> *user_interrupt) {
    // If we shifted the cols and we haven't asked for any row_stats, then we don't have to iterate
    // through the whole dense matrix
    if (row_stats == Stats::None) {
        StatsResult res = loader->computeMatrixStats(row_stats, col_stats, user_interrupt);
        if (res.col_stats.rows() > 0) {
            res.col_stats.row(0) = rows(); // All entries count as non-zero
        }
        if (res.col_stats.rows() > 1) {
            res.col_stats.row(1) += fit.col_params.row(0); // Means offset by the row shifts
        }
        // Don't need to adjust variance at all
        return res;
    } else {
        return MatrixTransformDense::computeMatrixStats(row_stats, col_stats, user_interrupt);
    }
}

} // end namespace BPCells
