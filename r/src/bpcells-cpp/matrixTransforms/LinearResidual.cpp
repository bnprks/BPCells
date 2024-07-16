// Copyright 2024 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "LinearResidual.h"
#include <atomic>

namespace BPCells {

bool LinearResidual::loadZeroSubtracted(MatrixLoader<double> &loader) {
    return loader.load();
}

void LinearResidual::loadZero(
    double *values, uint32_t count, uint32_t start_row, uint32_t col
) {
    Eigen::Map<Eigen::VectorXd> values_vec(values, count);
    // For prediction_axis == "row" (default)
    // col_params = t(Q)
    // row_params = Qty
    // val = X - t(Qty) %*% t(Q)
    // val = X - t(row_params) %*% col_params
    values_vec = -fit.row_params.middleCols(start_row, count).matrix().transpose() * 
      fit.col_params.col(col).matrix();
}

Eigen::MatrixXd LinearResidual::denseMultiplyRight(
    const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt
) {
    Eigen::MatrixXd res = loader->denseMultiplyRight(B, user_interrupt);
    res -= fit.row_params.matrix().transpose() * (fit.col_params.matrix() * B);
    return res;
}
Eigen::MatrixXd LinearResidual::denseMultiplyLeft(
    const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt
) {
    Eigen::MatrixXd res = loader->denseMultiplyLeft(B, user_interrupt);
    res -= (B * fit.row_params.matrix().transpose()) * fit.col_params.matrix();
    return res;
}
Eigen::VectorXd LinearResidual::vecMultiplyRight(
    const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt
) {
    Eigen::VectorXd res = loader->vecMultiplyRight(v, user_interrupt);
    res -= fit.row_params.matrix().transpose() * (fit.col_params.matrix() * v);
    return res;
}
Eigen::VectorXd LinearResidual::vecMultiplyLeft(
    const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt
) {
    Eigen::VectorXd res = loader->vecMultiplyLeft(v, user_interrupt);
    res -= (v.transpose() * fit.row_params.matrix().transpose()) * fit.col_params.matrix();
    return res;
}

} // namespace BPCells
