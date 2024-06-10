// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include <atomic>
#include "MatrixTransform.h"

namespace BPCells {

MatrixTransform::MatrixTransform(std::unique_ptr<MatrixLoader<double>> &&loader)
    : MatrixLoaderWrapper<double>(std::move(loader)) {}
MatrixTransform::MatrixTransform(std::unique_ptr<MatrixLoader<double>> &&loader, TransformFit fit)
    : MatrixLoaderWrapper<double>(std::move(loader))
    , fit(fit) {
    // Basic checks for dimensions
    if (fit.row_params.cols() != 0 && fit.row_params.cols() != this->loader->rows()) {
        throw std::runtime_error("C++ error constructing MatrixTransform: fit.row_params.cols() != loader->rows()");
    }
    if (fit.col_params.cols() != 0 && fit.col_params.cols() != this->loader->cols()) {
        throw std::runtime_error("C++ error constructing MatrixTransform: fit.col_params.cols() != loader->cols()");
    }
}

TransformFit MatrixTransform::getFit() { return fit; }

MatrixTransformDense::MatrixTransformDense(std::unique_ptr<MatrixLoader<double>> &&mat, TransformFit fit)
    : MatrixTransform(std::make_unique<OrderRows<double>>(std::move(mat)), fit) {}

void MatrixTransformDense::restart() {
    loader->restart();
    loader_idx = UINT32_MAX;
    loader_col = UINT32_MAX;
    loader_capacity = 0;
    current_col = UINT32_MAX;
    current_row = UINT32_MAX;
}

void MatrixTransformDense::seekCol(uint32_t col) {
    current_col = col;
    current_row = 0;
    loader->seekCol(col);
    loader_col = loader->currentCol();
    loader_capacity = UINT32_MAX;
    loader_idx = UINT32_MAX;
}

bool MatrixTransformDense::nextCol() {
    loader_capacity = 0;
    loader_idx = 0;
    current_row = 0;
    if (current_col == loader_col) {
        if (loader->nextCol()) {
            loader_col = loader->currentCol();
            loader_idx = UINT32_MAX;
            loader_capacity = UINT32_MAX;
        } else {
            loader_col = UINT32_MAX;
        }
    }
    current_col += 1;
    if (current_col >= loader->cols()) {
        current_col = loader->cols();
        current_row = UINT32_MAX;
        return false;
    }
    return true;
}

uint32_t MatrixTransformDense::currentCol() const { return current_col; }

bool MatrixTransformDense::load() {
    if (current_row >= loader->rows()) return false;
    uint32_t load_size = std::min(buf_size, loader->rows() - current_row);

    // Load the dense values assuming the underlying data is all 0
    for (uint32_t i = 0; i < load_size; i++)
        row_data[i] = current_row + i;
    loadZero(val_data.data(), load_size, current_row, current_col);

    // Correct the values at entries that are not actually 0
    uint32_t *loader_row = loader->rowData();
    double *loader_val = loader->valData();
    while (loader_capacity > 0) {
        if (loader_idx >= loader_capacity) {
            if (!loadZeroSubtracted(*loader)) loader_capacity = 0;
            else loader_capacity = loader->capacity();
            loader_idx = 0;
            if (loader_capacity == 0) break;
            loader_row = loader->rowData();
            loader_val = loader->valData();
        }
        if (loader_row[loader_idx] >= current_row + load_size) break;
        val_data[loader_row[loader_idx] - current_row] += loader_val[loader_idx];
        loader_idx += 1;
    }
    current_row += load_size;
    return true;
}

uint32_t MatrixTransformDense::capacity() const {
    // Always buf_size capacity unless we're at the end of a column, in which case whatever the
    // remainder is of buf_size
    return current_row != loader->rows() ? buf_size
                                        : current_row - ((current_row - 1) / buf_size) * buf_size;
}
uint32_t *MatrixTransformDense::rowData() { return row_data.data(); }
double *MatrixTransformDense::valData() { return val_data.data(); }

Eigen::MatrixXd MatrixTransformDense::denseMultiplyRight(
    const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt
) {
    // Perform the denseMultiplyRight operation like in MatrixOps.h, but using
    // loadZeroSubtracted in place of load
    if (cols() != B.rows()) throw std::runtime_error("Incompatible dimensions for matrix multiply");
    Eigen::MatrixXd res(B.cols(), rows());
    res.setZero();
    MatrixLoader<double> &unordered_loader = *(dynamic_cast<OrderRows<double> &>(*(this->loader.get())).loader.get());
    restart();
    while (nextCol()) {
        const uint32_t col = currentCol();
        if (user_interrupt != NULL && *user_interrupt) return res;
        // Don't need ordered loads here
        while (loadZeroSubtracted(unordered_loader)) {
            const double *val_data = unordered_loader.valData();
            const uint32_t *row_data = unordered_loader.rowData();
            const uint32_t count = unordered_loader.capacity();
            for (uint32_t i = 0; i < count; i++) {
                res.col(row_data[i]) += ((double)val_data[i]) * B.row(col);
            }
        }
    }

    // Make adjustments for the zero entries
    res.transposeInPlace();
    denseMultiplyRightZero(res, B, user_interrupt);
    return res;
}

Eigen::MatrixXd MatrixTransformDense::denseMultiplyLeft(
    const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt
) {
    // Perform the denseMultiplyRight operation like in MatrixOps.h, but using
    // loadZeroSubtracted in place of load
    if (rows() != B.cols()) throw std::runtime_error("Incompatible dimensions for matrix multiply");
    Eigen::MatrixXd res(B.rows(), cols());
    res.setZero();
    MatrixLoader<double> &unordered_loader = *(dynamic_cast<OrderRows<double> &>(*(this->loader.get())).loader.get());
    restart();
    while (nextCol()) {
        const uint32_t col = currentCol();
        if (user_interrupt != NULL && *user_interrupt) return res;
        // Don't need ordered loads here
        while (loadZeroSubtracted(unordered_loader)) {
            const double *val_data = unordered_loader.valData();
            const uint32_t *row_data = unordered_loader.rowData();
            const uint32_t count = unordered_loader.capacity();
            for (uint32_t i = 0; i < count; i++) {
                res.col(col) += ((double)val_data[i]) * B.col(row_data[i]);
            }
        }
    }

    // Make adjustments for the zero entries
    denseMultiplyLeftZero(res, B, user_interrupt);
    return res;
}

// Calculate matrix-vector product A*v where A (this) is sparse and B is a dense matrix.
Eigen::VectorXd MatrixTransformDense::vecMultiplyRight(
    const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt
) {
    // Perform the vecMultiplyRight operation like in MatrixOps.h, but using
    // loadZeroSubtracted in place of load
    if (cols() != v.rows()) throw std::runtime_error("Incompatible dimensions for vector multiply");
    Eigen::VectorXd res(rows());
    res.setZero();
    MatrixLoader<double> &unordered_loader = *(dynamic_cast<OrderRows<double> &>(*(this->loader.get())).loader.get());
    restart();
    while (nextCol()) {
        const uint32_t col = currentCol();
        if (user_interrupt != NULL && *user_interrupt) return res;
        // Don't need ordered loads here
        while (loadZeroSubtracted(unordered_loader)) {
            const double *val_data = unordered_loader.valData();
            const uint32_t *row_data = unordered_loader.rowData();
            const uint32_t count = unordered_loader.capacity();
            for (uint32_t i = 0; i < count; i++) {
                res(row_data[i]) += ((double)val_data[i]) * v(col);
            }
        }
    }

    // Make adjustments for the zero entries
    vecMultiplyRightZero(res, v, user_interrupt);
    return res;
}

Eigen::VectorXd MatrixTransformDense::vecMultiplyLeft(
    const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt
) {
    // Perform the vecMultiplyLeft operation like in MatrixOps.h, but using
    // loadZeroSubtracted in place of load
    if (rows() != v.rows()) throw std::runtime_error("Incompatible dimensions for vector multiply");
    Eigen::VectorXd res(cols());
    res.setZero();
    MatrixLoader<double> &unordered_loader = *(dynamic_cast<OrderRows<double> &>(*(this->loader.get())).loader.get());
    restart();
    while (nextCol()) {
        const uint32_t col = currentCol();
        if (user_interrupt != NULL && *user_interrupt) return res;
        // Don't need ordered loads here
        while (loadZeroSubtracted(unordered_loader)) {
            const double *val_data = unordered_loader.valData();
            const uint32_t *row_data = unordered_loader.rowData();
            const uint32_t count = unordered_loader.capacity();
            for (uint32_t i = 0; i < count; i++) {
                res(col) += ((double)val_data[i]) * v(row_data[i]);
            }
        }
    }

    vecMultiplyLeftZero(res, v, user_interrupt);
    return res;
}

void MatrixTransformDense::denseMultiplyRightZero(
    Eigen::MatrixXd &out, const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt
) {
    // Key invariants for this block processing: for L*R = O: colL=rowR, rowO=rowL, colO=colR
    Eigen::Matrix<double, buf_size, 1> values;
    uint32_t nrows = rows();
    uint32_t ncols = cols();
    restart();

    for (uint32_t col = 0; col < ncols; col++) {
        if (user_interrupt != NULL && *user_interrupt) return;
        uint32_t row;
        for (row = 0; row + buf_size <= nrows; row += buf_size) {
            loadZero(values.data(), buf_size, row, col);
            out.middleRows<buf_size>(row) += values * B.row(col);
        }
        if (row < nrows) {
            loadZero(values.data(), nrows - row, row, col);
            out.middleRows(row, nrows - row) += values.topRows(nrows - row) * B.row(col);
        }
    }
}

void MatrixTransformDense::denseMultiplyLeftZero(
    Eigen::MatrixXd &out, const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt
) {
    // Key invariants for this block processing: for L*R = O: colL=rowR, rowO=rowL, colO=colR
    Eigen::Matrix<double, buf_size, 1> values;
    uint32_t nrows = rows();
    uint32_t ncols = cols();
    restart();

    for (uint32_t col = 0; col < ncols; col++) {
        if (user_interrupt != NULL && *user_interrupt) return;
        uint32_t row;
        for (row = 0; row + buf_size <= nrows; row += buf_size) {
            loadZero(values.data(), buf_size, row, col);
            out.col(col) += B.middleCols<buf_size>(row) * values;
        }
        if (row < nrows) {
            loadZero(values.data(), nrows - row, row, col);
            out.col(col) += B.middleCols(row, nrows - row) * values.topRows(nrows - row);
        }
    }
}

void MatrixTransformDense::vecMultiplyRightZero(
    Eigen::VectorXd &out, const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt
) {
    // Key invariants for this block processing: for L*R = O: colL=rowR, rowO=rowL, colO=colR
    Eigen::Matrix<double, buf_size, 1> values;
    uint32_t nrows = rows();
    uint32_t ncols = cols();
    restart();

    for (uint32_t col = 0; col < ncols; col++) {
        if (user_interrupt != NULL && *user_interrupt) return;
        uint32_t row;
        for (row = 0; row + buf_size <= nrows; row += buf_size) {
            loadZero(values.data(), buf_size, row, col);
            out.middleRows<buf_size>(row) += values * v(col);
        }
        if (row < nrows) {
            loadZero(values.data(), nrows - row, row, col);
            out.middleRows(row, nrows - row) += values.topRows(nrows - row) * v(col);
        }
    }
}

void MatrixTransformDense::vecMultiplyLeftZero(
    Eigen::VectorXd &out, const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt
) {
    Eigen::Matrix<double, buf_size, 1> values;
    uint32_t nrows = rows();
    uint32_t ncols = cols();
    restart();

    for (uint32_t col = 0; col < ncols; col++) {
        if (user_interrupt != NULL && *user_interrupt) return;
        uint32_t row;
        for (row = 0; row + buf_size <= nrows; row += buf_size) {
            loadZero(values.data(), buf_size, row, col);
            out(col) += v.middleRows<buf_size>(row).transpose() * values;
        }
        if (row < nrows) {
            loadZero(values.data(), nrows - row, row, col);
            out(col) += v.middleRows(row, nrows - row).transpose() * values.topRows(nrows - row);
        }
    }
}

std::vector<double> MatrixTransformDense::colSums(std::atomic<bool> *user_interrupt) {
    std::vector<double> out(cols());

    Eigen::VectorXd v(rows());
    v.setOnes();
    Eigen::VectorXd res =
        vecMultiplyLeft(Eigen::Map<Eigen::VectorXd>(v.data(), v.rows(), v.cols()), user_interrupt);
    for (uint32_t i = 0; i < out.size(); i++) {
        out[i] = res(i);
    }
    return out;
}
std::vector<double> MatrixTransformDense::rowSums(std::atomic<bool> *user_interrupt) {
    std::vector<double> out(rows());

    Eigen::VectorXd v(cols());
    v.setOnes();
    Eigen::VectorXd res =
        vecMultiplyRight(Eigen::Map<Eigen::VectorXd>(v.data(), v.rows(), v.cols()), user_interrupt);
    for (uint32_t i = 0; i < out.size(); i++) {
        out[i] = res(i);
    }
    return out;
}

} // end namespace BPCells
