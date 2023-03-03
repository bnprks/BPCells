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
Eigen::MatrixXd
ShiftRows::denseMultiplyRight(const Eigen::Map<Eigen::MatrixXd> B, void (*checkInterrupt)(void)) {
    Eigen::MatrixXd res = loader->denseMultiplyRight(B, checkInterrupt);
    res += fit.row_params.row(0).transpose().matrix() * B.colwise().sum();
    return res;
}
Eigen::MatrixXd
ShiftRows::denseMultiplyLeft(const Eigen::Map<Eigen::MatrixXd> B, void (*checkInterrupt)(void)) {
    Eigen::MatrixXd res = loader->denseMultiplyLeft(B, checkInterrupt);
    res.colwise() += B * fit.row_params.row(0).transpose().matrix();
    return res;
}
// Calculate matrix-vector product A*v where A=this and B is a dense matrix.
Eigen::VectorXd
ShiftRows::vecMultiplyRight(const Eigen::Map<Eigen::VectorXd> v, void (*checkInterrupt)(void)) {
    Eigen::VectorXd res = loader->vecMultiplyRight(v, checkInterrupt);
    res += fit.row_params.row(0).transpose().matrix() * v.sum();
    return res;
}
Eigen::VectorXd
ShiftRows::vecMultiplyLeft(const Eigen::Map<Eigen::VectorXd> v, void (*checkInterrupt)(void)) {
    Eigen::VectorXd res = loader->vecMultiplyLeft(v, checkInterrupt);
    res.rowwise() += fit.row_params.row(0).matrix() * v;
    return res;
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
Eigen::MatrixXd
ShiftCols::denseMultiplyRight(const Eigen::Map<Eigen::MatrixXd> B, void (*checkInterrupt)(void)) {
    Eigen::MatrixXd res = loader->denseMultiplyRight(B, checkInterrupt);
    res.rowwise() += fit.col_params.row(0).matrix() * B;
    return res;
}
Eigen::MatrixXd
ShiftCols::denseMultiplyLeft(const Eigen::Map<Eigen::MatrixXd> B, void (*checkInterrupt)(void)) {
    Eigen::MatrixXd res = loader->denseMultiplyLeft(B, checkInterrupt);
    res += B.rowwise().sum() * fit.col_params.row(0).matrix();
    return res;
}
// Calculate matrix-vector product A*v where A=this and B is a dense matrix.
Eigen::VectorXd
ShiftCols::vecMultiplyRight(const Eigen::Map<Eigen::VectorXd> v, void (*checkInterrupt)(void)) {
    Eigen::VectorXd res = loader->vecMultiplyRight(v, checkInterrupt);
    res.rowwise() += fit.col_params.row(0).matrix() * v;
    return res;
}
Eigen::VectorXd
ShiftCols::vecMultiplyLeft(const Eigen::Map<Eigen::VectorXd> v, void (*checkInterrupt)(void)) {
    Eigen::VectorXd res = loader->vecMultiplyLeft(v, checkInterrupt);
    res += fit.col_params.row(0).transpose().matrix() * v.sum();
    return res;
}

} // end namespace BPCells
