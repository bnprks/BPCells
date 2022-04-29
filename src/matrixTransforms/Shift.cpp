#include "Shift.h"

namespace BPCells {

// Shift rows a matrix
// out[i,j] = in[i,j] + row_params[0,i]
ShiftRows::ShiftRows(MatrixLoader<double> &mat, TransformFit fit)
    : MatrixTransformDense(mat, fit) {}

bool ShiftRows::loadZeroSubtracted() {
    if (!loader.load()) return false;
    // The transform(data) - transform(0) = data by definition
    return true;
}

void ShiftRows::loadZero(double *values, uint32_t count, uint32_t start_row, uint32_t col) {
    for (uint32_t i = 0; i < count; i++) {
        values[i] = fit.row_params(0, start_row + i);
    }
}



// Shift rows a matrix
// out[i,j] = in[i,j] + col_params[0,j]
ShiftCols::ShiftCols(MatrixLoader<double> &mat, TransformFit fit)
    : MatrixTransformDense(mat, fit) {}

bool ShiftCols::loadZeroSubtracted() {
    if (!loader.load()) return false;
    // The transform(data) - transform(0) = data by definition
    return true;
}

void ShiftCols::loadZero(double *values, uint32_t count, uint32_t start_row, uint32_t col) {
    double val = fit.col_params(0, col);
    for (uint32_t i = 0; i < count; i++) {
        values[i] = val;
    }
}


} // end namespace BPCells