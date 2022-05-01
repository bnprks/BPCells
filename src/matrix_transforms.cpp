#include <Rcpp.h>
#include <RcppEigen.h>

#include "matrixTransforms/MatrixTransform.h"
#include "matrixTransforms/Log1p.h"
#include "matrixTransforms/Scale.h"
#include "matrixTransforms/Shift.h"

using namespace BPCells;
using namespace Rcpp;

// [[Rcpp::export]]
SEXP iterate_matrix_log1p_cpp(SEXP matrix) {
    XPtr<MatrixLoader<double>>input(matrix);
    return Rcpp::wrap(
        XPtr<MatrixLoader<double>>(new Log1p(*input))
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_scale_cpp(SEXP matrix, Eigen::Map<Eigen::ArrayXXd> row_scale, Eigen::Map<Eigen::ArrayXXd> col_scale) {
    XPtr<MatrixLoader<double>>input(matrix);
    return Rcpp::wrap(
        XPtr<MatrixLoader<double>>(new Scale(*input, TransformFit{row_scale, col_scale}))
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_row_shift_cpp(SEXP matrix, Eigen::Map<Eigen::ArrayXXd> row_shift) {
    XPtr<MatrixLoader<double>>input(matrix);
    return Rcpp::wrap(
        XPtr<MatrixLoader<double>>(new ShiftRows(*input, TransformFit{row_shift}))
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_col_shift_cpp(SEXP matrix, Eigen::Map<Eigen::ArrayXXd> col_shift) {
    XPtr<MatrixLoader<double>>input(matrix);
    return Rcpp::wrap(
        XPtr<MatrixLoader<double>>(new ShiftCols(*input, TransformFit{{NULL, 0, 0}, col_shift}))
    );
}