#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <Rcpp.h>
#include <RcppEigen.h>

#include "matrixTransforms/Log1p.h"
#include "matrixTransforms/Pow.h"
#include "matrixTransforms/MatrixTransform.h"
#include "matrixTransforms/Min.h"
#include "matrixTransforms/Scale.h"
#include "matrixTransforms/Shift.h"

using namespace BPCells;
using namespace Rcpp;

// [[Rcpp::export]]
SEXP iterate_matrix_log1p_cpp(SEXP matrix) {
    XPtr<MatrixLoader<double>> input(matrix);
    return Rcpp::wrap(XPtr<MatrixLoader<double>>(new Log1p(*input)));
}

// [[Rcpp::export]]
SEXP iterate_matrix_log1psimd_cpp(SEXP matrix) {
    XPtr<MatrixLoader<double>> input(matrix);
    return Rcpp::wrap(XPtr<MatrixLoader<double>>(new Log1pSIMD(*input)));
}

// [[Rcpp::export]]
SEXP iterate_matrix_expm1_cpp(SEXP matrix) {
    XPtr<MatrixLoader<double>> input(matrix);
    return Rcpp::wrap(XPtr<MatrixLoader<double>>(new Expm1(*input)));
}

// [[Rcpp::export]]
SEXP iterate_matrix_expm1simd_cpp(SEXP matrix) {
    XPtr<MatrixLoader<double>> input(matrix);
    return Rcpp::wrap(XPtr<MatrixLoader<double>>(new Expm1SIMD(*input)));
}

// [[Rcpp::export]]
SEXP iterate_matrix_min_cpp(SEXP matrix, double min_val) {
    XPtr<MatrixLoader<double>> input(matrix);
    Eigen::ArrayXd global_params(1);
    global_params = min_val;
    return Rcpp::wrap(
        XPtr<MatrixLoader<double>>(new Min(*input, TransformFit{{}, {}, global_params}))
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_pow_cpp(SEXP matrix, double exponent) {
    XPtr<MatrixLoader<double>> input(matrix);
    Eigen::ArrayXd global_params(1);
    global_params = exponent;
    return Rcpp::wrap(XPtr<MatrixLoader<double>>(new Pow(*input, TransformFit{{}, {}, global_params})));
}

// [[Rcpp::export]]
SEXP iterate_matrix_powsimd_cpp(SEXP matrix, double exponent) {
    XPtr<MatrixLoader<double>> input(matrix);
    Eigen::ArrayXd global_params(1);
    global_params = exponent;
    return Rcpp::wrap(XPtr<MatrixLoader<double>>(new PowSIMD(*input, TransformFit{{}, {}, global_params})));
}

// [[Rcpp::export]]
SEXP iterate_matrix_square_cpp(SEXP matrix) {
    XPtr<MatrixLoader<double>> input(matrix);
    return Rcpp::wrap(XPtr<MatrixLoader<double>>(new Square(*input)));
}

// [[Rcpp::export]]
SEXP iterate_matrix_scale_cpp(
    SEXP matrix, Eigen::Map<Eigen::ArrayXXd> row_scale, Eigen::Map<Eigen::ArrayXXd> col_scale
) {
    XPtr<MatrixLoader<double>> input(matrix);
    return Rcpp::wrap(
        XPtr<MatrixLoader<double>>(new Scale(*input, TransformFit{row_scale, col_scale}))
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_row_shift_cpp(SEXP matrix, Eigen::Map<Eigen::ArrayXXd> row_shift) {
    XPtr<MatrixLoader<double>> input(matrix);
    return Rcpp::wrap(XPtr<MatrixLoader<double>>(new ShiftRows(*input, TransformFit{row_shift})));
}

// [[Rcpp::export]]
SEXP iterate_matrix_col_shift_cpp(SEXP matrix, Eigen::Map<Eigen::ArrayXXd> col_shift) {
    XPtr<MatrixLoader<double>> input(matrix);
    return Rcpp::wrap(XPtr<MatrixLoader<double>>(
        new ShiftCols(*input, TransformFit{Eigen::ArrayXXd(), col_shift})
    ));
}
