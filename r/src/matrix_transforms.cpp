// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <Rcpp.h>
#include <RcppEigen.h>

#include "bpcells-cpp/matrixTransforms/Binarize.h"
#include "bpcells-cpp/matrixTransforms/Log1p.h"
#include "bpcells-cpp/matrixTransforms/MatrixTransform.h"
#include "bpcells-cpp/matrixTransforms/Min.h"
#include "bpcells-cpp/matrixTransforms/Pow.h"
#include "bpcells-cpp/matrixTransforms/Round.h"
#include "bpcells-cpp/matrixTransforms/SCTransform.h"
#include "bpcells-cpp/matrixTransforms/Scale.h"
#include "bpcells-cpp/matrixTransforms/Shift.h"

#include "R_xptr_wrapper.h"

using namespace BPCells;
using namespace Rcpp;

// [[Rcpp::export]]
SEXP iterate_matrix_log1p_cpp(SEXP matrix) {
    return make_unique_xptr<Log1p>(take_unique_xptr<MatrixLoader<double>>(matrix));
}

// [[Rcpp::export]]
SEXP iterate_matrix_log1psimd_cpp(SEXP matrix) {
    return make_unique_xptr<Log1pSIMD>(take_unique_xptr<MatrixLoader<double>>(matrix));
}

// [[Rcpp::export]]
SEXP iterate_matrix_expm1_cpp(SEXP matrix) {
    return make_unique_xptr<Expm1>(take_unique_xptr<MatrixLoader<double>>(matrix));
}

// [[Rcpp::export]]
SEXP iterate_matrix_expm1simd_cpp(SEXP matrix) {
    return make_unique_xptr<Expm1SIMD>(take_unique_xptr<MatrixLoader<double>>(matrix));
}

// [[Rcpp::export]]
SEXP iterate_matrix_pow_cpp(SEXP matrix, double exponent) {
    Eigen::ArrayXd global_params(1);
    global_params = exponent;
    return make_unique_xptr<Pow>(
        take_unique_xptr<MatrixLoader<double>>(matrix), TransformFit{{}, {}, global_params}
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_square_cpp(SEXP matrix) {
    return make_unique_xptr<Square>(take_unique_xptr<MatrixLoader<double>>(matrix));
}

// [[Rcpp::export]]
SEXP iterate_matrix_squaresimd_cpp(SEXP matrix) {
    return make_unique_xptr<SquareSIMD>(take_unique_xptr<MatrixLoader<double>>(matrix));
}

// [[Rcpp::export]]
SEXP iterate_matrix_min_cpp(SEXP matrix, double min_val) {
    Eigen::ArrayXd global_params(1);
    global_params = min_val;
    return make_unique_xptr<Min>(
        take_unique_xptr<MatrixLoader<double>>(matrix), TransformFit{{}, {}, global_params}
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_min_by_row_cpp(SEXP matrix, Eigen::ArrayXXd row_min) {
    return make_unique_xptr<MinByRow>(
        take_unique_xptr<MatrixLoader<double>>(matrix), TransformFit{row_min}
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_min_by_col_cpp(SEXP matrix, Eigen::ArrayXXd col_min) {
    return make_unique_xptr<MinByCol>(
        take_unique_xptr<MatrixLoader<double>>(matrix), TransformFit{{}, col_min}
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_binarize_cpp(SEXP matrix, double threshold, uint32_t strict_inequality) {
    Eigen::ArrayXd global_params(2);
    global_params[0] = threshold;
    global_params[1] = strict_inequality;
    return make_unique_xptr<Binarize>(
        take_unique_xptr<MatrixLoader<double>>(matrix), TransformFit{{}, {}, global_params}
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_round_cpp(SEXP matrix, uint32_t digits) {
    Eigen::ArrayXd global_params(1);
    global_params[0] = digits;
    return make_unique_xptr<Round>(
        take_unique_xptr<MatrixLoader<double>>(matrix), TransformFit{{}, {}, global_params}
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_sctransform_pearson_cpp(
    SEXP matrix,
    Eigen::Map<Eigen::ArrayXXd> gene_params,
    Eigen::Map<Eigen::ArrayXXd> cell_params,
    Eigen::Map<Eigen::ArrayXd> global_params
) {
    return make_unique_xptr<SCTransformPearson>(
        take_unique_xptr<MatrixLoader<double>>(matrix),
        TransformFit{gene_params, cell_params, global_params}
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_sctransform_pearson_transpose_cpp(
    SEXP matrix,
    Eigen::Map<Eigen::ArrayXXd> cell_params,
    Eigen::Map<Eigen::ArrayXXd> gene_params,
    Eigen::Map<Eigen::ArrayXd> global_params
) {
    return make_unique_xptr<SCTransformPearsonTranspose>(
        take_unique_xptr<MatrixLoader<double>>(matrix),
        TransformFit{cell_params, gene_params, global_params}
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_sctransform_pearson_simd_cpp(
    SEXP matrix,
    Eigen::Map<Eigen::ArrayXXd> gene_params,
    Eigen::Map<Eigen::ArrayXXd> cell_params,
    Eigen::Map<Eigen::ArrayXd> global_params
) {
    return make_unique_xptr<SCTransformPearsonSIMD>(
        take_unique_xptr<MatrixLoader<double>>(matrix),
        TransformFit{gene_params, cell_params, global_params}
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_sctransform_pearson_transpose_simd_cpp(
    SEXP matrix,
    Eigen::Map<Eigen::ArrayXXd> cell_params,
    Eigen::Map<Eigen::ArrayXXd> gene_params,
    Eigen::Map<Eigen::ArrayXd> global_params
) {
    return make_unique_xptr<SCTransformPearsonTransposeSIMD>(
        take_unique_xptr<MatrixLoader<double>>(matrix),
        TransformFit{cell_params, gene_params, global_params}
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_scale_cpp(
    SEXP matrix, Eigen::Map<Eigen::ArrayXXd> row_scale, Eigen::Map<Eigen::ArrayXXd> col_scale
) {
    return make_unique_xptr<Scale>(
        take_unique_xptr<MatrixLoader<double>>(matrix), TransformFit{row_scale, col_scale}
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_row_shift_cpp(SEXP matrix, Eigen::Map<Eigen::ArrayXXd> row_shift) {
    return make_unique_xptr<ShiftRows>(
        take_unique_xptr<MatrixLoader<double>>(matrix), TransformFit{row_shift}
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_col_shift_cpp(SEXP matrix, Eigen::Map<Eigen::ArrayXXd> col_shift) {
    return make_unique_xptr<ShiftCols>(
        take_unique_xptr<MatrixLoader<double>>(matrix), TransformFit{Eigen::ArrayXXd(), col_shift}
    );
}
