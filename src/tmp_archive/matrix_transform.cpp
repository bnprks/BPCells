#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <Rcpp.h>
#include <RcppEigen.h>

#include "matrixIterators/MatrixIterator.h"
#include "matrixTransforms/MatrixTransform.h"
#include "matrixTransforms/TFIDF.h"

using namespace BPCells;
using namespace Rcpp;

// [[Rcpp::export]]
List transform_get_fit_cpp(SEXP r_transform) {
    MatrixTransform *transform = &(*XPtr<MatrixTransform>(r_transform));
    TransformFit fit = transform->getFit();
    return List::create(
        Named("row_params") = fit.row_params,
        Named("col_params") = fit.col_params,
        Named("global_params") = fit.global_params
    );
}

// [[Rcpp::export]]
Eigen::VectorXd transform_vec_multiply_left_cpp(SEXP r_transform, Eigen::Map<Eigen::VectorXd> v) {
    MatrixTransform *transform = &(*XPtr<MatrixTransform>(r_transform));
    return transform->vecMultiplyLeft(v);
}

// [[Rcpp::export]]
Eigen::VectorXd transform_vec_multiply_right_cpp(SEXP r_transform, Eigen::Map<Eigen::VectorXd> v) {
    MatrixTransform *transform = &(*XPtr<MatrixTransform>(r_transform));
    return transform->vecMultiplyRight(v);
}

// [[Rcpp::export]]
Eigen::MatrixXd transform_dense_multiply_left_cpp(SEXP r_transform, Eigen::Map<Eigen::MatrixXd> v) {
    MatrixTransform *transform = &(*XPtr<MatrixTransform>(r_transform));
    return transform->denseMultiplyLeft(v);
}

// [[Rcpp::export]]
Eigen::MatrixXd
transform_dense_multiply_right_cpp(SEXP r_transform, Eigen::Map<Eigen::MatrixXd> v) {
    MatrixTransform *transform = &(*XPtr<MatrixTransform>(r_transform));
    return transform->denseMultiplyRight(v);
}

// [[Rcpp::export]]
SEXP transform_tfidf_cpp(SEXP r_matrix, double scale_to, bool transpose) {
    MatrixLoader<double> *matrix = &(*XPtr<MatrixLoader<double>>(r_matrix));

    return Rcpp::wrap(
        XPtr<MatrixTransform>(new TFIDF(*matrix, Eigen::ArrayXd::Constant(1, scale_to), transpose))
    );
}

// [[Rcpp::export]]
SEXP transform_project_tfidf_cpp(
    SEXP r_matrix, S4 r_transform, int recalculate_mode, bool transpose
) {
    MatrixLoader<double> *matrix = &(*XPtr<MatrixLoader<double>>(r_matrix));
    TransformFit fit{
        r_transform.slot("row_params"),
        r_transform.slot("col_params"),
        r_transform.slot("global_params")};
    auto mode = (MatrixTransform::RecalculateFit)recalculate_mode;
    return Rcpp::wrap(XPtr<MatrixTransform>(new TFIDF(*matrix, fit, mode, transpose)));
}