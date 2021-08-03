#include <Rcpp.h>
#include <RcppEigen.h>

#include "matrixIterators/MatrixIterator.h"
#include "matrixIterators/TSparseMatrixWriter.h"
#include "matrixIterators/CSparseMatrix.h"
#include "matrixIterators/MatrixOps.h"

#include "matrixTransforms/MatrixStats.h"

using namespace BPCells;
using namespace Rcpp;

// [[Rcpp::export]]
SEXP iterate_csparse_matrix_cpp(SEXP matrix) {
    auto eigen_mat = Rcpp::as<Eigen::Map<Eigen::SparseMatrix<double>>>(matrix);
    return Rcpp::wrap(
        XPtr<MatrixLoader<double>>(new CSparseMatrix(eigen_mat))
    );
}

// [[Rcpp::export]]
List build_tsparse_matrix_uint32_t_cpp(SEXP matrix) {
    
    MatrixLoader<uint32_t> *loader = &(*XPtr<MatrixLoader<uint32_t> >(matrix));
    MatrixIterator<uint32_t> iter(*loader);

    TSparseMatrixWriter<uint32_t> out;
    out.write(iter, &Rcpp::checkUserInterrupt);
    
    return List::create(
        Named("row") = Rcpp::IntegerVector((int32_t *) &out.rows[0], (int32_t *) &out.rows[out.rows.size()]),
        Named("col") = Rcpp::IntegerVector((int32_t *) &out.cols[0], (int32_t *) &out.cols[out.cols.size()]),
        Named("val") = Rcpp::IntegerVector((int32_t *) &out.vals[0], (int32_t *) &out.vals[out.vals.size()])
    );
}

// [[Rcpp::export]]
SEXP convert_matrix_uint32_t_double_cpp(SEXP matrix) {
    MatrixLoader<uint32_t> *input = &(*XPtr<MatrixLoader<uint32_t> >(matrix));

    return Rcpp::wrap(
        XPtr<MatrixLoader<double>>(new MatrixConverterLoader<uint32_t, double>(*input))
    );
}

// [[Rcpp::export]]
SEXP convert_matrix_double_uint32_t_cpp(SEXP matrix) {
    MatrixLoader<double> *input = &(*XPtr<MatrixLoader<double> >(matrix));

    return Rcpp::wrap(
        XPtr<MatrixLoader<uint32_t>>(new MatrixConverterLoader<double, uint32_t>(*input))
    );
}

// [[Rcpp::export]]
SEXP build_csparse_matrix_double_cpp(SEXP matrix) {
    MatrixLoader<double> *input = &(*XPtr<MatrixLoader<double> >(matrix));
    MatrixIterator<double> iter(*input);
    
    CSparseMatrixWriter writer;
    writer.write(iter, &Rcpp::checkUserInterrupt);
    return Rcpp::wrap(writer.getMat());
}

// [[Rcpp::export]]
NumericVector scan_matrix_uint32_t_cpp(SEXP matrix) {
    MatrixLoader<uint32_t> *loader = &(*XPtr<MatrixLoader<uint32_t> >(matrix));
    MatrixIterator<uint32_t> iter(*loader);
    uint64_t entries = 0;
    uint64_t val_sum = 0;
    uint64_t row_sum = 0;
    uint64_t col_sum = 0;
    while (iter.nextCol()) {
        while (iter.nextValue()) {
            entries ++;
            val_sum += iter.val();
            row_sum += iter.row();
            col_sum += iter.col();
            if (entries % 10000 == 0) Rcpp::checkUserInterrupt();
        }
    }
    return NumericVector({(double) entries,(double) val_sum,(double) row_sum,(double) col_sum});
}

// [[Rcpp::export]]
Eigen::MatrixXd dense_multiply_right_cpp(SEXP matrix, Eigen::Map<Eigen::MatrixXd> B) {
    MatrixLoader<double> *loader = &(*XPtr<MatrixLoader<double> >(matrix));
    MatrixIterator<double> A(*loader);
    return denseMultiplyRight(A, B);
}

// [[Rcpp::export]]
Eigen::MatrixXd dense_multiply_left_cpp(SEXP matrix, Eigen::Map<Eigen::MatrixXd> B) {
    MatrixLoader<double> *loader = &(*XPtr<MatrixLoader<double> >(matrix));
    MatrixIterator<double> A(*loader);
    return denseMultiplyLeft(A, B);
}

// [[Rcpp::export]]
Eigen::VectorXd vec_multiply_right_cpp(SEXP matrix, Eigen::Map<Eigen::VectorXd> v) {
    MatrixLoader<double> *loader = &(*XPtr<MatrixLoader<double> >(matrix));
    MatrixIterator<double> A(*loader);
    return vecMultiplyRight(A, v);
}

// [[Rcpp::export]]
Eigen::VectorXd vec_multiply_left_cpp(SEXP matrix, Eigen::Map<Eigen::VectorXd> v) {
    MatrixLoader<double> *loader = &(*XPtr<MatrixLoader<double> >(matrix));
    MatrixIterator<double> A(*loader);
    return vecMultiplyLeft(A, v);
}

// [[Rcpp::export]]
std::vector<double> row_sums_cpp(SEXP matrix) {
    MatrixLoader<double> *loader = &(*XPtr<MatrixLoader<double> >(matrix));
    MatrixIterator<double> A(*loader);
    return rowSums<double>(A);
}

// [[Rcpp::export]]
std::vector<double> col_sums_cpp(SEXP matrix) {
    MatrixLoader<double> *loader = &(*XPtr<MatrixLoader<double> >(matrix));
    MatrixIterator<double> A(*loader);
    return colSums<double>(A);
}

// [[Rcpp::export]]
List matrix_stats_cpp(SEXP matrix, int row_stats, int col_stats) {
    MatrixLoader<double> *loader = &(*XPtr<MatrixLoader<double> >(matrix));
    
    StatsResult res = computeMatrixStats(*loader, (Stats) row_stats, (Stats) col_stats);

    return List::create(
        Named("row_stats") = res.row_stats,
        Named("col_stats") = res.col_stats
    );
}