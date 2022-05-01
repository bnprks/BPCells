#include <Rcpp.h>
#include <RcppEigen.h>

#include "matrixIterators/MatrixIterator.h"
#include "matrixIterators/TSparseMatrixWriter.h"
#include "matrixIterators/CSparseMatrix.h"
#include "matrixIterators/MatrixIndexSelect.h"
#include "matrixIterators/MatrixOps.h"
#include "matrixIterators/MatrixStats.h"


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
SEXP convert_matrix_uint32_t_double_cpp(SEXP matrix) {
    XPtr<MatrixLoader<uint32_t>>input(matrix);

    return Rcpp::wrap(
        XPtr<MatrixLoader<double>>(new MatrixConverterLoader<uint32_t, double>(*input))
    );
}

// [[Rcpp::export]]
SEXP convert_matrix_double_uint32_t_cpp(SEXP matrix) {
    XPtr<MatrixLoader<double>>input(matrix);

    return Rcpp::wrap(
        XPtr<MatrixLoader<uint32_t>>(new MatrixConverterLoader<double, uint32_t>(*input))
    );
}

// [[Rcpp::export]]
SEXP build_csparse_matrix_double_cpp(SEXP matrix) {
    XPtr<MatrixLoader<double>>input(matrix);
    MatrixIterator<double> iter(*input);
    
    CSparseMatrixWriter writer;
    writer.write(iter, &Rcpp::checkUserInterrupt);
    return Rcpp::wrap(writer.getMat());
}

// [[Rcpp::export]]
SEXP iterate_matrix_col_select_double_cpp(SEXP matrix, std::vector<uint32_t> col_selection) {
    XPtr<MatrixLoader<double>> loader(matrix);
    return Rcpp::wrap(
        XPtr<MatrixLoader<double>>(new MatrixColSelect(*loader, col_selection))
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_col_select_uint32_t_cpp(SEXP matrix, std::vector<uint32_t> col_selection) {
    XPtr<MatrixLoader<uint32_t>> loader(matrix);
    return Rcpp::wrap(
        XPtr<MatrixLoader<uint32_t>>(new MatrixColSelect(*loader, col_selection))
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_row_select_double_cpp(SEXP matrix, std::vector<uint32_t> row_selection) {
    XPtr<MatrixLoader<double>> loader(matrix);
    return Rcpp::wrap(
        XPtr<MatrixLoader<double>>(new MatrixRowSelect(*loader, row_selection))
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_row_select_uint32_t_cpp(SEXP matrix, std::vector<uint32_t> row_selection) {
    XPtr<MatrixLoader<uint32_t>> loader(matrix);
    return Rcpp::wrap(
        XPtr<MatrixLoader<uint32_t>>(new MatrixRowSelect(*loader, row_selection))
    );
}

// [[Rcpp::export]]
NumericVector scan_matrix_uint32_t_cpp(SEXP matrix) {
    XPtr<MatrixLoader<uint32_t>>loader(matrix);
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
    XPtr<MatrixLoader<double>>loader(matrix);
    return loader->denseMultiplyRight(B, &Rcpp::checkUserInterrupt);
}

// [[Rcpp::export]]
Eigen::MatrixXd dense_multiply_left_cpp(SEXP matrix, Eigen::Map<Eigen::MatrixXd> B) {
    XPtr<MatrixLoader<double>>loader(matrix);
    return loader->denseMultiplyLeft(B, &Rcpp::checkUserInterrupt);
}

// [[Rcpp::export]]
Eigen::VectorXd vec_multiply_right_cpp(SEXP matrix, Eigen::Map<Eigen::VectorXd> v) {
    XPtr<MatrixLoader<double>>loader(matrix);
    return loader->vecMultiplyRight(v, &Rcpp::checkUserInterrupt);
}

// [[Rcpp::export]]
Eigen::VectorXd vec_multiply_left_cpp(SEXP matrix, Eigen::Map<Eigen::VectorXd> v) {
    XPtr<MatrixLoader<double>>loader(matrix);
    return loader->vecMultiplyLeft(v, &Rcpp::checkUserInterrupt);
}

// [[Rcpp::export]]
std::vector<double> row_sums_double_cpp(SEXP matrix) {
    XPtr<MatrixLoader<double>>loader(matrix);
    return loader->rowSums(&Rcpp::checkUserInterrupt);
}

// [[Rcpp::export]]
std::vector<double> col_sums_double_cpp(SEXP matrix) {
    XPtr<MatrixLoader<double>>loader(matrix);
    return loader->colSums(&Rcpp::checkUserInterrupt);
}

// [[Rcpp::export]]
List matrix_stats_cpp(SEXP matrix, int row_stats, int col_stats) {
    XPtr<MatrixLoader<double>>loader(matrix);
    StatsResult res = loader->computeMatrixStats((Stats) row_stats, (Stats) col_stats, &Rcpp::checkUserInterrupt);

    return List::create(
        Named("row_stats") = res.row_stats,
        Named("col_stats") = res.col_stats
    );
}

// [[Rcpp::export]]
bool matrix_identical_uint32_t_cpp(SEXP mat1, SEXP mat2) {
    XPtr<MatrixLoader<uint32_t>> l1(mat1);
    XPtr<MatrixLoader<uint32_t>> l2(mat2);
    l1->restart(); l2->restart();
    MatrixIterator<uint32_t> i1(*l1);
    MatrixIterator<uint32_t> i2(*l2);

    while(true) {
        bool res1 = i1.nextCol();
        bool res2 = i2.nextCol();
        if(res1 != res2) {
            Rcerr << "Different number of columns." << std::endl;
            return false;
        }
        if (!res1) break;
        if(i1.currentCol() != i2.currentCol()) {
            Rcerr << "Different column loaded" << std::endl;
            return false;
        }
        while(true) {
            bool res1 = i1.nextValue();
            bool res2 = i2.nextValue();
            if (res1 != res2) {
                Rcerr << "Different number of entries in column." << std::endl;
                return false;
            }
            if (!res1) break;
            if (i1.row() != i2.row() || i1.col() != i2.col() || i1.val() != i2.val()) {
                REprintf("Mismatched entries: (%d,%d=%d) vs. (%d,%d=%d)\n",
                    i1.row(), i1.col(), i1.val(),
                    i2.row(), i2.col(), i2.val()
                );
                return false;
            }
        }
    }
    // TODO: Check row/col names
    return true;
}