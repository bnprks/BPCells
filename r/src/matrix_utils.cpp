// Copyright 2021 BPCells contributors
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
#include "bpcells-cpp/matrixIterators/CSparseMatrix.h"
#include "bpcells-cpp/matrixIterators/ColwiseRank.h"
#include "bpcells-cpp/matrixIterators/ConcatenateMatrix.h"
#include "bpcells-cpp/matrixIterators/Mask.h"
#include "bpcells-cpp/matrixIterators/MatrixIndexSelect.h"
#include "bpcells-cpp/matrixIterators/MatrixIterator.h"
#include "bpcells-cpp/matrixIterators/MatrixMultiply.h"
#include "bpcells-cpp/matrixIterators/MatrixOps.h"
#include "bpcells-cpp/matrixIterators/MatrixStats.h"
#include "bpcells-cpp/matrixIterators/RenameDims.h"
#include "bpcells-cpp/matrixIterators/SVD.h"
#include "bpcells-cpp/matrixIterators/TSparseMatrixWriter.h"
#include "bpcells-cpp/matrixUtils/WilcoxonRankSum.h"
#include "bpcells-cpp/matrixUtils/Pseudobulk.h"
#include "bpcells-cpp/matrixUtils/Quantile.h"
#include "bpcells-cpp/matrixIterators/MatrixAccumulators.h"
#include "R_array_io.h"
#include "R_interrupts.h"
#include "R_xptr_wrapper.h"
#include <md5/md5.h>
#include <cstdio>
#include <string>

using namespace BPCells;
using namespace Rcpp;

// Add a protected handle to the underlying R data to prevent unwanted garbage collection
template <typename T> class RMatrixWrapper : public BPCells::MatrixLoaderWrapper<T> {
  private:
    Rcpp::RObject preserved_object;

  public:
    RMatrixWrapper(SEXP preserved_object, std::unique_ptr<MatrixLoader<T>> &&loader)
        : BPCells::MatrixLoaderWrapper<T>(std::move(loader))
        , preserved_object(preserved_object) {}
};

// [[Rcpp::export]]
SEXP iterate_csparse_matrix_cpp(
    SEXP matrix, const StringVector row_names, const StringVector col_names
) {
    auto eigen_mat = Rcpp::as<Eigen::Map<Eigen::SparseMatrix<double>>>(matrix);
    auto loader = std::make_unique<CSparseMatrix<double>>(
        eigen_mat,
        std::make_unique<RcppStringReader>(row_names),
        std::make_unique<RcppStringReader>(col_names)
    );

    return make_unique_xptr<RMatrixWrapper<double>>(Rcpp::RObject(matrix), std::move(loader));
}

template <typename From, typename To> SEXP convert_matrix_cpp(SEXP matrix) {
    return make_unique_xptr<MatrixConverterLoader<From, To>>(
        take_unique_xptr<MatrixLoader<From>>(matrix)
    );
}
// [[Rcpp::export]]
SEXP convert_matrix_uint32_t_double_cpp(SEXP matrix) {
    return convert_matrix_cpp<uint32_t, double>(matrix);
}
// [[Rcpp::export]]
SEXP convert_matrix_uint32_t_float_cpp(SEXP matrix) {
    return convert_matrix_cpp<uint32_t, float>(matrix);
}
// [[Rcpp::export]]
SEXP convert_matrix_double_uint32_t_cpp(SEXP matrix) {
    return convert_matrix_cpp<double, uint32_t>(matrix);
}
// [[Rcpp::export]]
SEXP convert_matrix_double_float_cpp(SEXP matrix) {
    return convert_matrix_cpp<double, float>(matrix);
}
// [[Rcpp::export]]
SEXP convert_matrix_float_uint32_t_cpp(SEXP matrix) {
    return convert_matrix_cpp<float, uint32_t>(matrix);
}
// [[Rcpp::export]]
SEXP convert_matrix_float_double_cpp(SEXP matrix) {
    return convert_matrix_cpp<float, double>(matrix);
}

// [[Rcpp::export]]
SEXP iterate_matrix_col_select_uint32_t_cpp(SEXP matrix, std::vector<uint32_t> col_selection) {
    return make_unique_xptr<MatrixColSelect<uint32_t>>(
        take_unique_xptr<MatrixLoader<uint32_t>>(matrix), col_selection
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_col_select_float_cpp(SEXP matrix, std::vector<uint32_t> col_selection) {
    return make_unique_xptr<MatrixColSelect<float>>(
        take_unique_xptr<MatrixLoader<float>>(matrix), col_selection
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_col_select_double_cpp(SEXP matrix, std::vector<uint32_t> col_selection) {
    return make_unique_xptr<MatrixColSelect<double>>(
        take_unique_xptr<MatrixLoader<double>>(matrix), col_selection
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_row_select_uint32_t_cpp(SEXP matrix, std::vector<uint32_t> row_selection) {
    return make_unique_xptr<MatrixRowSelect<uint32_t>>(
        take_unique_xptr<MatrixLoader<uint32_t>>(matrix), row_selection
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_row_select_float_cpp(SEXP matrix, std::vector<uint32_t> row_selection) {
    return make_unique_xptr<MatrixRowSelect<float>>(
        take_unique_xptr<MatrixLoader<float>>(matrix), row_selection
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_row_select_double_cpp(SEXP matrix, std::vector<uint32_t> row_selection) {
    return make_unique_xptr<MatrixRowSelect<double>>(
        take_unique_xptr<MatrixLoader<double>>(matrix), row_selection
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_rename_dims_uint32_t_cpp(
    SEXP matrix,
    std::vector<std::string> row_names,
    std::vector<std::string> col_names,
    bool clear_row_names,
    bool clear_col_names
) {
    return make_unique_xptr<RenameDims<uint32_t>>(
        take_unique_xptr<MatrixLoader<uint32_t>>(matrix),
        row_names,
        col_names,
        clear_row_names,
        clear_col_names
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_rename_dims_float_cpp(
    SEXP matrix,
    std::vector<std::string> row_names,
    std::vector<std::string> col_names,
    bool clear_row_names,
    bool clear_col_names
) {
    return make_unique_xptr<RenameDims<float>>(
        take_unique_xptr<MatrixLoader<float>>(matrix),
        row_names,
        col_names,
        clear_row_names,
        clear_col_names
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_rename_dims_double_cpp(
    SEXP matrix,
    std::vector<std::string> row_names,
    std::vector<std::string> col_names,
    bool clear_row_names,
    bool clear_col_names
) {
    return make_unique_xptr<RenameDims<double>>(
        take_unique_xptr<MatrixLoader<double>>(matrix),
        row_names,
        col_names,
        clear_row_names,
        clear_col_names
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_row_bind_uint32_t_cpp(SEXP matrix_list, int threads) {
    std::vector<std::unique_ptr<MatrixLoader<uint32_t>>> matrix_vec;
    List l = matrix_list;
    for (uint32_t i = 0; i < l.size(); i++) {
        SEXP elem = l[i];
        matrix_vec.push_back(take_unique_xptr<MatrixLoader<uint32_t>>(elem));
    }

    return make_unique_xptr<ConcatRows<uint32_t>>(std::move(matrix_vec), threads);
}

// [[Rcpp::export]]
SEXP iterate_matrix_row_bind_float_cpp(SEXP matrix_list, int threads) {
    std::vector<std::unique_ptr<MatrixLoader<float>>> matrix_vec;
    List l = matrix_list;
    for (uint32_t i = 0; i < l.size(); i++) {
        SEXP elem = l[i];
        matrix_vec.push_back(take_unique_xptr<MatrixLoader<float>>(elem));
    }

    return make_unique_xptr<ConcatRows<float>>(std::move(matrix_vec), threads);
}

// [[Rcpp::export]]
SEXP iterate_matrix_row_bind_double_cpp(SEXP matrix_list, int threads) {
    std::vector<std::unique_ptr<MatrixLoader<double>>> matrix_vec;
    List l = matrix_list;
    for (uint32_t i = 0; i < l.size(); i++) {
        SEXP elem = l[i];
        matrix_vec.push_back(take_unique_xptr<MatrixLoader<double>>(elem));
    }

    return make_unique_xptr<ConcatRows<double>>(std::move(matrix_vec), threads);
}

// [[Rcpp::export]]
SEXP iterate_matrix_col_bind_uint32_t_cpp(SEXP matrix_list, int threads) {
    std::vector<std::unique_ptr<MatrixLoader<uint32_t>>> matrix_vec;
    List l = matrix_list;
    for (uint32_t i = 0; i < l.size(); i++) {
        SEXP elem = l[i];
        matrix_vec.push_back(take_unique_xptr<MatrixLoader<uint32_t>>(elem));
    }

    return make_unique_xptr<ConcatCols<uint32_t>>(std::move(matrix_vec), threads);
}

// [[Rcpp::export]]
SEXP iterate_matrix_col_bind_float_cpp(SEXP matrix_list, int threads) {
    std::vector<std::unique_ptr<MatrixLoader<float>>> matrix_vec;
    List l = matrix_list;
    for (uint32_t i = 0; i < l.size(); i++) {
        SEXP elem = l[i];
        matrix_vec.push_back(take_unique_xptr<MatrixLoader<float>>(elem));
    }

    return make_unique_xptr<ConcatCols<float>>(std::move(matrix_vec), threads);
}

// [[Rcpp::export]]
SEXP iterate_matrix_col_bind_double_cpp(SEXP matrix_list, int threads) {
    std::vector<std::unique_ptr<MatrixLoader<double>>> matrix_vec;
    List l = matrix_list;
    for (uint32_t i = 0; i < l.size(); i++) {
        SEXP elem = l[i];
        matrix_vec.push_back(take_unique_xptr<MatrixLoader<double>>(elem));
    }

    return make_unique_xptr<ConcatCols<double>>(std::move(matrix_vec), threads);
}

// [[Rcpp::export]]
SEXP iterate_matrix_multiply_uint32_t_cpp(SEXP left, SEXP right) {
    return make_unique_xptr<SparseMultiply<uint32_t>>(
        take_unique_xptr<MatrixLoader<uint32_t>>(left),
        take_unique_xptr<MatrixLoader<uint32_t>>(right)
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_multiply_float_cpp(SEXP left, SEXP right) {
    return make_unique_xptr<SparseMultiply<float>>(
        take_unique_xptr<MatrixLoader<float>>(left), take_unique_xptr<MatrixLoader<float>>(right)
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_multiply_double_cpp(SEXP left, SEXP right) {
    return make_unique_xptr<SparseMultiply<double>>(
        take_unique_xptr<MatrixLoader<double>>(left), take_unique_xptr<MatrixLoader<double>>(right)
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_mask_uint32_t_cpp(SEXP mat, SEXP mask, bool invert) {
    if (invert)
        return make_unique_xptr<Mask<uint32_t, true>>(
            take_unique_xptr<MatrixLoader<uint32_t>>(mat),
            take_unique_xptr<MatrixLoader<uint32_t>>(mask)
        );
    else
        return make_unique_xptr<Mask<uint32_t, false>>(
            take_unique_xptr<MatrixLoader<uint32_t>>(mat),
            take_unique_xptr<MatrixLoader<uint32_t>>(mask)
        );
}

// [[Rcpp::export]]
SEXP iterate_matrix_mask_float_cpp(SEXP mat, SEXP mask, bool invert) {
    if (invert)
        return make_unique_xptr<Mask<float, true>>(
            take_unique_xptr<MatrixLoader<float>>(mat),
            take_unique_xptr<MatrixLoader<uint32_t>>(mask)
        );
    else
        return make_unique_xptr<Mask<float, false>>(
            take_unique_xptr<MatrixLoader<float>>(mat),
            take_unique_xptr<MatrixLoader<uint32_t>>(mask)
        );
}

// [[Rcpp::export]]
SEXP iterate_matrix_mask_double_cpp(SEXP mat, SEXP mask, bool invert) {
    if (invert)
        return make_unique_xptr<Mask<double, true>>(
            take_unique_xptr<MatrixLoader<double>>(mat),
            take_unique_xptr<MatrixLoader<uint32_t>>(mask)
        );
    else
        return make_unique_xptr<Mask<double, false>>(
            take_unique_xptr<MatrixLoader<double>>(mat),
            take_unique_xptr<MatrixLoader<uint32_t>>(mask)
        );
}

// [[Rcpp::export]]
SEXP iterate_matrix_rank_uint32_t_cpp(SEXP matrix) {
    return make_unique_xptr<ColwiseRank<uint32_t>>(take_unique_xptr<MatrixLoader<uint32_t>>(matrix)
    );
}

// [[Rcpp::export]]
SEXP iterate_matrix_rank_float_cpp(SEXP matrix) {
    return make_unique_xptr<ColwiseRank<float>>(take_unique_xptr<MatrixLoader<float>>(matrix));
}

// [[Rcpp::export]]
SEXP iterate_matrix_rank_double_cpp(SEXP matrix) {
    return make_unique_xptr<ColwiseRank<double>>(take_unique_xptr<MatrixLoader<double>>(matrix));
}

// [[Rcpp::export]]
Eigen::MatrixXd dense_multiply_right_cpp(SEXP matrix, Eigen::Map<Eigen::MatrixXd> B) {
    auto mat = take_unique_xptr<MatrixLoader<double>>(matrix);
    return run_with_R_interrupt_check(&MatrixLoader<double>::denseMultiplyRight, mat.get(), B);
}

// [[Rcpp::export]]
Eigen::MatrixXd dense_multiply_left_cpp(SEXP matrix, Eigen::Map<Eigen::MatrixXd> B) {
    auto mat = take_unique_xptr<MatrixLoader<double>>(matrix);
    return run_with_R_interrupt_check(&MatrixLoader<double>::denseMultiplyLeft, mat.get(), B);
}

// [[Rcpp::export]]
Eigen::VectorXd vec_multiply_right_cpp(SEXP matrix, Eigen::Map<Eigen::VectorXd> v) {
    auto mat = take_unique_xptr<MatrixLoader<double>>(matrix);
    return run_with_R_interrupt_check(&MatrixLoader<double>::vecMultiplyRight, mat.get(), v);
}

// [[Rcpp::export]]
Eigen::VectorXd vec_multiply_left_cpp(SEXP matrix, Eigen::Map<Eigen::VectorXd> v) {
    auto mat = take_unique_xptr<MatrixLoader<double>>(matrix);
    return run_with_R_interrupt_check(&MatrixLoader<double>::vecMultiplyLeft, mat.get(), v);
}

// [[Rcpp::export]]
Eigen::MatrixXd
dense_multiply_right_preserve_loader_cpp(SEXP matrix, Eigen::Map<Eigen::MatrixXd> B) {
    return run_with_R_interrupt_check(
        &MatrixLoader<double>::denseMultiplyRight, peek_unique_xptr<MatrixLoader<double>>(matrix), B
    );
}

// [[Rcpp::export]]
Eigen::MatrixXd
dense_multiply_left_preserve_loader_cpp(SEXP matrix, Eigen::Map<Eigen::MatrixXd> B) {
    return run_with_R_interrupt_check(
        &MatrixLoader<double>::denseMultiplyLeft, peek_unique_xptr<MatrixLoader<double>>(matrix), B
    );
}

// [[Rcpp::export]]
Eigen::VectorXd vec_multiply_right_preserve_loader_cpp(SEXP matrix, Eigen::Map<Eigen::VectorXd> v) {
    return run_with_R_interrupt_check(
        &MatrixLoader<double>::vecMultiplyRight, peek_unique_xptr<MatrixLoader<double>>(matrix), v
    );
}

// [[Rcpp::export]]
Eigen::VectorXd vec_multiply_left_preserve_loader_cpp(SEXP matrix, Eigen::Map<Eigen::VectorXd> v) {
    return run_with_R_interrupt_check(
        &MatrixLoader<double>::vecMultiplyLeft, peek_unique_xptr<MatrixLoader<double>>(matrix), v
    );
}

// [[Rcpp::export]]
std::vector<double> row_sums_double_cpp(SEXP matrix) {
    auto mat = take_unique_xptr<MatrixLoader<double>>(matrix);
    return run_with_R_interrupt_check(&MatrixLoader<double>::rowSums, mat.get());
}

// [[Rcpp::export]]
std::vector<double> col_sums_double_cpp(SEXP matrix) {
    auto mat = take_unique_xptr<MatrixLoader<double>>(matrix);
    return run_with_R_interrupt_check(&MatrixLoader<double>::colSums, mat.get());
}

// [[Rcpp::export]]
List matrix_stats_cpp(SEXP matrix, int row_stats, int col_stats) {
    auto mat = take_unique_xptr<MatrixLoader<double>>(matrix);
    StatsResult res = run_with_R_interrupt_check(
        &MatrixLoader<double>::computeMatrixStats, mat.get(), (Stats)row_stats, (Stats)col_stats
    );

    return List::create(Named("row_stats") = res.row_stats, Named("col_stats") = res.col_stats);
}

// [[Rcpp::export]]
Eigen::MatrixXd wilcoxon_rank_sum_pval_uint32_t_cpp(SEXP matrix, std::vector<uint32_t> groups) {
    return run_with_R_interrupt_check(
        &wilcoxon_rank_sum<uint32_t>,
        take_unique_xptr<MatrixLoader<uint32_t>>(matrix),
        std::cref(groups)
    );
}

// [[Rcpp::export]]
Eigen::MatrixXd wilcoxon_rank_sum_pval_float_cpp(SEXP matrix, std::vector<uint32_t> groups) {
    return run_with_R_interrupt_check(
        &wilcoxon_rank_sum<float>, take_unique_xptr<MatrixLoader<float>>(matrix), std::cref(groups)
    );
}

// [[Rcpp::export]]
Eigen::MatrixXd wilcoxon_rank_sum_pval_double_cpp(SEXP matrix, std::vector<uint32_t> groups) {
    return run_with_R_interrupt_check(
        &wilcoxon_rank_sum<double>,
        take_unique_xptr<MatrixLoader<double>>(matrix),
        std::cref(groups)
    );
}

// [[Rcpp::export]]
SEXP svds_cpp(SEXP matrix, int k, int n_cv, int maxit, double tol, bool transpose) {
    auto mat = take_unique_xptr<MatrixLoader<double>>(matrix);
    SVDResult res = run_with_R_interrupt_check(svd, mat.get(), k, n_cv, maxit, tol, transpose);
    if (!res.success) warning("SVD calculation did not converge");
    return List::create(
        Named("d") = res.d,
        Named("u") = res.u,
        Named("v") = res.v,
        Named("niter") = res.num_iterations,
        Named("nops") = res.num_operations,
        Named("nconv") = res.num_converged
    );
}

// Compute a histogram of the values in a matrix.
// [[Rcpp::export]]
NumericVector matrix_value_histogram_cpp(SEXP matrix, uint32_t max_value) {
    MatrixIterator<uint32_t> it(take_unique_xptr<MatrixLoader<uint32_t>>(matrix));
    std::vector<double> result(max_value + 1, 0.0);
    while (it.nextCol()) {
        Rcpp::checkUserInterrupt();
        while (it.nextValue()) {
            if (it.val() == 0) continue;
            result[std::min(it.val(), max_value + 1) - 1]++;
        }
    }
    return Rcpp::wrap(result);
}

// Check if we're currently running on an
bool is_little_endian() {
    int iv = 1;
    return *((char *)&iv) == 1;
}

// Compute the MD5 checksum of a double precision matrix.
// [[Rcpp::export]]
StringVector checksum_double_cpp(SEXP matrix) {
    MatrixIterator<double> it(take_unique_xptr<MatrixLoader<double>>(matrix));

    uint32_t nrow, ncol;
    md5_byte_t buffer[16];
    md5_byte_t digest[16];
    char hash[48];
    md5_state_t md5s;

    md5_init(&md5s);

    nrow = it.rows();
    ncol = it.cols();

    if (is_little_endian()) {
        // Little endian architecture.
        double *pval;
        uint32_t *prow;
        uint32_t *pcol;
        pval = (double *)buffer;
        prow = (uint32_t *)&(buffer[8]);
        pcol = (uint32_t *)&(buffer[12]);
        while (it.nextCol()) {
            Rcpp::checkUserInterrupt();
            while (it.nextValue()) {
                *pval = it.val();
                *prow = it.row();
                *pcol = it.col();
                md5_append(&md5s, (md5_byte_t *)buffer, (size_t)16);
            }
        }

        // Include numbers of rows and columns.
        md5_append(&md5s, (md5_byte_t *)&nrow, sizeof(nrow));
        md5_append(&md5s, (md5_byte_t *)&ncol, sizeof(ncol));
    } else {
        // Big endian architecture. Untested on 20240404.
        double dval;
        uint32_t ival;
        uint8_t *ptr;
        while (it.nextCol()) {
            Rcpp::checkUserInterrupt();
            while (it.nextValue()) {
                dval = it.val();
                ptr = (uint8_t *)&dval;
                buffer[7] = *(ptr);
                buffer[6] = *(++ptr);
                buffer[5] = *(++ptr);
                buffer[4] = *(++ptr);
                buffer[3] = *(++ptr);
                buffer[2] = *(++ptr);
                buffer[1] = *(++ptr);
                buffer[0] = *(++ptr);

                ival = it.row();
                ptr = (uint8_t *)&ival;
                buffer[11] = *(ptr);
                buffer[10] = *(++ptr);
                buffer[9] = *(++ptr);
                buffer[8] = *(++ptr);

                ival = it.col();
                ptr = (uint8_t *)&ival;
                buffer[15] = *(ptr);
                buffer[14] = *(++ptr);
                buffer[13] = *(++ptr);
                buffer[12] = *(++ptr);

                md5_append(&md5s, (md5_byte_t *)buffer, (size_t)16);
            }
        }
        // Include numbers of rows and columns.
        ival = (uint32_t)nrow;
        ptr = (uint8_t *)&ival;
        buffer[3] = *(ptr);
        buffer[2] = *(++ptr);
        buffer[1] = *(++ptr);
        buffer[0] = *(++ptr);
        md5_append(&md5s, (md5_byte_t *)buffer, (size_t)4);

        ival = (uint32_t)ncol;
        ptr = (uint8_t *)&ival;
        buffer[3] = *(ptr);
        buffer[2] = *(++ptr);
        buffer[1] = *(++ptr);
        buffer[0] = *(++ptr);
        md5_append(&md5s, (md5_byte_t *)buffer, (size_t)4);
    }

    // Include row and column names.
    for (uint32_t i = 0; i < nrow; ++i) {
      if (it.rowNames(i) != NULL) {
        md5_append(&md5s, (md5_byte_t *)it.rowNames(i), strlen(it.rowNames(i)));
      }
    }
    for (uint32_t i = 0; i < ncol; ++i) {
      if (it.colNames(i) != NULL) {
        md5_append(&md5s, (md5_byte_t *)it.colNames(i), strlen(it.colNames(i)));
      }
    }

    md5_finish(&md5s, digest);

    for (int i = 0; i < 16; ++i) {
        snprintf(&(hash[i * 2]), (size_t)3, "%02x", digest[i]);
    }
    StringVector svec(hash);
    return Rcpp::wrap(svec);
}

// [[Rcpp::export]]
List apply_matrix_double_cpp(SEXP mat_sexp, Function f, bool row_major) {
    std::unique_ptr<MatrixLoader<double>> mat = take_unique_xptr<MatrixLoader<double>>(mat_sexp);
    List ret(mat->cols());
    std::vector<char> called_col(mat->cols(), 0);

    std::vector<double> val;
    std::vector<int> row;
    while (mat->nextCol()) {
        int col = mat->currentCol();
        val.clear();
        row.clear();
        while (mat->load()) {
            row.insert(row.end(), mat->rowData(), mat->rowData() + mat->capacity());
            val.insert(val.end(), mat->valData(), mat->valData() + mat->capacity());
        }
        for (auto &r : row) r += 1;
        if (row_major) {
            ret[col] = f(val, col + 1, row);
        } else {
            ret[col] = f(val, row, col + 1);
        }
        called_col[col] = 1;
    }

    // Check any columns that got skipped and call them
    val.clear();
    row.clear();
    for (size_t i = 0; i < mat->cols(); i++) {
        if (!called_col[i]) {
            if (row_major) {
                ret[i] = f(val, i+1, row);
            } else {
                ret[i] = f(val, row, i+1);
            }
        }
    }

    return ret;
}

// Compute the max of each row
// [[Rcpp::export]]
NumericVector matrix_max_per_row_cpp(SEXP matrix) {
    MatrixIterator<double> it(take_unique_xptr<MatrixLoader<double>>(matrix));
    std::vector<double> result(it.rows());
    std::vector<uint32_t> row_count(it.rows(), 0);
    std::fill(result.begin(), result.end(), -INFINITY);
    // keep track of the number of times we've seen each row
    while (it.nextCol()) {
        Rcpp::checkUserInterrupt();
        while (it.nextValue()) {
            result[it.row()] = std::max(result[it.row()], it.val());
            row_count[it.row()]++;
        }
    }
    // If we've seen a row less than the num of cols, we know there's at least one 0
    for (size_t i = 0; i < it.rows(); i++) {
        if (row_count[i] < it.cols()) {
            result[i] = std::max(0.0, result[i]);
        }
    }
    return Rcpp::wrap(result);
}


// Compute the max of each col
// [[Rcpp::export]]
NumericVector matrix_max_per_col_cpp(SEXP matrix) {
    MatrixIterator<double> it(take_unique_xptr<MatrixLoader<double>>(matrix));
    std::vector<double> result(it.cols(), -INFINITY);
    // std::fill(result.begin(), result.end(), 0);
    // keep track of the number of times we've seen each col
    std::vector<uint32_t> col_count(it.cols(), 0);
    while (it.nextCol()) {
        Rcpp::checkUserInterrupt();
        while (it.nextValue()) {
            result[it.col()] = std::max(result[it.col()], it.val());
            col_count[it.col()]++;
        }
    }
    // If we've seen a col less than the num of rows, we know there's at least one 0
    for (size_t i = 0; i < it.cols(); i++) {
        if (col_count[i] < it.rows()) {
            result[i] = std::max(0.0, result[i]);
        }
    }
    return Rcpp::wrap(result);
}


// [[Rcpp::export]]
List pseudobulk_matrix_cpp(SEXP mat,
                           std::vector<uint32_t> cell_groups,
                           std::vector<std::string> method,
                           bool transpose) {
    PseudobulkStatsMethod methodFlags = static_cast<PseudobulkStatsMethod>(0);
    for (std::string &m : method) {
        if (m == "nonzeros") {
            methodFlags = methodFlags | PseudobulkStatsMethod::NonZeros;
        } else if (m == "sum") {
            methodFlags = methodFlags | PseudobulkStatsMethod::Sum;
        } else if (m == "mean") {
            methodFlags = methodFlags | PseudobulkStatsMethod::NonZeros | PseudobulkStatsMethod::Mean;
        } else if (m == "variance") {
            methodFlags = methodFlags | PseudobulkStatsMethod::NonZeros | PseudobulkStatsMethod::Mean | PseudobulkStatsMethod::Variance;
        }
    }
    
    PseudobulkStats res = run_with_R_interrupt_check(
        &pseudobulk_matrix<double>,
        take_unique_xptr<MatrixLoader<double>>(mat),
        std::cref(cell_groups),
        (PseudobulkStatsMethod)methodFlags,
        transpose
    );
    return List::create(
        Named("nonzeros") = res.non_zeros,
        Named("sum") = res.sum,
        Named("mean") = res.mean,
        Named("variance") = res.var
    );
}

// [[Rcpp::export]]
Eigen::ArrayXXd matrix_quantile_per_col_cpp(SEXP mat, std::vector<double> quantile, double alpha, double beta) {
    return run_with_R_interrupt_check(
        &matrix_quantile_per_col<double>,
        take_unique_xptr<MatrixLoader<double>>(mat),
        quantile,
        alpha,
        beta
    );
}


// [[Rcpp::export]]
bool matrix_identical_uint32_t_cpp(SEXP mat1, SEXP mat2) {
    MatrixIterator<uint32_t> i1(take_unique_xptr<MatrixLoader<uint32_t>>(mat1));
    MatrixIterator<uint32_t> i2(take_unique_xptr<MatrixLoader<uint32_t>>(mat2));
    i1.restart();
    i2.restart();

    while (true) {
        bool res1 = i1.nextCol();
        bool res2 = i2.nextCol();
        if (res1 != res2) {
            Rcerr << "Different number of columns." << std::endl;
            return false;
        }
        if (!res1) break;
        if (i1.currentCol() != i2.currentCol()) {
            Rcerr << "Different column loaded" << std::endl;
            return false;
        }
        uint32_t count = 0;
        while (true) {
            if (count++ % (1 << 14) == 0) Rcpp::checkUserInterrupt();
            bool res1 = i1.nextValue();
            bool res2 = i2.nextValue();
            if (res1 != res2) {
                Rcerr << "Different number of entries in column." << std::endl;
                return false;
            }
            if (!res1) break;
            if (i1.row() != i2.row() || i1.col() != i2.col() || i1.val() != i2.val()) {
                REprintf(
                    "Mismatched entries: (%d,%d=%d) vs. (%d,%d=%d)\n",
                    i1.row(),
                    i1.col(),
                    i1.val(),
                    i2.row(),
                    i2.col(),
                    i2.val()
                );
                return false;
            }
        }
    }
    // TODO: Check row/col names
    return true;
}
