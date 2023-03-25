#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <Rcpp.h>
#include <RcppEigen.h>

#include "matrixIterators/CSparseMatrix.h"
#include "matrixIterators/ColwiseRank.h"
#include "matrixIterators/ConcatenateMatrix.h"
#include "matrixIterators/Mask.h"
#include "matrixIterators/MatrixIndexSelect.h"
#include "matrixIterators/MatrixIterator.h"
#include "matrixIterators/MatrixMultiply.h"
#include "matrixIterators/MatrixOps.h"
#include "matrixIterators/MatrixStats.h"
#include "matrixIterators/TSparseMatrixWriter.h"
#include "matrixIterators/WilcoxonRankSum.h"

#include "R_array_io.h"
#include "R_xptr_wrapper.h"
#include "R_interrupts.h"

using namespace BPCells;
using namespace Rcpp;

// [[Rcpp::export]]
SEXP iterate_csparse_matrix_cpp(
    SEXP matrix, const StringVector row_names, const StringVector col_names
) {
    auto eigen_mat = Rcpp::as<Eigen::Map<Eigen::SparseMatrix<double>>>(matrix);
    return make_unique_xptr<CSparseMatrix>(
        eigen_mat,
        std::make_unique<RcppStringReader>(row_names),
        std::make_unique<RcppStringReader>(col_names)
    );
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
SEXP build_csparse_matrix_double_cpp(SEXP matrix) {
    CSparseMatrixWriter writer;
    run_with_R_interrupt_check(&CSparseMatrixWriter::write, &writer, std::ref(*peek_unique_xptr<MatrixLoader<double>>(matrix)));
    return Rcpp::wrap(writer.getMat());
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
SEXP iterate_matrix_row_bind_uint32_t_cpp(SEXP matrix_list) {
    std::vector<std::unique_ptr<MatrixLoader<uint32_t>>> matrix_vec;
    List l = matrix_list;
    for (uint32_t i = 0; i < l.size(); i++) {
        SEXP elem = l[i];
        matrix_vec.push_back(take_unique_xptr<MatrixLoader<uint32_t>>(elem));
    }

    return make_unique_xptr<ConcatRows<uint32_t>>(std::move(matrix_vec));
}

// [[Rcpp::export]]
SEXP iterate_matrix_row_bind_float_cpp(SEXP matrix_list) {
    std::vector<std::unique_ptr<MatrixLoader<float>>> matrix_vec;
    List l = matrix_list;
    for (uint32_t i = 0; i < l.size(); i++) {
        SEXP elem = l[i];
        matrix_vec.push_back(take_unique_xptr<MatrixLoader<float>>(elem));
    }

    return make_unique_xptr<ConcatRows<float>>(std::move(matrix_vec));
}

// [[Rcpp::export]]
SEXP iterate_matrix_row_bind_double_cpp(SEXP matrix_list) {
    std::vector<std::unique_ptr<MatrixLoader<double>>> matrix_vec;
    List l = matrix_list;
    for (uint32_t i = 0; i < l.size(); i++) {
        SEXP elem = l[i];
        matrix_vec.push_back(take_unique_xptr<MatrixLoader<double>>(elem));
    }

    return make_unique_xptr<ConcatRows<double>>(std::move(matrix_vec));
}

// [[Rcpp::export]]
SEXP iterate_matrix_col_bind_uint32_t_cpp(SEXP matrix_list) {
    std::vector<std::unique_ptr<MatrixLoader<uint32_t>>> matrix_vec;
    List l = matrix_list;
    for (uint32_t i = 0; i < l.size(); i++) {
        SEXP elem = l[i];
        matrix_vec.push_back(take_unique_xptr<MatrixLoader<uint32_t>>(elem));
    }

    return make_unique_xptr<ConcatCols<uint32_t>>(std::move(matrix_vec));
}

// [[Rcpp::export]]
SEXP iterate_matrix_col_bind_float_cpp(SEXP matrix_list) {
    std::vector<std::unique_ptr<MatrixLoader<float>>> matrix_vec;
    List l = matrix_list;
    for (uint32_t i = 0; i < l.size(); i++) {
        SEXP elem = l[i];
        matrix_vec.push_back(take_unique_xptr<MatrixLoader<float>>(elem));
    }

    return make_unique_xptr<ConcatCols<float>>(std::move(matrix_vec));
}

// [[Rcpp::export]]
SEXP iterate_matrix_col_bind_double_cpp(SEXP matrix_list) {
    std::vector<std::unique_ptr<MatrixLoader<double>>> matrix_vec;
    List l = matrix_list;
    for (uint32_t i = 0; i < l.size(); i++) {
        SEXP elem = l[i];
        matrix_vec.push_back(take_unique_xptr<MatrixLoader<double>>(elem));
    }

    return make_unique_xptr<ConcatCols<double>>(std::move(matrix_vec));
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
    return run_with_R_interrupt_check(
        &MatrixLoader<double>::denseMultiplyRight, peek_unique_xptr<MatrixLoader<double>>(matrix), B
    );
}

// [[Rcpp::export]]
Eigen::MatrixXd dense_multiply_left_cpp(SEXP matrix, Eigen::Map<Eigen::MatrixXd> B) {
    return run_with_R_interrupt_check(
        &MatrixLoader<double>::denseMultiplyLeft, peek_unique_xptr<MatrixLoader<double>>(matrix), B
    );
}

// [[Rcpp::export]]
Eigen::VectorXd vec_multiply_right_cpp(SEXP matrix, Eigen::Map<Eigen::VectorXd> v) {
    return run_with_R_interrupt_check(
        &MatrixLoader<double>::vecMultiplyRight, peek_unique_xptr<MatrixLoader<double>>(matrix), v
    );
}

// [[Rcpp::export]]
Eigen::VectorXd vec_multiply_left_cpp(SEXP matrix, Eigen::Map<Eigen::VectorXd> v) {
    return run_with_R_interrupt_check(
        &MatrixLoader<double>::vecMultiplyLeft, peek_unique_xptr<MatrixLoader<double>>(matrix), v
    );
}

// [[Rcpp::export]]
std::vector<double> row_sums_double_cpp(SEXP matrix) {
    return run_with_R_interrupt_check(&MatrixLoader<double>::rowSums, peek_unique_xptr<MatrixLoader<double>>(matrix));
}

// [[Rcpp::export]]
std::vector<double> col_sums_double_cpp(SEXP matrix) {
    return run_with_R_interrupt_check(&MatrixLoader<double>::colSums, peek_unique_xptr<MatrixLoader<double>>(matrix));
}

// [[Rcpp::export]]
List matrix_stats_cpp(SEXP matrix, int row_stats, int col_stats) {
    StatsResult res = run_with_R_interrupt_check(
        &MatrixLoader<double>::computeMatrixStats, peek_unique_xptr<MatrixLoader<double>>(matrix),
        (Stats)row_stats,
        (Stats)col_stats
    );

    return List::create(Named("row_stats") = res.row_stats, Named("col_stats") = res.col_stats);
}

// [[Rcpp::export]]
Eigen::MatrixXd wilcoxon_rank_sum_pval_uint32_t_cpp(SEXP matrix, std::vector<uint32_t> groups) {
    return wilcoxon_rank_sum(take_unique_xptr<MatrixLoader<uint32_t>>(matrix), groups);
}

// [[Rcpp::export]]
Eigen::MatrixXd wilcoxon_rank_sum_pval_float_cpp(SEXP matrix, std::vector<uint32_t> groups) {
    return wilcoxon_rank_sum(take_unique_xptr<MatrixLoader<float>>(matrix), groups);
}

// [[Rcpp::export]]
Eigen::MatrixXd wilcoxon_rank_sum_pval_double_cpp(SEXP matrix, std::vector<uint32_t> groups) {
    return wilcoxon_rank_sum(take_unique_xptr<MatrixLoader<double>>(matrix), groups);
}

// [[Rcpp::export]]
bool matrix_identical_uint32_t_cpp(SEXP mat1, SEXP mat2) {
    MatrixIterator<uint32_t> i1(
        std::unique_ptr<MatrixLoader<uint32_t>>(peek_unique_xptr<MatrixLoader<uint32_t>>(mat1))
    );
    MatrixIterator<uint32_t> i2(
        std::unique_ptr<MatrixLoader<uint32_t>>(peek_unique_xptr<MatrixLoader<uint32_t>>(mat2))
    );
    i1.preserve_input_loader();
    i2.preserve_input_loader();
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
        while (true) {
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
