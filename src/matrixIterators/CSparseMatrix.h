#include <atomic>

#include "../arrayIO/array_interfaces.h"
#include "MatrixIterator.h"

#ifndef RCPP_EIGEN
#include <Eigen/SparseCore>
#else
#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <RcppEigen.h>
#endif
// [[Rcpp::depends(RcppEigen)]]

namespace BPCells {

// Get Eigen and iterate over an Eigen sparse matrix
class CSparseMatrix : public MatrixLoader<double> {
    typedef Eigen::Map<Eigen::SparseMatrix<double>> EigenMat;

  private:
    const EigenMat mat;
    std::vector<uint32_t> row_buf;
    std::vector<double> val_buf;
    std::unique_ptr<StringReader> row_names, col_names;
    uint32_t idx;
    uint32_t load_size;
    uint32_t num_loaded = 0;
    uint32_t col;

  public:
    CSparseMatrix(
        const EigenMat mat,
        std::unique_ptr<StringReader> &&row_names = NULL,
        std::unique_ptr<StringReader> &&col_names = NULL,
        uint32_t load_size = 1024
    )
        : mat(mat)
        , row_names(std::move(row_names))
        , col_names(std::move(col_names))
        , load_size(load_size) {

        row_buf.resize(load_size);
        val_buf.resize(load_size);
        restart();
    }

    // Return the count of rows and columns
    uint32_t rows() const override { return mat.rows(); };
    uint32_t cols() const override { return mat.cols(); };

    // Return name for a given row or column.
    // If a matrix doesn't have assigned names this will return NULL
    const char *rowNames(uint32_t row) override {
        if (!row_names) return NULL;
        return row_names->get(row);
    }
    const char *colNames(uint32_t col) override {
        if (!col_names) return NULL;
        return col_names->get(col);
    }

    // Reset the iterator to start from the beginning
    void restart() override {
        col = UINT32_MAX;
        num_loaded = 0;
        idx = UINT32_MAX;
    };

    void seekCol(uint32_t new_col) override {
        col = new_col;
        idx = mat.outerIndexPtr()[col];
        num_loaded = 0;
    }

    bool nextCol() override {
        if (col + 1 >= mat.cols()) {
            idx = UINT32_MAX;
            return false;
        }
        col++;
        idx = mat.outerIndexPtr()[col];
        num_loaded = 0;
        return true;
    }

    uint32_t currentCol() const override { return col; }

    bool load() override {
        idx += capacity();
        if (idx >= (uint32_t)mat.outerIndexPtr()[col + 1]) {
            num_loaded = 0;
            return false;
        }
        num_loaded = std::min(load_size, mat.outerIndexPtr()[col + 1] - idx);
        std::memmove(row_buf.data(), mat.innerIndexPtr() + idx, sizeof(uint32_t) * num_loaded);
        std::memmove(val_buf.data(), mat.valuePtr() + idx, sizeof(double) * num_loaded);
        return true;
    };

    uint32_t capacity() const override { return num_loaded; }

    uint32_t *rowData() override { return row_buf.data(); }
    double *valData() override { return val_buf.data(); }
};

class CSparseMatrixWriter : public MatrixWriter<double> {
  private:
    Eigen::SparseMatrix<double> eigen_mat;

  public:
    void write(MatrixLoader<double> &loader, std::atomic<bool> *user_interrupt = NULL) override {
        MatrixIterator<double> mat((std::unique_ptr<MatrixLoader<double>>(&loader)));
        // Don't take ownership of our input loader
        mat.preserve_input_loader();

        uint32_t count = 0;
        std::vector<Eigen::Triplet<double>> triplets;

        while (mat.nextCol()) {
            while (mat.nextValue()) {
                triplets.push_back(Eigen::Triplet<double>(mat.row(), mat.col(), mat.val()));
                if (count++ % 8192 == 0 && user_interrupt != NULL && *user_interrupt) return;
            }
        }

        eigen_mat = Eigen::SparseMatrix<double>(mat.rows(), mat.cols());
        eigen_mat.setFromTriplets(triplets.begin(), triplets.end());
    };

    const Eigen::SparseMatrix<double> getMat() { return eigen_mat; }
};

class CSparseTransposeMatrixWriter : public MatrixWriter<double> {
  private:
    Eigen::SparseMatrix<double> eigen_mat;

  public:
    void write(MatrixLoader<double> &loader, std::atomic<bool> *user_interrupt = NULL) override {
        MatrixIterator<double> mat((std::unique_ptr<MatrixLoader<double>>(&loader)));
        // Don't take ownership of our input loader
        mat.preserve_input_loader();
        uint32_t count = 0;
        std::vector<Eigen::Triplet<double>> triplets;

        while (mat.nextCol()) {
            while (mat.nextValue()) {
                triplets.push_back(Eigen::Triplet<double>(mat.col(), mat.row(), mat.val()));
                if (count++ % 8192 == 0 && user_interrupt != NULL && *user_interrupt) return;
            }
        }

        eigen_mat = Eigen::SparseMatrix<double>(mat.cols(), mat.rows());
        eigen_mat.setFromTriplets(triplets.begin(), triplets.end());
    };

    const Eigen::SparseMatrix<double> getMat() { return eigen_mat; }
};

} // end namespace BPCells