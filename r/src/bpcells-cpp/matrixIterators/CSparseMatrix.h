// Copyright 2021 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

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
template<typename T>
class CSparseMatrix : public MatrixLoader<T> {
    typedef Eigen::Map<Eigen::SparseMatrix<T>> EigenMat;

  private:
    const EigenMat mat;
    std::vector<uint32_t> row_buf;
    std::vector<T> val_buf;
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
        std::memmove(val_buf.data(), mat.valuePtr() + idx, sizeof(T) * num_loaded);
        return true;
    };

    uint32_t capacity() const override { return num_loaded; }

    uint32_t *rowData() override { return row_buf.data(); }
    T *valData() override { return val_buf.data(); }
};

template<typename T>
class CSparseMatrixWriter : public MatrixWriter<T> {
  private:
    Eigen::SparseMatrix<T> eigen_mat;

  public:
    void write(MatrixLoader<T> &loader, std::atomic<bool> *user_interrupt = NULL) override {
        MatrixIterator<T> mat((std::unique_ptr<MatrixLoader<T>>(&loader)));
        // Don't take ownership of our input loader
        mat.preserve_input_loader();

        uint32_t count = 0;
        std::vector<Eigen::Triplet<T>> triplets;

        while (mat.nextCol()) {
            while (mat.nextValue()) {
                triplets.push_back(Eigen::Triplet<T>(mat.row(), mat.col(), mat.val()));
                if (count++ % 8192 == 0 && user_interrupt != NULL && *user_interrupt) return;
            }
        }

        eigen_mat = Eigen::SparseMatrix<T>(mat.rows(), mat.cols());
        eigen_mat.setFromTriplets(triplets.begin(), triplets.end());
    };

    const Eigen::SparseMatrix<T> getMat() { return eigen_mat; }
};

template<typename T>
class CSparseTransposeMatrixWriter : public MatrixWriter<T> {
  private:
    Eigen::SparseMatrix<T> eigen_mat;

  public:
    void write(MatrixLoader<T> &loader, std::atomic<bool> *user_interrupt = NULL) override {
        MatrixIterator<T> mat((std::unique_ptr<MatrixLoader<T>>(&loader)));
        // Don't take ownership of our input loader
        mat.preserve_input_loader();
        uint32_t count = 0;
        std::vector<Eigen::Triplet<T>> triplets;

        while (mat.nextCol()) {
            while (mat.nextValue()) {
                triplets.push_back(Eigen::Triplet<T>(mat.col(), mat.row(), mat.val()));
                if (count++ % 8192 == 0 && user_interrupt != NULL && *user_interrupt) return;
            }
        }

        eigen_mat = Eigen::SparseMatrix<T>(mat.cols(), mat.rows());
        eigen_mat.setFromTriplets(triplets.begin(), triplets.end());
    };

    const Eigen::SparseMatrix<T> getMat() { return eigen_mat; }
};

} // end namespace BPCells