// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once
#include "../utils/radix_sort.h"
#include "MatrixIterator.h"

namespace BPCells {

namespace {

// Generic subsetting over
template <int Rows, int Cols>
Eigen::Matrix<double, Rows, Cols> subset_eigen_rows(
    const Eigen::Map<Eigen::Matrix<double, Rows, Cols>> m, const std::vector<uint32_t> &indices
) {
    Eigen::Matrix<double, Rows, Cols> ret(indices.size(), m.cols());
    for (size_t i = 0; i < indices.size(); i++) {
        ret.row(i) = m.row(indices[i]);
    }
    return ret;
}

template <int Rows, int Cols>
Eigen::Matrix<double, Rows, Cols> subset_eigen_cols(
    const Eigen::Map<Eigen::Matrix<double, Rows, Cols>> m, const std::vector<uint32_t> &indices
) {
    Eigen::Matrix<double, Rows, Cols> ret(m.rows(), indices.size());
    for (size_t i = 0; i < indices.size(); i++) {
        ret.col(i) = m.col(indices[i]);
    }
    return ret;
}

// Non-map variants
template <int Rows, int Cols>
Eigen::Matrix<double, Rows, Cols>
subset_eigen_rows(const Eigen::Matrix<double, Rows, Cols> m, const std::vector<uint32_t> &indices) {
    Eigen::Matrix<double, Rows, Cols> ret(indices.size(), m.cols());
    for (size_t i = 0; i < indices.size(); i++) {
        ret.row(i) = m.row(indices[i]);
    }
    return ret;
}

template <int Rows, int Cols>
Eigen::Matrix<double, Rows, Cols>
subset_eigen_cols(const Eigen::Matrix<double, Rows, Cols> m, const std::vector<uint32_t> &indices) {
    Eigen::Matrix<double, Rows, Cols> ret(m.rows(), indices.size());
    for (size_t i = 0; i < indices.size(); i++) {
        ret.col(i) = m.col(indices[i]);
    }
    return ret;
}

} // namespace

// Select specific columns from a dataset
template <class T> class MatrixColSelect : public MatrixLoaderWrapper<T> {
  private:
    uint32_t current_col = UINT32_MAX;
    const std::vector<uint32_t> col_indices;
    bool is_reorder = false;
    std::vector<uint32_t> reverse_indices; // Only populated if is_reorder == true;
                                           // reverse_indices[col_indices[i]] = i

  public:
    // col_indices -- vector of columns to select
    MatrixColSelect(
        std::unique_ptr<MatrixLoader<T>> &&loader, const std::vector<uint32_t> col_indices
    )
        : MatrixLoaderWrapper<T>(std::move(loader))
        , col_indices(col_indices) {

        if (col_indices.size() == this->loader->cols()) {
            std::vector<uint8_t> seen(col_indices.size(), 0);
            is_reorder = true;
            for (const auto &c : col_indices) {
                if (seen[c]) {
                    is_reorder = false;
                    break;
                }
                seen[c] = 1;
            }
        }
        if (is_reorder) {
            reverse_indices.resize(col_indices.size(), UINT32_MAX);
            for (uint32_t i = 0; i < col_indices.size(); i++) {
                if (reverse_indices[col_indices[i]] != UINT32_MAX)
                    throw std::runtime_error("Error constructing reverse_indices in MatrixColSelect"
                    );
                reverse_indices[col_indices[i]] = i;
            }
        }
    }

    ~MatrixColSelect() = default;

    uint32_t cols() const override { return col_indices.size(); }
    const char *colNames(uint32_t col) override {
        if (col < col_indices.size()) return this->loader->colNames(col_indices[col]);
        return NULL;
    }

    void restart() override {
        current_col = UINT32_MAX;
        this->loader->restart();
    }
    void seekCol(uint32_t col) override {
        if (col >= col_indices.size())
            throw std::runtime_error("Requested column is greater than number of columns");
        this->loader->seekCol(col_indices[col]);
        current_col = col;
    }

    bool nextCol() override {
        current_col += 1;
        if (current_col >= col_indices.size()) {
            current_col -= 1;
            return false;
        } else if (current_col > 0 && col_indices[current_col - 1] == col_indices[current_col] - 1) {
            return this->loader->nextCol();
        } else {
            this->loader->seekCol(col_indices[current_col]);
            return true;
        }
    }
    uint32_t currentCol() const override { return current_col; }

    // Override sparse-dense multiply ops
    // This proved needed in cases of subsetting/reordering rows after a concatenation.
    // It will potentially increase overall work in case of a small subset, so the
    // long-term fix is to separate out permuting from filtering (since filtering
    // can be pushed down and permuting can be pulled up)
    Eigen::MatrixXd denseMultiplyRight(
        const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt = NULL
    ) override {
        if (is_reorder) {
            Eigen::MatrixXd B2 = subset_eigen_rows(B, reverse_indices);
            Eigen::Map<Eigen::MatrixXd> B2_map(B2.data(), B2.rows(), B2.cols());
            return this->loader->denseMultiplyRight(B2_map, user_interrupt);
        } else {
            return MatrixLoaderWrapper<T>::denseMultiplyRight(B, user_interrupt);
        }
    }

    Eigen::MatrixXd denseMultiplyLeft(
        const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt = NULL
    ) override {
        if (is_reorder) {
            Eigen::MatrixXd R1 = this->loader->denseMultiplyLeft(B, user_interrupt);
            return subset_eigen_cols(R1, col_indices);
        } else {
            return MatrixLoaderWrapper<T>::denseMultiplyLeft(B, user_interrupt);
        }
    }

    Eigen::VectorXd vecMultiplyRight(
        const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt = NULL
    ) override {
        if (is_reorder) {
            Eigen::VectorXd v2 = subset_eigen_rows(v, reverse_indices);
            Eigen::Map<Eigen::VectorXd> v2_map(v2.data(), v2.rows());
            return this->loader->vecMultiplyRight(v2_map, user_interrupt);
        } else {
            return MatrixLoaderWrapper<T>::vecMultiplyRight(v, user_interrupt);
        }
    }
    Eigen::VectorXd vecMultiplyLeft(
        const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt = NULL
    ) override {
        if (is_reorder) {
            Eigen::VectorXd r1 = this->loader->vecMultiplyLeft(v, user_interrupt);
            return subset_eigen_rows(r1, col_indices);
        } else {
            return MatrixLoaderWrapper<T>::vecMultiplyLeft(v, user_interrupt);
        }
    }
};

// Select range of columns from a dataset
template <class T> class MatrixColSlice : public MatrixLoaderWrapper<T> {
  private:
    uint32_t start_col, end_col;

  public:
    // [start_col, end_col) are the range of columns that will be output
    MatrixColSlice(
        std::unique_ptr<MatrixLoader<T>> &&loader, uint32_t start_col, uint32_t end_col
    )
        : MatrixLoaderWrapper<T>(std::move(loader))
        , start_col(start_col)
        , end_col(end_col) {

        if (end_col <= start_col) {
            throw std::runtime_error("MatrixColSlice: end_col must be > start_col");
        }
        if (end_col > this->loader->cols()) {
            throw std::runtime_error("MatrixColSlice: end_col must be <= loader.cols()");
        }
    }

    ~MatrixColSlice() = default;

    
    uint32_t cols() const override { return end_col - start_col; }
    const char *colNames(uint32_t col) override {
        if (col < cols()) return this->loader->colNames(col + start_col);
        return NULL;
    }

    void seekCol(uint32_t col) override {
        if (col >= cols())
            throw std::runtime_error("Requested column is greater than number of columns");
        this->loader->seekCol(col + start_col);
    }

    bool nextCol() override {
        if (!this->loader->nextCol()) return false;
        if (this->loader->currentCol() >= end_col) return false;
        if (this->loader->currentCol() < start_col) {
            this->loader->seekCol(start_col);
        }
        return true;
    }
    uint32_t currentCol() const override { return this->loader->currentCol() - start_col; }
};

// Select specific rows from a dataset
// Note: unlike the MatrixColSelect, MatrixRowSelect does not support duplicating rows,
// only filtering and/or reordering
template <class T> class MatrixRowSelect : public MatrixLoaderWrapper<T> {
  private:
    uint32_t loaded = 0;

    // Reverse lookup for row indices -- reverse_indices[i] gives the output row_id
    // for input row_id i
    std::vector<uint32_t> reverse_indices;
    const std::vector<uint32_t> row_indices;
    bool is_reorder;

  public:
    // cell_names -- vector with length <= the number of chromosomes in the input
    //     FragmentLoader-> The output cell `i` will come from input cell
    //     `chr_assignments[i]`. The entries of cell_names must be unique
    MatrixRowSelect(
        std::unique_ptr<MatrixLoader<T>> &&loader, const std::vector<uint32_t> row_indices
    )
        : MatrixLoaderWrapper<T>(std::move(loader))
        , reverse_indices(this->loader->rows(), UINT32_MAX)
        , row_indices(row_indices) {
        for (uint32_t i = 0; i < row_indices.size(); i++) {
            if (reverse_indices[row_indices[i]] != UINT32_MAX)
                throw std::runtime_error("Cannot duplicate rows using MatrixRowSelect");
            reverse_indices[row_indices[i]] = i;
        }
        is_reorder = row_indices.size() == this->loader->rows();
    }

    ~MatrixRowSelect() = default;

    uint32_t rows() const override { return row_indices.size(); }
    const char *rowNames(uint32_t row) override {
        if (row < row_indices.size()) return this->loader->rowNames(row_indices[row]);
        return NULL;
    }

    void restart() override {
        loaded = 0;
        this->loader->restart();
    }
    void seekCol(uint32_t col) override {
        loaded = 0;
        this->loader->seekCol(col);
    }

    bool nextCol() override {
        loaded = 0;
        return this->loader->nextCol();
    }

    bool load() override {
        // Just perform a straight filter and load incrementally
        loaded = 0;
        while (loaded == 0) {
            if (!this->loader->load()) return false;

            uint32_t cap = this->loader->capacity();
            uint32_t *row_ptr = this->loader->rowData();
            T *val_ptr = this->loader->valData();

            for (uint32_t i = 0; i < cap; i++) {
                uint32_t new_row = reverse_indices[row_ptr[i]];
                row_ptr[loaded] = new_row;
                val_ptr[loaded] = val_ptr[i];
                loaded += new_row != UINT32_MAX;
            }
        }
        return true;
    }

    uint32_t capacity() const override { return loaded; }
    uint32_t *rowData() override { return this->loader->rowData(); }
    T *valData() override { return this->loader->valData(); }

    // Override sparse-dense multiply ops
    // This proved needed in cases of subsetting/reordering rows after a concatenation.
    // It will potentially increase overall work in case of a small subset, so the
    // long-term fix is to separate out permuting from filtering (since filtering
    // can be pushed down and permuting can be pulled up)
    Eigen::MatrixXd denseMultiplyRight(
        const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt = NULL
    ) override {
        if (is_reorder) {
            Eigen::MatrixXd R1 = this->loader->denseMultiplyRight(B, user_interrupt);
            return subset_eigen_rows(R1, row_indices);
        } else {
            return MatrixLoaderWrapper<T>::denseMultiplyRight(B, user_interrupt);
        }
    }

    Eigen::MatrixXd denseMultiplyLeft(
        const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt = NULL
    ) override {
        if (is_reorder) {
            Eigen::MatrixXd B2 = subset_eigen_cols(B, reverse_indices);
            Eigen::Map<Eigen::MatrixXd> B2_map(B2.data(), B2.rows(), B2.cols());
            return this->loader->denseMultiplyLeft(B2_map, user_interrupt);
        } else {
            return MatrixLoaderWrapper<T>::denseMultiplyLeft(B, user_interrupt);
        }
    }

    Eigen::VectorXd vecMultiplyRight(
        const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt = NULL
    ) override {
        if (is_reorder) {
            Eigen::VectorXd r1 = this->loader->vecMultiplyRight(v, user_interrupt);
            return subset_eigen_rows(r1, row_indices);
        } else {
            return MatrixLoaderWrapper<T>::vecMultiplyRight(v, user_interrupt);
        }
    }

    Eigen::VectorXd vecMultiplyLeft(
        const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt = NULL
    ) override {
        if (is_reorder) {
            Eigen::VectorXd v2 = subset_eigen_rows(v, reverse_indices);
            Eigen::Map<Eigen::VectorXd> v2_map(v2.data(), v2.rows());
            return this->loader->vecMultiplyLeft(v2_map, user_interrupt);
        } else {
            return MatrixLoaderWrapper<T>::vecMultiplyLeft(v, user_interrupt);
        }
    }
};

} // end namespace BPCells