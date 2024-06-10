// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once
#include "MatrixIterator.h"

namespace BPCells {

// Rename rows/cols in a matrix
template <class T> class RenameDims : public MatrixLoaderWrapper<T> {
  private:
    const std::vector<std::string> row_names, col_names;
    const bool clear_row_names;
    const bool clear_col_names;

  public:
    // To chanage row names, provide row_names and clear_row_names=false
    // To preserve row names, make row_names length 0 and clear_row_names=false
    // To clear row names, make clear_row_names=true
    // Same goes for cols
    RenameDims(
        std::unique_ptr<MatrixLoader<T>> &&loader, const std::vector<std::string> &row_names, const std::vector<std::string> &col_names,
        bool clear_row_names = false, bool clear_col_names = false
    )
        : MatrixLoaderWrapper<T>(std::move(loader))
        , row_names(row_names)
        , col_names(col_names)
        , clear_row_names(clear_row_names)
        , clear_col_names(clear_col_names) {

        if (row_names.size() > 0 && row_names.size() != this->loader->rows()) {
            throw std::runtime_error("RenameDims: Row names must be length 0 or equal to number of input rows");
        }

        if (col_names.size() > 0 && col_names.size() != this->loader->cols()) {
            throw std::runtime_error("RenameDims: Col names must be length 0 or equal to number of input cols");
        }

        if (clear_row_names && row_names.size() != 0) {
            throw std::runtime_error("RenameDims: if clear_row_names is true, row names must be length 0");
        }

        if (clear_col_names && col_names.size() != 0) {
            throw std::runtime_error("RenameDims: if clear_col_names is true, col names must be length 0");
        }
    }

    ~RenameDims() = default;

    const char *colNames(uint32_t col) override {
        if (clear_col_names) return NULL;
        if (col_names.size() == 0) return this->loader->colNames(col);
        if (col < col_names.size()) return col_names[col].c_str();
        return NULL;
    }

    const char *rowNames(uint32_t row) override {
        if (clear_row_names) return NULL;
        if (row_names.size() == 0) return this->loader->rowNames(row);
        if (row < row_names.size()) return row_names[row].c_str();
        return NULL;
    }

    // Defer through all the math ops to tihe input matrix
    Eigen::MatrixXd denseMultiplyRight(
        const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt = NULL
    ) override {
        return this->loader->denseMultiplyRight(B, user_interrupt);
    }
    Eigen::MatrixXd denseMultiplyLeft(
        const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt = NULL
    ) override {
        return this->loader->denseMultiplyLeft(B, user_interrupt);
    }
    Eigen::VectorXd vecMultiplyRight(
        const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt = NULL
    ) override {
        return this->loader->vecMultiplyRight(v, user_interrupt);
    }
    Eigen::VectorXd vecMultiplyLeft(
        const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt = NULL
    ) override {
        return this->loader->vecMultiplyLeft(v, user_interrupt);
    }
    std::vector<T> colSums(std::atomic<bool> *user_interrupt = NULL) override {
        return this->loader->colSums();
    }
    std::vector<T> rowSums(std::atomic<bool> *user_interrupt = NULL) override {
        return this->loader->rowSums();
    }
    StatsResult computeMatrixStats(Stats row_stats, Stats col_stats, std::atomic<bool> *user_interrupt = NULL) override {
        return this->loader->computeMatrixStats(row_stats, col_stats, user_interrupt);
    }

};

} // namespace BPCells