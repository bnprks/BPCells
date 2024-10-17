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
#include <unordered_map>

namespace BPCells {

// Re-order a MatixIterator so that each column gives entries sorted by row index
template <class T> class OrderRows : public MatrixLoaderWrapper<T> {
    friend class MatrixTransformDense;

  private:
    std::vector<uint32_t> row_data, row_buf;
    std::vector<T> val_data, val_buf;

    uint32_t idx = 0;
    uint32_t cap = 0;
    uint32_t load_size;

  public:
    OrderRows(std::unique_ptr<MatrixLoader<T>> &&loader, uint32_t load_size = 1024)
        : MatrixLoaderWrapper<T>(std::move(loader))
        , load_size(load_size) {
        row_data.resize(this->loader->rows());
        row_buf.resize(this->loader->rows());
        val_data.resize(this->loader->rows());
        val_buf.resize(this->loader->rows());
    }

    ~OrderRows() = default;

    void restart() override {
        idx = 0;
        cap = 0;
        this->loader->restart();
    }
    void seekCol(uint32_t col) override {
        idx = 0;
        cap = 0;
        this->loader->seekCol(col);
    }

    bool nextCol() override {
        idx = 0;
        cap = 0;
        return this->loader->nextCol();
    }

    bool load() override {
        if (idx == 0 && cap == 0) {
            bool needs_reorder = false;
            uint32_t prev_row = 0;
            // Load all the values for the column and sort
            while (this->loader->load()) {
                uint32_t loaded = this->loader->capacity();

                uint32_t *row_ptr = this->loader->rowData();
                T *val_ptr = this->loader->valData();

                std::memmove(val_data.data() + cap, val_ptr, sizeof(T) * loaded);

                for (uint32_t i = 0; i < loaded; i++) {
                    row_data[cap + i] = row_ptr[i];
                    // Assume no duplicate row indices so we don't need <=
                    if (row_ptr[i] < prev_row) needs_reorder = true;
                    prev_row = row_ptr[i];
                }
                cap += loaded;
            }
            if (needs_reorder) {
                // Sort the entries by increasing row ID
                lsdRadixSortArrays<uint32_t, T>(cap, row_data, val_data, row_buf, val_buf);
            }
        } else {
            idx += std::min(cap - idx, load_size);
        }
        return idx < cap;
    }

    uint32_t capacity() const override { return std::min(cap - idx, load_size); }
    uint32_t *rowData() override { return row_data.data() + idx; }
    T *valData() override { return val_data.data() + idx; }

    // Override the mat-vec ops, since sorting is not necessary for these and it showed up
    // prominently in a profile I ran (patterns from a Scale+Shift matrix with a subset assignment.
    // The input to Scale got row-reordered in the process)
    Eigen::MatrixXd
    denseMultiplyRight(const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt) override {
        return this->loader->denseMultiplyRight(B, user_interrupt);
    }
    Eigen::MatrixXd
    denseMultiplyLeft(const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt) override {
        return this->loader->denseMultiplyLeft(B, user_interrupt);
    }
    // Calculate matrix-vector product A*v where A=this and B is a dense matrix.
    Eigen::VectorXd
    vecMultiplyRight(const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt) override {
        return this->loader->vecMultiplyRight(v, user_interrupt);
    }
    Eigen::VectorXd
    vecMultiplyLeft(const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt) override {
        return this->loader->vecMultiplyLeft(v, user_interrupt);
    }
};

} // end namespace BPCells
