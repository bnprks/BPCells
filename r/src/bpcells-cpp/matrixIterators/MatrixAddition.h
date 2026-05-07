// Copyright 2026 BPCells contributors
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include "MatrixIterator.h"
#include "OrderRows.h"

namespace BPCells {

// Element-wise addition of two input matrices with matching dimensions via 2-way merge.
template <typename T> class MatrixAddition : public MatrixLoader<T> {
  protected:
    // Running invariants:
    // - left and right are always about to read the next column for output
    // - left and right are always positioned at the same column (lockstep)
    // - [*_idx, *_cap) are the unconsumed left/right buffer ranges
    // - *_done are true when the column is exhausted
    // - OrderRows guarantees that row indices are strictly increasing
    OrderRows<T> left, right;
    std::vector<T> val_buffer;
    std::vector<uint32_t> row_buffer;
    uint32_t loaded = 0;
    uint32_t max_load_size;
    uint32_t l_idx = 0, l_cap = 0;
    uint32_t r_idx = 0, r_cap = 0;
    bool l_done = true, r_done = true;

    bool refill(OrderRows<T> &side, uint32_t &idx, uint32_t &cap, bool &done) {
        idx = 0;
        if (side.load()) {
            cap = side.capacity();
            return true;
        }
        cap = 0;
        done = true;
        return false;
    }

  public:
    MatrixAddition(std::unique_ptr<MatrixLoader<T>> &&left, std::unique_ptr<MatrixLoader<T>> &&right, uint32_t load_size = 1024)
        : left(std::move(left))
        , right(std::move(right))
        , max_load_size(load_size) {

        if (this->left.rows() != this->right.rows() ||
            this->left.cols() != this->right.cols())
            throw std::runtime_error("Matrices have incompatible dimensions for element-wise addition");

        row_buffer.resize(load_size);
        val_buffer.resize(load_size);
    }

    uint32_t rows() const override { return left.rows(); }
    uint32_t cols() const override { return left.cols(); }

    const char *rowNames(uint32_t row) override { return left.rowNames(row); }
    const char *colNames(uint32_t col) override { return left.colNames(col); }

    // Reset the iterators to start from the beginning
    void restart() override {
        left.restart();
        right.restart();
        loaded = 0;
        l_idx = l_cap = 0;
        r_idx = r_cap = 0;
        l_done = r_done = true; // Stays done till nextCol/seekCol
    }

    // Seek to a specific column without reading data
    // Next call should be to load(). col must be < cols()
    void seekCol(uint32_t col) override {
        left.seekCol(col);
        right.seekCol(col);
        loaded = 0;
        l_idx = l_cap = 0;
        r_idx = r_cap = 0;
        l_done = r_done = false;
    }

    // Advance both matrices to the next column,
    // return false if there are no more columns,
    // or error if columns mismatch
    bool nextCol() override {
        bool l = left.nextCol();
        bool r = right.nextCol();
        if (!l && !r) return false;
        if (l != r) throw std::runtime_error("MatrixAddition: left and right nextCol() differs. Matrices desynchronized");
        loaded = 0;
        l_idx = l_cap = 0;
        r_idx = r_cap = 0;
        l_done = r_done = false;
        return true;
    }

    // Return the index of the current column
    uint32_t currentCol() const override { return left.currentCol(); }

    // Return false if there are no more entries to load
    bool load() override {
        loaded = 0;

        // Refill column buffers
        if (!l_done && l_idx >= l_cap) refill(left, l_idx, l_cap, l_done);
        if (!r_done && r_idx >= r_cap) refill(right, r_idx, r_cap, r_done);

        while (loaded < max_load_size) {
            bool l_has = !l_done && l_idx < l_cap;
            bool r_has = !r_done && r_idx < r_cap;
            if (!l_has && !r_has) break;

            if (!r_has) {
                // Drain left in bulk
                uint32_t take = std::min<uint32_t>(max_load_size - loaded, l_cap - l_idx);
                const uint32_t *lr = left.rowData();
                const T *lv = left.valData();
                std::copy(lr + l_idx, lr + l_idx + take, row_buffer.data() + loaded);
                std::copy(lv + l_idx, lv + l_idx + take, val_buffer.data() + loaded);
                l_idx += take;
                loaded += take;
                if (l_idx >= l_cap) refill(left, l_idx, l_cap, l_done);
            } else if (!l_has) {
                // Drain right in bulk
                uint32_t take = std::min<uint32_t>(max_load_size - loaded, r_cap - r_idx);
                const uint32_t *rr = right.rowData();
                const T *rv = right.valData();
                std::copy(rr + r_idx, rr + r_idx + take, row_buffer.data() + loaded);
                std::copy(rv + r_idx, rv + r_idx + take, val_buffer.data() + loaded);
                r_idx += take;
                loaded += take;
                if (r_idx >= r_cap) refill(right, r_idx, r_cap, r_done);
            } else {
                const uint32_t *lr = left.rowData();
                const uint32_t *rr = right.rowData();
                const T *lv = left.valData();
                const T *rv = right.valData();
                while (loaded < max_load_size && l_idx < l_cap && r_idx < r_cap) {
                    uint32_t lr_v = lr[l_idx];
                    uint32_t rr_v = rr[r_idx];
                    T lv_v = lv[l_idx];
                    T rv_v = rv[r_idx];

                    bool l_pick = lr_v <= rr_v;
                    bool r_pick = lr_v >= rr_v;
                    uint32_t out_row = l_pick ? lr_v : rr_v;
                    T out_val = l_pick ? (r_pick ? lv_v + rv_v : lv_v) : rv_v;

                    bool emit = out_val != T(0);
                    row_buffer[loaded] = out_row;
                    val_buffer[loaded] = out_val;
                    loaded += emit;
                    l_idx += l_pick;
                    r_idx += r_pick;
                }
                if (!l_done && l_idx >= l_cap) refill(left, l_idx, l_cap, l_done);
                if (!r_done && r_idx >= r_cap) refill(right, r_idx, r_cap, r_done);
            }
        }

        return loaded > 0;
    }

    // Number of loaded entries available
    uint32_t capacity() const override { return loaded; }

    // Pointers to the loaded entries
    uint32_t *rowData() override { return row_buffer.data(); }
    T *valData() override { return val_buffer.data(); }

    // MATH OPERATIONS: utilize distributive property that (A+B)*v = A*v + B*v
    Eigen::VectorXd vecMultiplyRight(
        const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt = NULL
    ) override {
        return left.vecMultiplyRight(v, user_interrupt) +
               right.vecMultiplyRight(v, user_interrupt);
    }
    Eigen::VectorXd vecMultiplyLeft(
        const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt = NULL
    ) override {
        return left.vecMultiplyLeft(v, user_interrupt) +
               right.vecMultiplyLeft(v, user_interrupt);
    }
    Eigen::MatrixXd denseMultiplyRight(
        const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt = NULL
    ) override {
        return left.denseMultiplyRight(B, user_interrupt) +
               right.denseMultiplyRight(B, user_interrupt);
    }
    Eigen::MatrixXd denseMultiplyLeft(
        const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt = NULL
    ) override {
        return left.denseMultiplyLeft(B, user_interrupt) +
               right.denseMultiplyLeft(B, user_interrupt);
    }
    // Linearity of sums (distribution over addition):
    //     rowSums(A+B) = rowSums(A) + rowSums(B) 
    //     colSums(A+B) = colSums(A) + colSums(B) 
    std::vector<T> colSums(std::atomic<bool> *user_interrupt = NULL) override {
        auto l = left.colSums(user_interrupt);
        auto r = right.colSums(user_interrupt);
        for (uint32_t i = 0; i < l.size(); i++) {
            l[i] += r[i];
        }
        return l;
    }
    std::vector<T> rowSums(std::atomic<bool> *user_interrupt = NULL) override {
        auto l = left.rowSums(user_interrupt);
        auto r = right.rowSums(user_interrupt);
        for (uint32_t i = 0; i < l.size(); i++) {
            l[i] += r[i];
        }
        return l;
    }
};

} // end namespace BPCells
