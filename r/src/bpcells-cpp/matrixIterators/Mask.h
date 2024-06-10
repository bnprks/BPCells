// Copyright 2023 BPCells contributors
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

// Mask the elements from a matrix.
// If Invert == false, for non-zero elements in the mask matrix, set to zero.
// If Invert == true, for zero elements in mask matrix, set mat to zero.
template <typename T, bool Invert = false> class Mask : public MatrixLoader<T> {
  protected:
    OrderRows<T> mat;
    OrderRows<uint32_t> mask;
    uint32_t mask_idx = UINT32_MAX;
    uint32_t loaded;
  public:
    Mask(std::unique_ptr<MatrixLoader<T>> &&mat, std::unique_ptr<MatrixLoader<uint32_t>> &&mask)
        : mat(std::move(mat))
        , mask(std::move(mask), this->mat.rows()) {

        if (this->mat.cols() != this->mask.cols() || this->mat.rows() != this->mask.rows())
            throw std::runtime_error("Matrices have mismatched dimensions");
    }

    uint32_t rows() const override { return mat.rows(); }
    uint32_t cols() const override { return mat.cols(); }

    const char *rowNames(uint32_t row) override { return mat.rowNames(row); }
    const char *colNames(uint32_t col) override { return mat.colNames(col); }

    // Reset the iterator to start from the beginning
    void restart() override {
        mat.restart();
        mask.restart();
        mask_idx = UINT32_MAX;
    }

    // Seek to a specific column without reading data
    // Next call should be to load(). col must be < cols()
    void seekCol(uint32_t col) override {
        mat.seekCol(col);
        mask.seekCol(col);
        mask_idx = UINT32_MAX;
    }

    // Advance to the next column, or return false if there
    // are no more columns
    bool nextCol() override {
        if (!mat.nextCol()) return false;
        if(!mask.nextCol())
            throw std::runtime_error("ElementwiseMultiply: Unexpected failure to read next column");
        mask_idx = UINT32_MAX;
        return true;
    }

    // Return the index of the current column
    uint32_t currentCol() const override { return mat.currentCol(); }

    // Return false if there are no more entries to load
    bool load() override {
        if (mask_idx == UINT32_MAX) {
            if (!mask.load()) {
                if constexpr (Invert) {
                    return false;
                }
                if constexpr (!Invert) {
                    if (!mat.load()) return false;
                    loaded = mat.capacity();
                    return true;
                }
            }
            mask_idx = 0;
        }
        uint32_t *mask_row = mask.rowData();
        uint32_t mask_cap = mask.capacity();

        loaded = 0;
        while (loaded == 0) {
            if (!mat.load()) return false;
            uint32_t *mat_row = mat.rowData();
            T *mat_val = mat.valData();
            uint32_t cap = mat.capacity();
        
            for (uint32_t i = 0; i < cap; i++) {
                while (mask_idx < mask_cap && mask_row[mask_idx] < mat_row[i]) mask_idx++;
                if (mask_idx >= mask_cap) {
                    if constexpr (Invert) return loaded > 0; // No more can be loaded
                    if constexpr (!Invert) {
                        // Copy the remaining values directly 
                        std::memmove(mat_row + loaded, mat_row + i, sizeof(uint32_t) * (cap - i));
                        std::memmove(mat_val + loaded, mat_val + i, sizeof(T) * (cap - i));
                        loaded += cap - i;
                        break;
                    }
                }
                mat_row[loaded] = mat_row[i];
                mat_val[loaded] = mat_val[i];
                if constexpr (!Invert) {
                    loaded += mask_row[mask_idx] != mat_row[i];
                }
                if constexpr (Invert) {
                    loaded += mask_row[mask_idx] == mat_row[i];
                }
            }
        }
        
        return true;
    }

    // Number of loaded entries available
    uint32_t capacity() const override { return loaded; }

    // Pointers to the loaded entries
    uint32_t *rowData() override { return mat.rowData(); }
    T *valData() override { return mat.valData(); }
};

} // end namespace BPCells