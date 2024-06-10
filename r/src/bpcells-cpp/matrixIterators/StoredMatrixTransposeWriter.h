// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <atomic>
#include <string>

#include "../arrayIO/array_interfaces.h"
#include "../arrayIO/binaryfile.h"
#include "../arrayIO/vector.h"
#include "../../vendor/dary_heap/dary_heap.hpp"
#include "../utils/radix_sort.h"
#include "MatrixIterator.h"
#include "StoredMatrix.h"
#include "StoredMatrixSorter.h"

namespace BPCells {

// Transpose the elements of an input column-major matrix and save a stored
// matrix on disk. Optionally, if saving as row-major the matrix will be logically
// the same as the input, though the storage order will physically be transposed.
template <typename T> class StoredMatrixTransposeWriter : public StoredMatrixSorter<T>, public MatrixWriter<T> {
    // Number of loaded matrix value that have been previously output via load_entries
    uint64_t previously_loaded = 0;
    MatrixLoader<T> *mat = NULL;

    size_t load_entries(
        std::vector<uint32_t> &row, std::vector<uint32_t> &col, std::vector<T> &val, std::atomic<bool> *user_interrupt
    ) override {
        uint64_t loaded = 0;
        while (loaded < row.size()) {
            if (user_interrupt != NULL && *user_interrupt) return loaded;
            // Load data (or re-use leftover data)
            if (previously_loaded == 0 && !mat->load()) {
                if (!mat->nextCol()) break;
                continue;
            }
            // Calculate the number of values to copy to outputs
            size_t available = mat->capacity() - previously_loaded;
            uint32_t copy_count = std::min<size_t>(row.size() - loaded, available);
            
            // Copy the values
            std::memmove(
                col.data() + loaded,
                mat->rowData() + previously_loaded,
                copy_count * sizeof(uint32_t)
            );
            std::memmove(
                val.data() + loaded, mat->valData() + previously_loaded, copy_count * sizeof(T)
            );
            uint32_t current_col = mat->currentCol();
            for (uint32_t i = 0; i < copy_count; i++) {
                row[i + loaded] = current_col;
            }
            loaded += copy_count;

            previously_loaded = available != copy_count ? previously_loaded + copy_count : 0; 
        }
        return loaded;
    }
    bool row_sorted() const override { return true; }

  public:
    using StoredMatrixSorter<T>::StoredMatrixSorter;

    void write(MatrixLoader<T> &mat, std::atomic<bool> *user_interrupt = NULL) override {
        this->mat = &mat;
        mat.restart();
        // Store row and col names. This probably incurs a few extra copies,
        // but it shouldn't matter since writing the actual matrix should dominate cost
        std::vector<std::string> col_names(0);
        for (int i = 0;; i++) {
            const char *col_name = mat.rowNames(i);
            if (col_name == NULL) break;
            col_names.push_back(std::string(col_name));
        }

        std::vector<std::string> row_names(0);
        for (int i = 0;; i++) {
            const char *row_name = mat.colNames(i);
            if (row_name == NULL) break;
            row_names.push_back(std::string(row_name));
        }

        mat.nextCol();
        this->writeValues(std::move(row_names), std::move(col_names), mat.cols(), mat.rows(), user_interrupt);
    }
};

} // end namespace BPCells