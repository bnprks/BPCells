// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include "../arrayIO/array_interfaces.h"
#include "../arrayIO/bp128.h"
#include "MatrixIterator.h"

namespace BPCells {

// Main class for accessing matrices stored on disk.
// Templated to help with compatibility reading 10x and AnnData matrix formats.
// Supports uint, ulong, float, and double types. Only uint supports bitpacking compression of
// values. Row-major matrices will be loaded as their column-major transpose
template <class T> class StoredMatrix : public MatrixLoader<T> {
  private:
    NumReader<uint32_t> row;
    NumReader<T> val;
    NumReader<uint64_t> col_ptr;
    std::unique_ptr<StringReader> row_names, col_names;
    uint32_t n_rows;
    uint32_t n_cols;
    uint32_t current_col = UINT32_MAX;
    uint64_t current_idx = 0;
    uint64_t next_col_ptr;
    uint64_t current_capacity = 0;

  public:
    StoredMatrix() = default;
    StoredMatrix(StoredMatrix &&other) = default;
    StoredMatrix &operator=(StoredMatrix &&other) = default;
    StoredMatrix(
        NumReader<uint32_t> &&row,
        NumReader<T> &&val,
        NumReader<uint64_t> &&col_ptr,
        uint32_t row_count,
        std::unique_ptr<StringReader> &&row_names,
        std::unique_ptr<StringReader> &&col_names
    )
        : row(std::move(row))
        , val(std::move(val))
        , col_ptr(std::move(col_ptr))
        , row_names(std::move(row_names))
        , col_names(std::move(col_names))
        , n_rows(row_count)
        , n_cols(this->col_ptr.size() - 1)
        , next_col_ptr(this->col_ptr.read_one()) {}

    static std::string versionString(bool packed, uint32_t version) {
        std::string ret = packed ? "packed-" : "unpacked-";
        // Static check that our type is one we expect
        static_assert(
            std::disjunction_v<
                std::is_same<T, uint32_t>,
                std::is_same<T, uint64_t>,
                std::is_same<T, double>,
                std::is_same<T, float>>,
            "Type must be one of uint32_t, uint64_t, double, or float"
        );
        if constexpr (std::is_same_v<T, uint32_t>) ret += "uint-";
        if constexpr (std::is_same_v<T, uint64_t>) ret += "ulong-";
        if constexpr (std::is_same_v<T, double>) ret += "double-";
        if constexpr (std::is_same_v<T, float>) ret += "float-";
        ret += "matrix-v";
        ret += std::to_string(version);
        return ret;
    }

    // Open an unpacked StoredMatrix from a ReaderBuilder in a column-major orientation
    static StoredMatrix<T> openUnpacked(
        ReaderBuilder &rb,
        std::unique_ptr<StringReader> &&row_names,
        std::unique_ptr<StringReader> &&col_names,
        uint32_t row_count
    ) {
        ULongReader col_ptr;
        if (rb.readVersion() == versionString(false, 1)) {
            col_ptr = rb.openUIntReader("idxptr").convert<uint64_t>();
        } else if (rb.readVersion() == versionString(false, 2)) {
            col_ptr = rb.openULongReader("idxptr");
        } else {
            throw std::runtime_error(
                std::string("Version does not match ") + versionString(false, 2) + ": " +
                rb.readVersion()
            );
        }

        return StoredMatrix<T>(
            rb.openUIntReader("index"),
            rb.open<T>("val"),
            std::move(col_ptr),
            row_count,
            std::move(row_names),
            std::move(col_names)
        );
    }

    // Open an unpacked StoredMatrix from a ReaderBuilder, converting row-major orientation to
    // column-major as needed
    static StoredMatrix<T> openUnpacked(ReaderBuilder &rb) {
        auto storage_order_reader = rb.openStringReader("storage_order");
        auto storage_order = storage_order_reader->get(0);

        bool row_major = false;
        if (std::string_view("row") == storage_order) row_major = true;
        else if (std::string("col") == storage_order) row_major = false;
        else
            throw std::runtime_error(
                std::string("storage_order must be either \"row\" or \"col\", found: \"") +
                storage_order + "\""
            );

        auto row_names = rb.openStringReader("row_names");
        auto col_names = rb.openStringReader("col_names");

        auto shape = rb.openUIntReader("shape");
        uint32_t row_count = shape.read_one();
        if (row_major) {
            row_count = shape.read_one();
            std::swap(row_names, col_names);
        }

        return StoredMatrix<T>::openUnpacked(
            rb, std::move(row_names), std::move(col_names), row_count
        );
    }

    // Open a packed StoredMatrix from a ReaderBuilder in a column-major orientation
    static StoredMatrix<T> openPacked(
        ReaderBuilder &rb,
        uint32_t load_size,
        std::unique_ptr<StringReader> &&row_names,
        std::unique_ptr<StringReader> &&col_names,
        uint32_t row_count
    ) {
        ULongReader col_ptr, index_idx_offsets;
        if (rb.readVersion() == versionString(true, 1)) {
            col_ptr = rb.openUIntReader("idxptr").convert<uint64_t>();
            index_idx_offsets = ConstNumReader<uint64_t>::create({0, UINT64_MAX});
        } else if (rb.readVersion() == versionString(true, 2)) {
            col_ptr = rb.openULongReader("idxptr");
            index_idx_offsets = rb.openULongReader("index_idx_offsets");
        } else {
            throw std::runtime_error(
                std::string("Version does not match ") + versionString(true, 2) + ": " +
                rb.readVersion()
            );
        }

        col_ptr.seek(col_ptr.size() - 1);
        uint64_t count = col_ptr.read_one();
        col_ptr.seek(0);

        NumReader<T> val;
        if constexpr (std::is_same_v<T, uint32_t>) {
            ULongReader val_idx_offsets;
            if (rb.readVersion() == versionString(true, 1)) {
                val_idx_offsets = ConstNumReader<uint64_t>::create({0, UINT64_MAX});
            } else if (rb.readVersion() == versionString(true, 2)) {
                val_idx_offsets = rb.openULongReader("val_idx_offsets");
            }
            val = UIntReader(
                std::make_unique<BP128_FOR_UIntReader>(
                    rb.openUIntReader("val_data"),
                    rb.openUIntReader("val_idx"),
                    std::move(val_idx_offsets),
                    count
                ),
                load_size,
                load_size
            );
        } else {
            val = rb.open<T>("val");
        }

        return StoredMatrix(
            UIntReader(
                std::make_unique<BP128_D1Z_UIntReader>(
                    rb.openUIntReader("index_data"),
                    rb.openUIntReader("index_idx"),
                    std::move(index_idx_offsets),
                    rb.openUIntReader("index_starts"),
                    count
                ),
                load_size,
                load_size
            ),
            std::move(val),
            std::move(col_ptr),
            row_count,
            std::move(row_names),
            std::move(col_names)
        );
    }

    // Open a packed StoredMatrix from a ReaderBuilder, converting row-major orientation to
    // column-major as needed
    static StoredMatrix<T> openPacked(ReaderBuilder &rb, uint32_t load_size = 1024) {
        auto storage_order_reader = rb.openStringReader("storage_order");
        auto storage_order = storage_order_reader->get(0);

        bool row_major = false;
        if (std::string_view("row") == storage_order) row_major = true;
        else if (std::string("col") == storage_order) row_major = false;
        else
            throw std::runtime_error(
                std::string("storage_order must be either \"row\" or \"col\", found: \"") +
                storage_order + "\""
            );

        auto row_names = rb.openStringReader("row_names");
        auto col_names = rb.openStringReader("col_names");

        auto shape = rb.openUIntReader("shape");
        uint32_t row_count = shape.read_one();
        if (row_major) {
            row_count = shape.read_one();
            std::swap(row_names, col_names);
        }

        return StoredMatrix<T>::openPacked(
            rb, load_size, std::move(row_names), std::move(col_names), row_count
        );
    }

    // Return the count of rows and columns
    uint32_t rows() const override { return n_rows; }
    uint32_t cols() const override { return n_cols; }

    const char *rowNames(uint32_t row) override { return row_names->get(row); }
    const char *colNames(uint32_t col) override { return col_names->get(col); }

    // Reset the iterator to start from the beginning
    void restart() override {
        current_col = UINT32_MAX;
        // Don't change current_idx so we will correctly seek when nextCol iscalled
        current_capacity = 0;
        col_ptr.seek(0);
        next_col_ptr = col_ptr.read_one();
    }

    // Seek to a specific column without reading data
    void seekCol(uint32_t col) override {
        current_col = col - 1;
        col_ptr.seek(col);
        next_col_ptr = col_ptr.read_one();
        nextCol();
    }

    // Advance to the next column, or return false if there
    // are no more columns
    bool nextCol() override {
        current_col += 1;
        if (current_col >= n_cols) {
            current_col -= 1;
            return false;
        }

        // Check if we need to perform any seeks, or if we're all good
        // since we just finished a column
        if (next_col_ptr != current_idx) {
            // We need to perform seeks to get to the right data
            // reading location
            current_idx = next_col_ptr;
            val.seek(current_idx);
            row.seek(current_idx);
        }

        next_col_ptr = col_ptr.read_one();
        current_capacity = 0;
        return true;
    }

    // Return the index of the current column
    uint32_t currentCol() const override { return current_col; }

    // Return false if there are no more entries to load
    bool load() override {
        val.advance(current_capacity);
        row.advance(current_capacity);

        if (current_idx >= next_col_ptr) {
            current_capacity = 0;
            return false;
        }

        // Load more data if necessary
        if (val.capacity() == 0) val.ensureCapacity(1);
        if (row.capacity() == 0) row.ensureCapacity(1);

        current_capacity = std::min({val.capacity(), row.capacity(), next_col_ptr - current_idx});
        current_idx += current_capacity;
        return true;
    }

    // Number of loaded entries available
    uint32_t capacity() const override { return current_capacity; }

    // Pointers to the loaded entries
    uint32_t *rowData() override { return row.data(); }
    T *valData() override { return val.data(); }
};

} // end namespace BPCells