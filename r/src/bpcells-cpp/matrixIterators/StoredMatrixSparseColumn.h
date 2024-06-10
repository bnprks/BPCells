// Copyright 2023 BPCells contributors
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include "../arrayIO/array_interfaces.h"
#include "../arrayIO/bp128.h"
#include "StoredMatrix.h"
#include "StoredMatrixWriter.h"

namespace BPCells {

// This contains experimental code to read/write stored matrices where the idxptr is
// stored with BP128-d1 encoding. This helps lower storage sizes of matrices with
// many empty columns, though this implementation will cause errors if there are
// >=2^32-1 non-zero entries in the matrix

template <typename T>
StoredMatrix<T> EXPERIMENTAL_openPackedSparseColumn(ReaderBuilder &rb, uint32_t load_size = 1024) {
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
    uint32_t col_count = shape.read_one();
    if (row_major) {
        std::swap(row_count, col_count);
        std::swap(row_names, col_names);
    }

    ULongReader col_ptr, index_idx_offsets;
    if (rb.readVersion() == StoredMatrix<T>::versionString(true, 9999)) {
        col_ptr = UIntReader(
                      std::make_unique<BP128_D1_UIntReader>(
                          rb.openUIntReader("idxptr_data"),
                          rb.openUIntReader("idxptr_idx"),
                          rb.openULongReader("idxptr_idx_offsets"),
                          rb.openUIntReader("idxptr_starts"),
                          col_count + 1
                      ),
                      load_size,
                      load_size
        )
                      .convert<uint64_t>();
        index_idx_offsets = rb.openULongReader("index_idx_offsets");
    } else {
        throw std::runtime_error(
            std::string("Version does not match ") + StoredMatrix<T>::versionString(true, 9999) +
            ": " + rb.readVersion()
        );
    }

    col_ptr.seek(col_ptr.size() - 1);
    uint64_t count = col_ptr.read_one();
    col_ptr.seek(0);

    NumReader<T> val;
    if constexpr (std::is_same_v<T, uint32_t>) {
        ULongReader val_idx_offsets;
        val = UIntReader(
            std::make_unique<BP128_FOR_UIntReader>(
                rb.openUIntReader("val_data"),
                rb.openUIntReader("val_idx"),
                rb.openULongReader("val_idx_offsets"),
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
                rb.openULongReader("index_idx_offsets"),
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

template <typename T>
StoredMatrixWriter<T> EXPERIMENTAL_createPackedSparseColumn(
    WriterBuilder &wb, bool row_major = false, uint32_t buffer_size = 1024
) {
    wb.writeVersion(StoredMatrix<T>::versionString(true, 9999));
    NumWriter<T> val;

    if constexpr (std::is_same_v<T, uint32_t>) {
        val = UIntWriter(
            std::make_unique<BP128_FOR_UIntWriter>(
                wb.createUIntWriter("val_data"),
                wb.createUIntWriter("val_idx"),
                wb.createULongWriter("val_idx_offsets")
            ),
            buffer_size
        );
    } else {
        val = wb.create<T>("val");
    }

    return StoredMatrixWriter(
        UIntWriter(
            std::make_unique<BP128_D1Z_UIntWriter>(
                wb.createUIntWriter("index_data"),
                wb.createUIntWriter("index_idx"),
                wb.createULongWriter("index_idx_offsets"),
                wb.createUIntWriter("index_starts")
            ),
            buffer_size
        ),
        std::move(val),
        UIntWriter(
            std::make_unique<BP128_D1_UIntWriter>(
                wb.createUIntWriter("idxptr_data"),
                wb.createUIntWriter("idxptr_idx"),
                wb.createULongWriter("idxptr_idx_offsets"),
                wb.createUIntWriter("idxptr_starts")
            ),
            buffer_size
        )
            .convert<uint64_t>(),
        wb.createUIntWriter("shape"),
        wb.createStringWriter("row_names"),
        wb.createStringWriter("col_names"),
        wb.createStringWriter("storage_order"),
        row_major
    );
}

} // end namespace BPCells