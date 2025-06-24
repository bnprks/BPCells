// Copyright 2024 BPCells contributors
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "ImportMatrix10xHDF5.h"
#include "../arrayIO/array_interfaces.h"
#include "../arrayIO/vector.h"

namespace BPCells {

// Reader interfaces for 10x and AnnData matrices
template <typename T>
std::unique_ptr<MatrixLoader<T>> open10xFeatureMatrix(
    std::string file,
    std::string group,
    uint32_t buffer_size,
    std::unique_ptr<StringReader> &&row_names,
    std::unique_ptr<StringReader> &&col_names,
    uint32_t read_size
) {
    H5ReaderBuilder rb(file, group, buffer_size, read_size);
    if (rb.getGroup().exist("features/id")) {
        // Most up-to-date matrix format (cellranger V3+)
        uint32_t rows = rb.openUIntReader("shape").read_one();
        if (row_names.get() == nullptr) {
            row_names = rb.openStringReader("features/id");
        }
        if (col_names.get() == nullptr) {
            col_names = rb.openStringReader("barcodes");
        }
        return std::make_unique<StoredMatrix<T>>(
            rb.openULongReader("indices").convert<uint32_t>(),
            rb.open<T>("data"),
            rb.openULongReader("indptr"),
            rows,
            std::move(row_names),
            std::move(col_names)
        );
    }

    // Older-style 10x matrix format
    uint32_t rows = rb.openUIntReader("shape").read_one();
    if (row_names.get() == nullptr) {
        row_names = rb.openStringReader("genes");
    }
    if (col_names.get() == nullptr) {
        col_names = rb.openStringReader("barcodes");
    }
    return std::make_unique<StoredMatrix<T>>(
        rb.openULongReader("indices").convert<uint32_t>(),
        rb.open<T>("data"),
        rb.openULongReader("indptr"),
        rows,
        std::move(row_names),
        std::move(col_names)
    );
}

template <typename T>
std::unique_ptr<MatrixLoader<T>> open10xFeatureMatrix(
    std::string file, std::string group, uint32_t buffer_size, uint32_t read_size
) {
    return open10xFeatureMatrix<T>(
        file,
        group,
        buffer_size,
        std::unique_ptr<StringReader>(nullptr),
        std::unique_ptr<StringReader>(nullptr),
        read_size
    );
}

template <typename T>
StoredMatrixWriter<T> create10xFeatureMatrix(
    std::string file_path,
    StringReader &&barcodes,
    StringReader &&feature_ids,
    StringReader &&feature_names,
    StringReader &&feature_types,
    const std::map<std::string, std::unique_ptr<StringReader>> &feature_metadata,
    uint32_t buffer_size,
    uint32_t chunk_size,
    uint32_t gzip_level
) {
    H5WriterBuilder wb(file_path, "matrix", buffer_size, chunk_size, false, gzip_level);

    wb.createStringWriter("barcodes")->write(barcodes);
    wb.createStringWriter("features/id")->write(feature_ids);
    wb.createStringWriter("features/name")->write(feature_names);
    wb.createStringWriter("features/feature_type")->write(feature_types);
    std::vector<std::string> tag_keys;
    for (const auto &[key, value] : feature_metadata) {
        wb.createStringWriter(std::string("features/") + key)->write(*value);
        tag_keys.push_back(key);
    }
    wb.createStringWriter("features/_all_tag_keys")->write(VecStringReader(tag_keys));

    return StoredMatrixWriter(
        wb.create<int64_t>("indices").convert<uint32_t>(),
        wb.create<T>("data"),
        wb.create<int64_t>("indptr").convert<uint64_t>(),
        wb.create<int32_t>("shape").convert<uint32_t>(),
        std::make_unique<NullStringWriter>(),
        std::make_unique<NullStringWriter>(),
        std::make_unique<NullStringWriter>()
    );
}

// Explicit template instantiations for 10x matrix
template std::unique_ptr<MatrixLoader<uint32_t>> open10xFeatureMatrix<uint32_t>(
    std::string file,
    std::string group,
    uint32_t buffer_size,
    std::unique_ptr<StringReader> &&row_names,
    std::unique_ptr<StringReader> &&col_names,
    uint32_t read_size
);
template std::unique_ptr<MatrixLoader<uint64_t>> open10xFeatureMatrix<uint64_t>(
    std::string file,
    std::string group,
    uint32_t buffer_size,
    std::unique_ptr<StringReader> &&row_names,
    std::unique_ptr<StringReader> &&col_names,
    uint32_t read_size
);
template std::unique_ptr<MatrixLoader<float>> open10xFeatureMatrix<float>(
    std::string file,
    std::string group,
    uint32_t buffer_size,
    std::unique_ptr<StringReader> &&row_names,
    std::unique_ptr<StringReader> &&col_names,
    uint32_t read_size
);
template std::unique_ptr<MatrixLoader<double>> open10xFeatureMatrix<double>(
    std::string file,
    std::string group,
    uint32_t buffer_size,
    std::unique_ptr<StringReader> &&row_names,
    std::unique_ptr<StringReader> &&col_names,
    uint32_t read_size
);

template std::unique_ptr<MatrixLoader<uint32_t>> open10xFeatureMatrix<uint32_t>(
    std::string file, std::string group, uint32_t buffer_size, uint32_t read_size
);
template std::unique_ptr<MatrixLoader<uint64_t>> open10xFeatureMatrix<uint64_t>(
    std::string file, std::string group, uint32_t buffer_size, uint32_t read_size
);
template std::unique_ptr<MatrixLoader<float>> open10xFeatureMatrix<float>(
    std::string file, std::string group, uint32_t buffer_size, uint32_t read_size
);
template std::unique_ptr<MatrixLoader<double>> open10xFeatureMatrix<double>(
    std::string file, std::string group, uint32_t buffer_size, uint32_t read_size
);

template StoredMatrixWriter<uint32_t> create10xFeatureMatrix<uint32_t>(
    std::string file_path,
    StringReader &&barcodes,
    StringReader &&feature_ids,
    StringReader &&feature_names,
    StringReader &&feature_types,
    const std::map<std::string, std::unique_ptr<StringReader>> &feature_metadata,
    uint32_t buffer_size,
    uint32_t chunk_size,
    uint32_t gzip_level
);
template StoredMatrixWriter<uint64_t> create10xFeatureMatrix<uint64_t>(
    std::string file_path,
    StringReader &&barcodes,
    StringReader &&feature_ids,
    StringReader &&feature_names,
    StringReader &&feature_types,
    const std::map<std::string, std::unique_ptr<StringReader>> &feature_metadata,
    uint32_t buffer_size,
    uint32_t chunk_size,
    uint32_t gzip_level
);
template StoredMatrixWriter<float> create10xFeatureMatrix<float>(
    std::string file_path,
    StringReader &&barcodes,
    StringReader &&feature_ids,
    StringReader &&feature_names,
    StringReader &&feature_types,
    const std::map<std::string, std::unique_ptr<StringReader>> &feature_metadata,
    uint32_t buffer_size,
    uint32_t chunk_size,
    uint32_t gzip_level
);
template StoredMatrixWriter<double> create10xFeatureMatrix<double>(
    std::string file_path,
    StringReader &&barcodes,
    StringReader &&feature_ids,
    StringReader &&feature_names,
    StringReader &&feature_types,
    const std::map<std::string, std::unique_ptr<StringReader>> &feature_metadata,
    uint32_t buffer_size,
    uint32_t chunk_size,
    uint32_t gzip_level
);
// End explicit template instantiations for 10x matrix

std::string get10xMatrixType(std::string file, std::string group) {
    H5ReaderBuilder rb(file, group, 1024, 1024);
    HighFive::SilenceHDF5 s;
    HighFive::DataType t = rb.getGroup().getDataSet("data").getDataType();
    if (t == HighFive::AtomicType<int>()) {
        return "uint32_t";
    } else if (t == HighFive::AtomicType<uint32_t>()) {
        return "uint32_t";
    } else if (t == HighFive::AtomicType<uint64_t>()) {
        return "uint64_t";
    } else if (t == HighFive::AtomicType<float>()) {
        return "float";
    } else if (t == HighFive::AtomicType<double>()) {
        return "double";
    }
    throw std::runtime_error(
        "get10xMatrixType: unrecognized type for group " + group + " in file " + file
    );
}

} // end namespace BPCells
