// Copyright 2024 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once
#include "../arrayIO/hdf5.h"
#include "StoredMatrix.h"
#include "StoredMatrixWriter.h"


namespace BPCells {

// Reader interfaces for 10x and AnnData matrices
// Open 10x matrix should not assume uint32_t
template <typename T>
std::unique_ptr<MatrixLoader<T>> open10xFeatureMatrix(
    std::string file, std::string group, uint32_t buffer_size, uint32_t read_size = 1024
);

template <typename T>
std::unique_ptr<MatrixLoader<T>> open10xFeatureMatrix(
    std::string file,
    std::string group,
    uint32_t buffer_size,
    std::unique_ptr<StringReader> &&row_names,
    std::unique_ptr<StringReader> &&col_names,
    uint32_t read_size = 1024
);

template <typename T>
StoredMatrixWriter<T> create10xFeatureMatrix(
    std::string file,
    StringReader &&barcodes,
    StringReader &&feature_ids,
    StringReader &&feature_names,
    StringReader &&feature_types,
    const std::map<std::string, std::unique_ptr<StringReader>> &feature_metadata,
    uint32_t buffer_size,
    uint32_t chunk_size,
    uint32_t gzip_level
);

std::string get10xMatrixType(std::string file, std::string group);

} // end namespace BPCells