// Copyright 2024 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once
#include <cstdint>
#include <memory>
#include <string>

#include "MatrixIterator.h"

namespace BPCells {

/**
 * @brief Write a matrix to an AnnData HDF5 file in dense format
 * 
 * @param file Path of the HDF5 file on disk
 * @param dataset Path of the output dataset within the HDF5 file
 * @param row_major If true, write output in transposed (row major) order. (Input is always interpeted as col major)
 * @param chunk_size Chunk size used for HDF5 array storage
 * @param gzip_level If non-zero, level of gzip compression to use while writing.
 */
template <typename T> std::unique_ptr<MatrixWriter<T>> createAnnDataDenseMatrix(
    std::string file,
    std::string dataset,
    bool row_major,
    uint32_t chunk_size,
    uint32_t gzip_level
);

} // namespace BPCells
