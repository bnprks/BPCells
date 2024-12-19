// Copyright 2024 BPCells contributors
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "H5DenseMatrixWriter.h"

namespace BPCells {

template <typename T>
H5DenseMatrixWriter<T> createAnnDataDenseMatrix(
    std::string file,
    std::string dataset,
    bool row_major,
    uint32_t buffer_size,
    uint32_t chunk_size,
    uint32_t gzip_level
) {
    HighFive::SilenceHDF5 s;

    // Create the HDF5 file
    std_fs::path path(file);
    if (path.has_parent_path() && !std_fs::exists(path.parent_path())) {
      std_fs::create_directories(path.parent_path());
    }
    HighFive::File h5file(file, HighFive::File::OpenOrCreate);
    
    return H5DenseMatrixWriter<T>(h5file, dataset, row_major, buffer_size, chunk_size, gzip_level);
}

// Explicit template instantiations
template H5DenseMatrixWriter<uint32_t> createAnnDataDenseMatrix<uint32_t>(
    std::string file,
    std::string dataset,
    bool row_major,
    uint32_t buffer_size,
    uint32_t chunk_size,
    uint32_t gzip_level
);
template H5DenseMatrixWriter<uint64_t> createAnnDataDenseMatrix<uint64_t>(
    std::string file,
    std::string dataset,
    bool row_major,
    uint32_t buffer_size,
    uint32_t chunk_size,
    uint32_t gzip_level
);
template H5DenseMatrixWriter<float> createAnnDataDenseMatrix<float>(
    std::string file,
    std::string dataset,
    bool row_major,
    uint32_t buffer_size,
    uint32_t chunk_size,
    uint32_t gzip_level
);
template H5DenseMatrixWriter<double> createAnnDataDenseMatrix<double>(
    std::string file,
    std::string dataset,
    bool row_major,
    uint32_t buffer_size,
    uint32_t chunk_size,
    uint32_t gzip_level
);

} // namespace BPCells
