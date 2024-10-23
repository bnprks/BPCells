// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once
#include "../arrayIO/array_interfaces.h"
#include "../arrayIO/hdf5.h"
#include "../arrayIO/vector.h"
#include "StoredMatrix.h"
#include "StoredMatrixWriter.h"

namespace BPCells {

// Read AnnData sparse matrix, with an implicit transpose to CSC format for
// any data stored in CSR format
template <typename T>
StoredMatrix<T> openAnnDataMatrix(
    std::string file, std::string group, uint32_t buffer_size, uint32_t read_size = 1024
);

template <typename T>
StoredMatrix<T> openAnnDataMatrix(
    std::string file,
    std::string group,
    uint32_t buffer_size,
    std::unique_ptr<StringReader> &&row_names,
    std::unique_ptr<StringReader> &&col_names,
    uint32_t read_size = 1024
);

// Write a Sparse Array to an AnnData file
template <typename T>
StoredMatrixWriter<T> createAnnDataMatrix(
    std::string file,
    std::string group,
    bool row_major,
    uint32_t buffer_size,
    uint32_t chunk_size,
    uint32_t gzip_level
);

// Add obs + var metadata to an anndata file if needed.
// Draw from the matrix row/col names if present, or otherwise
// insert dummy identifiers.
template <typename T>
void createAnnDataObsVarIfMissing(
    MatrixLoader<T> &mat, std::string file, bool row_major, uint32_t gzip_level
);

std::string getAnnDataMatrixType(std::string file, std::string group);

bool isRowOrientedAnnDataMatrix(std::string file, std::string group);

// Read a single member from a struct-type dataset
template <typename T>
void readMember(HighFive::DataSet &&dataset, std::string name, std::vector<T> &out);

} // end namespace BPCells