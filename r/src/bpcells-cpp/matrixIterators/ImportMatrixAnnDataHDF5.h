// Copyright 2022 BPCells contributors
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

// Read AnnData sparse matrix, with an implicit transpose to CSC format for
// any data stored in CSR format
template <typename T>
std::unique_ptr<MatrixLoader<T>> openAnnDataMatrix(
    std::string file, std::string group, uint32_t buffer_size, uint32_t read_size = 1024
);

// Read AnnData sparse matrix, with an implicit transpose to CSC format for
// any data stored in CSR format.
// Row/col names are handled as follows:
//   - If row_names or col_names are provided, they are assumed to already have taken into
//     account the internal transpose that will hapen for row-major matrices. i.e. for
//     a row-major X matrix, row_names should be same length as `var/_index`
//   - If row_names or col_names are not provided, they are optimistically inferred by
//     length matching dimensions with the `obs` and `var` index names. This is a compromise
//     to better infer `obsm`, `varm`, `obsp` and `varp` matrix names without having to
//     inspect sub-objects in a way that might break embedded AnnData e.g. with MuData.
//     If `len(obs)` == `len(var)` we just don't infer names for safety.
template <typename T>
std::unique_ptr<MatrixLoader<T>> openAnnDataMatrix(
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