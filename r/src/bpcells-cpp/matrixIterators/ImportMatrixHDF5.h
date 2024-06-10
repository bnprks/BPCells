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

// Helper class for when we need to fake reading a single integer from memory
template <class T> class SingletonNumReader : public BulkNumReader<T> {
  private:
    T num;
    bool read = false;

  public:
    SingletonNumReader(T num) : num(num) {}
    uint64_t size() const override { return 1; }
    void seek(uint64_t pos) override { read = pos > 0; }
    uint64_t load(T *out, uint64_t count) override {
        if (read) return 0;
        out[0] = num;
        return 1;
    }
};

// Reader interfaces for 10x and AnnData matrices
// Open 10x matrix should not assume uint32_t
template <typename T>
StoredMatrix<T> open10xFeatureMatrix(
    std::string file, std::string group, uint32_t buffer_size, uint32_t read_size = 1024
);

template <typename T>
StoredMatrix<T> open10xFeatureMatrix(
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
std::string get10xMatrixType(std::string file, std::string group);

bool isRowOrientedAnnDataMatrix(std::string file, std::string group);

// Read a single member from a struct-type dataset
template <typename T>
void readMember(HighFive::DataSet &&dataset, std::string name, std::vector<T> &out);

} // end namespace BPCells