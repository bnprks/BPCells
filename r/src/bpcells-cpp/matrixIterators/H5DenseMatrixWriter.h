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
#include "../arrayIO/hdf5.h"
#include "../utils/filesystem_compat.h"
#include "MatrixIterator.h"
#include "OrderRows.h"

namespace BPCells {

template <typename T> class H5DenseMatrixWriter : public MatrixWriter<T> {
  private:
    HighFive::DataSet dataset;
    HighFive::DataType datatype = HighFive::create_datatype<T>();
    uint64_t max_capacity;
    bool row_major;

    static HighFive::DataSet createH5Matrix(
        HighFive::File file,
        std::string dataset_path,
        uint64_t nrow,  // nrow of input matrix
        uint64_t ncol,  // ncol of input matrix
        bool row_major,
        uint64_t chunk_size,
        uint32_t gzip_level
    ) {
      HighFive::SilenceHDF5 s;
      // Create a dataspace with initial shape and max shape
      uint64_t nrow_h5 = nrow;
      uint64_t ncol_h5 = ncol;
      if (row_major) {
          nrow_h5 = ncol;
          ncol_h5 = nrow;
      }
      HighFive::DataSpace dataspace({nrow_h5, ncol_h5}, {HighFive::DataSpace::UNLIMITED, HighFive::DataSpace::UNLIMITED});

      // Use chunking
      HighFive::DataSetCreateProps props;
      props.add(HighFive::Chunking(std::vector<hsize_t>{std::min(chunk_size, nrow_h5), std::min(chunk_size, ncol_h5)}));
      if (gzip_level > 0) {
          props.add(HighFive::Shuffle());
          props.add(HighFive::Deflate(gzip_level));      
      }

      if (file.exist(dataset_path)) {
          file.unlink(dataset_path);
      }
      // Create the dataset
      return file.createDataSet<T>(dataset_path, dataspace, props);
    }

  public:
    H5DenseMatrixWriter(
        const HighFive::File &file,
        std::string path,
        uint64_t nrow,  // nrow of input matrix
        uint64_t ncol,  // ncol of input matrix
        bool row_major,
        uint64_t max_capacity = 16384,
        uint64_t chunk_size = 1024,
        uint32_t gzip_level = 0
    )
        : dataset(createH5Matrix(file, path, nrow, ncol, row_major, chunk_size, gzip_level))
        , max_capacity(max_capacity) 
        , row_major(row_major) {}

    void write(MatrixLoader<T> &mat_in, std::atomic<bool> *user_interrupt = NULL) override {

        // Ensure that we write matrices sorted by row
        OrderRows<T> mat((std::unique_ptr<MatrixLoader<T>>(&mat_in)));
        // Don't delete our original matrix
        mat.preserve_input_loader();

        uint32_t col = 0;
        uint64_t idx = 0; // Index of for col_ptr array
        std::vector<T> zero_buf(mat_in.rows());
        std::vector<T> val_buf = zero_buf;  // buffer for each column
        while (mat.nextCol()) {
            if (user_interrupt != NULL && *user_interrupt) return;
            if (mat.currentCol() < col)
                throw std::runtime_error("H5DenseMatrixWriter encountered out-of-order columns");

            while (mat.load()) {
                uint32_t i = 0;
                while (i < mat.capacity()) {
                    uint32_t capacity = std::min((uint32_t)max_capacity, mat.capacity() - i);
                    for (uint32_t ii = 0; ii < capacity; ii++) {
                        val_buf[*(mat.rowData() + i + ii)] = *(mat.valData() + i + ii);
                    }
                    idx += capacity;
                    i += capacity;
                }
                if (user_interrupt != NULL && *user_interrupt) return;
            }
            if (col == mat.currentCol() && idx > 0) {
                if (row_major) {
                    dataset.select({(uint64_t)mat.currentCol(), 0}, {1, val_buf.size()}).write_raw(val_buf.data(), datatype);
                } else {
                    dataset.select({0, (uint64_t)mat.currentCol()}, {val_buf.size(), 1}).write_raw(val_buf.data(), datatype);
                }
            }
            col += 1;
            val_buf = zero_buf;
        }
        dataset.createAttribute("encoding-type", std::string("array"));
        dataset.createAttribute("encoding-version", std::string("0.2.0"));
    }
};

// Write a Dense Array to an AnnData file
template <typename T>
H5DenseMatrixWriter<T> createAnnDataDenseMatrix(
    std::string file,
    std::string dataset,
    uint32_t nrow,
    uint32_t ncol,
    bool row_major,
    uint32_t buffer_size,
    uint32_t chunk_size,
    uint32_t gzip_level
) {
    HighFive::SilenceHDF5 s;

    std_fs::path path(file);
    if (path.has_parent_path() && !std_fs::exists(path.parent_path())) {
        std_fs::create_directories(path.parent_path());
    }

    HighFive::File h5file(file, HighFive::File::OpenOrCreate);

    return H5DenseMatrixWriter<T>(h5file, dataset, nrow, ncol, row_major, (uint64_t)buffer_size, chunk_size, gzip_level);
}

} // namespace BPCells
