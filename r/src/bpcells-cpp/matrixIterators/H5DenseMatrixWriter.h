// Copyright 2024 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once
#include "../arrayIO/hdf5.h"
#include "../utils/filesystem_compat.h"
#include "OrderRows.h"

namespace BPCells {

template <typename T> class H5DenseMatrixWriter : public MatrixWriter<T> {
  private:
    HighFive::File h5file;
    std::string dataset_path;
    uint64_t chunk_size;
    uint32_t gzip_level;

    HighFive::DataType datatype = HighFive::create_datatype<T>();
    bool row_major;

    HighFive::DataSet createH5Matrix(uint64_t nrow,  uint64_t ncol) {
      HighFive::SilenceHDF5 s;
      
      // Create a dataspace with initial shape and max shape
      uint64_t nrow_h5 = nrow;
      uint64_t ncol_h5 = ncol;
      if (row_major) {
          nrow_h5 = ncol;
          ncol_h5 = nrow;
      }
      HighFive::DataSpace dataspace({nrow_h5, ncol_h5});

      // Use chunking
      HighFive::DataSetCreateProps props;
      props.add(HighFive::Chunking(std::vector<hsize_t>{std::min(chunk_size, nrow_h5), std::min(chunk_size, ncol_h5)}));
      if (gzip_level > 0) {
          props.add(HighFive::Shuffle());
          props.add(HighFive::Deflate(gzip_level));      
      }
      // Create the dataset
      return h5file.createDataSet<T>(dataset_path, dataspace, props);
    }

  public:
    H5DenseMatrixWriter(
        HighFive::File h5file,
        std::string dataset,
        bool row_major,
        uint64_t chunk_size = 1024,
        uint32_t gzip_level = 0
    )
        : h5file(h5file)
        , dataset_path(dataset)
        , chunk_size(chunk_size)
        , gzip_level(gzip_level)
        , row_major(row_major) {}

    void write(MatrixLoader<T> &mat_in, std::atomic<bool> *user_interrupt = NULL) override {

        HighFive::DataSet h5dataset = createH5Matrix(mat_in.rows(), mat_in.cols());

        bool loaded = false; // Any non-zero values has been loaded.
        std::vector<T> val_buf(mat_in.rows());  // buffer for each column
        while (mat_in.nextCol()) {
            if (user_interrupt != NULL && *user_interrupt) return;
            while (mat_in.load()) {
                loaded = true;
                uint32_t *row_data = mat_in.rowData();
                T *val_data = mat_in.valData();
                uint32_t capacity = mat_in.capacity();
                for (uint32_t i = 0; i < capacity; i++) {
                    val_buf[row_data[i]] = val_data[i];
                }
                if (user_interrupt != NULL && *user_interrupt) return;
            }
            if (loaded) {
                if (row_major) {
                    h5dataset.select({(uint64_t)mat_in.currentCol(), 0}, {1, val_buf.size()}).write_raw(val_buf.data(), datatype);
                } else {
                    h5dataset.select({0, (uint64_t)mat_in.currentCol()}, {val_buf.size(), 1}).write_raw(val_buf.data(), datatype);
                }
            }
            for (auto &x : val_buf) {
                x = 0;
            }
            loaded = false;
        }
        h5dataset.createAttribute("encoding-type", std::string("array"));
        h5dataset.createAttribute("encoding-version", std::string("0.2.0"));
    }
};

// Write a Dense Array to an AnnData file
template <typename T> H5DenseMatrixWriter<T> createAnnDataDenseMatrix(
    std::string file,
    std::string dataset,
    bool row_major,
    uint32_t chunk_size,
    uint32_t gzip_level
);

} // namespace BPCells
