// Copyright 2024 BPCells contributors
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "../utils/filesystem_compat.h"
#include "../arrayIO/hdf5_threadsafe.h"

#include "H5DenseMatrixWriter.h"

namespace BPCells {

template <typename T> class H5DenseMatrixWriter : public MatrixWriter<T> {
  private:
    std::string file_path;
    std::string dataset_path;
    uint64_t chunk_size;
    uint32_t gzip_level;

    bool row_major;

    H5DataSet2D<T> createH5Matrix(uint64_t nrow,  uint64_t ncol) {
        H5Group group {H5Group::create(file_path, dataset_path)}
      // Create a dataspace with initial shape and max shape
      uint64_t nrow_h5 = nrow;
      uint64_t ncol_h5 = ncol;
      if (row_major) {
          return 
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
        const std::string &file_path,
        const std::string &dataset,
        bool row_major,
        uint64_t chunk_size = 1024,
        uint32_t gzip_level = 0
    )
        : file_path(file_path)
        , dataset_path(dataset)
        , chunk_size(chunk_size)
        , gzip_level(gzip_level)
        , row_major(row_major) {}

    void write(MatrixLoader<T> &mat_in, std::atomic<bool> *user_interrupt = NULL) override {
        H5Group group(file_path, "/", BPCells::H5Group::OpenMode::WriteOrCreate);
        
        size_t nrow = mat_in.rows();
        size_t ncol = mat_in.cols();
        if (row_major) std::swap(nrow, ncol);
        

        H5DataSet2D<T> dataset = group.createDataSet2D(dataset_path, nrow, ncol, chunk_size, gzip_level);

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
                    dataset.write_row(mat_in.currentCol(), val_buf);
                } else {
                    dataset.write_col(mat_in.currentCol(), val_buf);
                }
            }
            for (auto &x : val_buf) {
                x = 0;
            }
            loaded = false;
        }
        
        dataset.setAttribute("encoding-type", "array");
        dataset.setAttribute("encoding-version", "0.2.0");
    }
};


template <typename T>
std::unique_ptr<MatrixWriter<T>> createAnnDataDenseMatrix(
    std::string file,
    std::string dataset,
    bool row_major,
    uint32_t chunk_size,
    uint32_t gzip_level
) {
    // Create the HDF5 file
    std_fs::path path(file);
    if (path.has_parent_path() && !std_fs::exists(path.parent_path())) {
      std_fs::create_directories(path.parent_path());
    }
    
    return std::make_unique<H5DenseMatrixWriter<T>>(file, dataset, row_major, chunk_size, gzip_level);
}

// Explicit template instantiations
template std::unique_ptr<MatrixWriter<uint32_t>> createAnnDataDenseMatrix<uint32_t>(
    std::string file,
    std::string dataset,
    bool row_major,
    uint32_t chunk_size,
    uint32_t gzip_level
);
template std::unique_ptr<MatrixWriter<uint64_t>> createAnnDataDenseMatrix<uint64_t>(
    std::string file,
    std::string dataset,
    bool row_major,
    uint32_t chunk_size,
    uint32_t gzip_level
);
template std::unique_ptr<MatrixWriter<float>> createAnnDataDenseMatrix<float>(
    std::string file,
    std::string dataset,
    bool row_major,
    uint32_t chunk_size,
    uint32_t gzip_level
);
template std::unique_ptr<MatrixWriter<double>> createAnnDataDenseMatrix<double>(
    std::string file,
    std::string dataset,
    bool row_major,
    uint32_t chunk_size,
    uint32_t gzip_level
);

} // namespace BPCells
