// Copyright 2022 BPCells contributors
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "ImportMatrixHDF5.h"

namespace BPCells {

// Reader interfaces for 10x and AnnData matrices
template <typename T>
StoredMatrix<T> open10xFeatureMatrix(
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
        return StoredMatrix(
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
    return StoredMatrix(
        rb.openULongReader("indices").convert<uint32_t>(),
        rb.open<T>("data"),
        rb.openULongReader("indptr"),
        rows,
        std::move(row_names),
        std::move(col_names)
    );
}

template <typename T>
StoredMatrix<T> open10xFeatureMatrix(
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
        wb.createULongWriter("indices").convert<uint32_t>(),
        wb.create<T>("data"),
        wb.createULongWriter("indptr"),
        wb.createUIntWriter("shape"),
        std::make_unique<NullStringWriter>(),
        std::make_unique<NullStringWriter>(),
        std::make_unique<NullStringWriter>()
    );
}

// Explicit template instantiations for 10x matrix
template StoredMatrix<uint32_t> open10xFeatureMatrix<uint32_t>(
    std::string file,
    std::string group,
    uint32_t buffer_size,
    std::unique_ptr<StringReader> &&row_names,
    std::unique_ptr<StringReader> &&col_names,
    uint32_t read_size
);
template StoredMatrix<uint64_t> open10xFeatureMatrix<uint64_t>(
    std::string file,
    std::string group,
    uint32_t buffer_size,
    std::unique_ptr<StringReader> &&row_names,
    std::unique_ptr<StringReader> &&col_names,
    uint32_t read_size
);
template StoredMatrix<float> open10xFeatureMatrix<float>(
    std::string file,
    std::string group,
    uint32_t buffer_size,
    std::unique_ptr<StringReader> &&row_names,
    std::unique_ptr<StringReader> &&col_names,
    uint32_t read_size
);
template StoredMatrix<double> open10xFeatureMatrix<double>(
    std::string file,
    std::string group,
    uint32_t buffer_size,
    std::unique_ptr<StringReader> &&row_names,
    std::unique_ptr<StringReader> &&col_names,
    uint32_t read_size
);

template StoredMatrix<uint32_t> open10xFeatureMatrix<uint32_t>(
    std::string file, std::string group, uint32_t buffer_size, uint32_t read_size
);
template StoredMatrix<uint64_t> open10xFeatureMatrix<uint64_t>(
    std::string file, std::string group, uint32_t buffer_size, uint32_t read_size
);
template StoredMatrix<float> open10xFeatureMatrix<float>(
    std::string file, std::string group, uint32_t buffer_size, uint32_t read_size
);
template StoredMatrix<double> open10xFeatureMatrix<double>(
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

// Helper class for reading AnnData dense matrix values in transposed order.
// Given a 2D dataset matrix, outputs values row-by-row
template <class T> class H5AnnDataDenseValReader : public BulkNumReader<T> {
  private:
    HighFive::DataSet d_;
    HighFive::DataType type_;
    const uint64_t rows_;
    const uint64_t cols_;
    uint64_t cur_row_ = 0, cur_col_ = 0;

  public:
    H5AnnDataDenseValReader(HighFive::DataSet &&d)
        : d_(std::move(d))
        , type_(d_.getDataType())
        , rows_(d_.getDimensions()[0])
        , cols_(d_.getDimensions()[1]) {}

    // Return total number of integers in the reader
    virtual uint64_t size() const { return rows_ * cols_; }

    // Change the next load to start at index pos
    virtual void seek(uint64_t pos) {
        cur_row_ = pos / cols_;
        cur_col_ = pos % cols_;
    }

    // Copy up to `count` integers into `out`, returning the actual number copied.
    // Will always load >0 if count is >0
    // Note: It is the caller's responsibility to ensure there is no data overflow, i.e.
    // that a load does not try to read past size() total elements
    virtual uint64_t load(T *out, uint64_t count) {
        if (cur_row_ >= rows_) return 0;
        count = std::min(count, cols_ - cur_col_);
        d_.select({cur_row_, cur_col_}, {1, count}).read_raw(out, type_);
        cur_col_ += count;
        if (cur_col_ == cols_) {
            cur_col_ = 0;
            cur_row_ += 1;
        }
        return count;
    }
};

// Helper class for reading AnnData dense matrix values in transposed order.
// Given dimensions of a 2D dataset matrix, output numbers [0,cols) rows times.
class H5AnnDataDenseRowReader : public BulkNumReader<uint32_t> {
  private:
    const uint64_t rows_;
    const uint64_t cols_;
    uint64_t cur_row_ = 0, cur_col_ = 0;

  public:
    H5AnnDataDenseRowReader(uint64_t rows, uint64_t cols) : rows_(rows), cols_(cols) {}
    // Return total number of integers in the reader
    virtual uint64_t size() const { return rows_ * cols_; }

    // Change the next load to start at index pos
    virtual void seek(uint64_t pos) {
        cur_row_ = pos / cols_;
        cur_col_ = pos % cols_;
    }

    // Copy up to `count` integers into `out`, returning the actual number copied.
    // Will always load >0 if count is >0
    // Note: It is the caller's responsibility to ensure there is no data overflow, i.e.
    // that a load does not try to read past size() total elements
    virtual uint64_t load(uint32_t *out, uint64_t count) {
        if (cur_row_ >= rows_) return 0;
        count = std::min(count, cols_ - cur_col_);
        for (uint64_t i = 0; i < count; i++) {
            out[i] = i + cur_col_;
        }
        cur_col_ += count;
        if (cur_col_ == cols_) {
            cur_col_ = 0;
            cur_row_ += 1;
        }
        return count;
    }
};

// Helper class for reading AnnData dense matrix values in transposed order.
// Given dimensions of a 2D dataset matrix, output numbers 0, cols, 2*cols, ..., rows*cols
class H5AnnDataDenseColPtrReader : public BulkNumReader<uint64_t> {
  private:
    const uint64_t rows_;
    const uint64_t cols_;
    uint64_t cur_row_ = 0;

  public:
    H5AnnDataDenseColPtrReader(uint64_t rows, uint64_t cols) : rows_(rows), cols_(cols) {}
    // Return total number of integers in the reader
    virtual uint64_t size() const { return rows_ + 1; }

    // Change the next load to start at index pos
    virtual void seek(uint64_t pos) { cur_row_ = pos; }

    // Copy up to `count` integers into `out`, returning the actual number copied.
    // Will always load >0 if count is >0
    // Note: It is the caller's responsibility to ensure there is no data overflow, i.e.
    // that a load does not try to read past size() total elements
    virtual uint64_t load(uint64_t *out, uint64_t count) {
        if (cur_row_ > rows_) return 0;
        count = std::min(count, rows_ + 1 - cur_row_);
        for (uint64_t i = 0; i < count; i++) {
            out[i] = cols_ * (cur_row_ + i);
        }
        cur_row_ += count;
        return count;
    }
};

void assertAnnDataSparse(std::string file, std::string group) {
    // Check for a dense matrix where we expect a sparse matrix
    HighFive::File f = openH5ForReading(file);
    auto node_type = f.getObjectType(group);
    if (node_type == HighFive::ObjectType::Dataset) {
        throw std::runtime_error(
            "Error in opening AnnData matrix: \"" + group +
            "\" is a dataset rather than a group. This likely indicates a dense matrix which "
            "is not yet supported by BPCells."
        );
    }
}

// Read either obs names or var names from AnnData file
// Args:
//   root: Root group object of file
//   axis: Name of axis to read (either "obs" or "var")
std::unique_ptr<StringReader> readAnnDataDimname(HighFive::Group &root, std::string axis) {
    if (root.getObjectType(axis) == HighFive::ObjectType::Dataset) {
        // Support for legacy format
        std::vector<std::string> name_data;
        readMember(root.getDataSet(axis), "index", name_data);
        return std::make_unique<VecStringReader>(name_data);
    } else {
        std::string dim_ids;
        root.getGroup(axis).getAttribute("_index").read(dim_ids);
        dim_ids = axis + "/" + dim_ids;
        return std::make_unique<H5StringReader>(root, dim_ids);
    }
}

// Read length of either obs or var dimension from AnnData file
// Args:
//   root: Root group object of file
//   axis: Name of axis to read (either "obs" or "var")
size_t readAnnDataDimnameLength(HighFive::Group &root, std::string axis) {
    std::vector<size_t> dims;
    if (root.getObjectType(axis) == HighFive::ObjectType::Dataset) {
        // Support for legacy format
        dims = root.getDataSet(axis).getDimensions();
    } else {
        std::string dim_ids;
        root.getGroup(axis).getAttribute("_index").read(dim_ids);
        dim_ids = axis + "/" + dim_ids;
        dims = root.getDataSet(dim_ids).getDimensions();
    }
    if (dims.size() != 1)
        throw std::runtime_error(
            std::string("readAnnDataDimname(): expected ") + axis +
            " index to be 1-dimensional aray."
        );
    return dims[0];
}

// Read dims and dimnames for an Anndata matrix
// Args:
//  file: File with read access
//  group: string path of matrix to read within file
//  dims: (out) dimensions of matrix
//  row_major: (out) storage order of matrix
//  is_sparse: (out) whether storage format is sparse
// Name matching is done based on dimension length. If obs and var have same length, name
// mis-assignments may occur.
void readAnnDataDims(
    HighFive::File &file,
    std::string group,
    std::vector<uint32_t> &dims,
    bool &row_major,
    bool &is_sparse
) {

    auto node_type = file.getObjectType(group);
    if (node_type == HighFive::ObjectType::Dataset) {
        HighFive::DataSet d = file.getDataSet(group);
        std::vector<size_t> dataset_dims = d.getDimensions();

        dims.resize(dataset_dims.size());
        for (size_t i = 0; i < dims.size(); i++)
            dims[i] = dataset_dims[i];
        row_major = true;
        is_sparse = false;
    } else {
        HighFive::Group g = file.getGroup(group);
        std::string encoding;
        if (g.hasAttribute("h5sparse_format")) {
            // Support for legacy format
            g.getAttribute("h5sparse_format").read(encoding);
            g.getAttribute("h5sparse_shape").read(dims);
            encoding += "_matrix";
        } else if (g.hasAttribute("encoding-type")) {
            g.getAttribute("encoding-type").read(encoding);
            g.getAttribute("shape").read(dims);
        } else {
            throw std::runtime_error(
                "readAnnDataDims(): Could not detect matrix encoding-type on group: " + group
            );
        }

        if (encoding == "csr_matrix") {
            row_major = true;
        } else if (encoding == "csc_matrix") {
            row_major = false;
        } else {
            throw std::runtime_error("readAnnDataDims(): Unsupported matrix encoding: " + encoding);
        }
        is_sparse = true;
    }

    if (dims.size() != 2) {
        throw std::runtime_error(
            "Error in opening AnnData matrix: \"" + group + "\" has " +
            std::to_string(dims.size()) + " dimensions rather than the required two."
        );
    }
}

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
StoredMatrix<T> openAnnDataMatrix(
    std::string file,
    std::string group,
    uint32_t buffer_size,
    std::unique_ptr<StringReader> &&row_names,
    std::unique_ptr<StringReader> &&col_names,
    uint32_t read_size
) {
    std::vector<uint32_t> dims;
    bool row_major;
    bool is_sparse;
    HighFive::SilenceHDF5 s;

    if (group == "") group = "/";

    // Read dims and dimnames
    {
        HighFive::File h5file = openH5ForReading(file);
        readAnnDataDims(h5file, group, dims, row_major, is_sparse);
        HighFive::Group root = h5file.getGroup("/");
        size_t obs_len = readAnnDataDimnameLength(root, "obs");
        size_t var_len = readAnnDataDimnameLength(root, "var");

        if (row_major) std::swap(row_names, col_names);

        // When searching for names, match based on length to obs and var dimensions if required
        if (row_names.get() == nullptr && obs_len != var_len) {
            if (obs_len == dims[0]) row_names = readAnnDataDimname(root, "obs");
            else if (var_len == dims[0]) row_names = readAnnDataDimname(root, "var");
            else row_names = std::make_unique<VecStringReader>(std::vector<std::string>());
        }

        if (col_names.get() == nullptr && obs_len != var_len) {
            if (obs_len == dims[1]) col_names = readAnnDataDimname(root, "obs");
            else if (var_len == dims[1]) col_names = readAnnDataDimname(root, "var");
            else col_names = std::make_unique<VecStringReader>(std::vector<std::string>());
        }

        if (row_major) std::swap(row_names, col_names);
    }

    if (is_sparse) {
        H5ReaderBuilder rb(file, group, 1024, 1024);
        return StoredMatrix<T>(
            rb.openUIntReader("indices"),
            rb.open<T>("data"),
            rb.openULongReader("indptr"),
            row_major ? dims[1] : dims[0],
            std::move(row_names),
            std::move(col_names)
        );
    } else {
        HighFive::File h5file = openH5ForReading(file);
        return StoredMatrix<T>(
            NumReader<uint32_t>(std::make_unique<H5AnnDataDenseRowReader>(dims[0], dims[1]), 1024, 1024),
            NumReader<T>(
                std::make_unique<H5AnnDataDenseValReader<T>>(h5file.getDataSet(group)), 1024, 1024
            ),
            NumReader<uint64_t>(std::make_unique<H5AnnDataDenseColPtrReader>(dims[0], dims[1]), 1024, 1024),
            row_major ? dims[1] : dims[0],
            std::move(row_names),
            std::move(col_names)
        );
    }
}

template <typename T>
StoredMatrix<T>
openAnnDataMatrix(std::string file, std::string group, uint32_t buffer_size, uint32_t read_size) {
    return openAnnDataMatrix<T>(
        file,
        group,
        buffer_size,
        std::unique_ptr<StringReader>(nullptr),
        std::unique_ptr<StringReader>(nullptr),
        read_size
    );
}

// Convenience class for writing shape attribute
template <class T> class H5AttributeNumWriter : public BulkNumWriter<T> {
  private:
    HighFive::Group g;
    std::string attribute_name;
    std::vector<T> data;

  public:
    H5AttributeNumWriter(HighFive::Group g, std::string attribute_name)
        : g(g)
        , attribute_name(attribute_name) {}

    uint64_t write(T *in, uint64_t count) override {
        data.insert(data.end(), in, in + count);
        return count;
    }

    void finalize() override { g.createAttribute(attribute_name, data); }
};

template <typename T>
StoredMatrixWriter<T> createAnnDataMatrix(
    std::string file,
    std::string group,
    bool row_major,
    uint32_t buffer_size,
    uint32_t chunk_size,
    uint32_t gzip_level
) {
    H5WriterBuilder wb(file, group, buffer_size, chunk_size, false, gzip_level);
    HighFive::Group &g = wb.getGroup();

    // See format docs:
    // https://anndata.readthedocs.io/en/latest/fileformat-prose.html#sparse-array-specification-v0-1-0
    g.createAttribute("encoding-type", std::string(row_major ? "csr_matrix" : "csc_matrix"));
    g.createAttribute("encoding-version", std::string("0.1.0"));

    return StoredMatrixWriter<T>(
        wb.createUIntWriter("indices"),
        wb.create<T>("data"),
        wb.createLongWriter("indptr").convert<uint64_t>(),
        UIntWriter(std::make_unique<H5AttributeNumWriter<uint32_t>>(g, "shape"), 16),
        std::make_unique<NullStringWriter>(),
        std::make_unique<NullStringWriter>(),
        std::make_unique<NullStringWriter>(),
        row_major
    );
}

// Add obs + var metadata to an anndata file if needed.
// Draw from the matrix row/col names if present, or otherwise
// insert dummy identifiers.
template <typename T>
void createAnnDataObsVarIfMissing(
    MatrixLoader<T> &mat, std::string file, bool row_major, uint32_t gzip_level
) {
    HighFive::SilenceHDF5 s;
    H5WriterBuilder wb(file, "/", 8192, 1024, true, gzip_level);
    HighFive::Group &g = wb.getGroup();

    const std::string index_name = "bpcells_name";

    std::string row_key = "obs";
    std::string col_key = "var";
    if (row_major) std::swap(row_key, col_key);

    if (!g.exist(row_key)) {
        HighFive::Group metadata = g.createGroup(row_key);
        metadata.createAttribute("encoding-type", std::string("dataframe"));
        metadata.createAttribute("encoding-version", std::string("0.2.0"));
        metadata.createAttribute("_index", index_name);
        metadata.createAttribute("column-order", std::vector<std::string>{index_name});

        std::vector<std::string> row_names(mat.rows());
        if (mat.rows() > 0 && mat.rowNames(0) != NULL) {
            for (uint32_t i = 0; i < mat.rows(); i++) {
                row_names[i] = mat.rowNames(i);
            }
        } else {
            for (uint32_t i = 0; i < mat.rows(); i++) {
                row_names[i] = std::to_string(i);
            }
        }

        wb.createStringWriter(row_key + "/" + index_name)->write(VecStringReader(row_names));
    }
    if (!g.exist(col_key)) {
        HighFive::Group metadata = g.createGroup(col_key);
        metadata.createAttribute("encoding-type", std::string("dataframe"));
        metadata.createAttribute("encoding-version", std::string("0.2.0"));
        metadata.createAttribute("_index", index_name);
        metadata.createAttribute("column-order", std::vector<std::string>{index_name});

        std::vector<std::string> col_names(mat.cols());
        if (mat.cols() > 0 && mat.colNames(0) != NULL) {
            for (uint32_t i = 0; i < mat.cols(); i++) {
                col_names[i] = mat.colNames(i);
            }
        } else {
            for (uint32_t i = 0; i < mat.cols(); i++) {
                col_names[i] = std::to_string(i);
            }
        }

        wb.createStringWriter(col_key + "/" + index_name)->write(VecStringReader(col_names));
    }
}

// Explicit template instantiations
template StoredMatrix<uint32_t> openAnnDataMatrix<uint32_t>(
    std::string file,
    std::string group,
    uint32_t buffer_size,
    std::unique_ptr<StringReader> &&row_names,
    std::unique_ptr<StringReader> &&col_names,
    uint32_t read_size
);
template StoredMatrix<uint64_t> openAnnDataMatrix<uint64_t>(
    std::string file,
    std::string group,
    uint32_t buffer_size,
    std::unique_ptr<StringReader> &&row_names,
    std::unique_ptr<StringReader> &&col_names,
    uint32_t read_size
);
template StoredMatrix<float> openAnnDataMatrix<float>(
    std::string file,
    std::string group,
    uint32_t buffer_size,
    std::unique_ptr<StringReader> &&row_names,
    std::unique_ptr<StringReader> &&col_names,
    uint32_t read_size
);
template StoredMatrix<double> openAnnDataMatrix<double>(
    std::string file,
    std::string group,
    uint32_t buffer_size,
    std::unique_ptr<StringReader> &&row_names,
    std::unique_ptr<StringReader> &&col_names,
    uint32_t read_size
);

template StoredMatrix<uint32_t> openAnnDataMatrix<uint32_t>(
    std::string file, std::string group, uint32_t buffer_size, uint32_t read_size
);
template StoredMatrix<uint64_t> openAnnDataMatrix<uint64_t>(
    std::string file, std::string group, uint32_t buffer_size, uint32_t read_size
);
template StoredMatrix<float> openAnnDataMatrix<float>(
    std::string file, std::string group, uint32_t buffer_size, uint32_t read_size
);
template StoredMatrix<double> openAnnDataMatrix<double>(
    std::string file, std::string group, uint32_t buffer_size, uint32_t read_size
);

template StoredMatrixWriter<uint32_t> createAnnDataMatrix<uint32_t>(
    std::string file,
    std::string group,
    bool row_major,
    uint32_t buffer_size,
    uint32_t chunk_size,
    uint32_t gzip_level
);
template StoredMatrixWriter<uint64_t> createAnnDataMatrix<uint64_t>(
    std::string file,
    std::string group,
    bool row_major,
    uint32_t buffer_size,
    uint32_t chunk_size,
    uint32_t gzip_level
);
template StoredMatrixWriter<float> createAnnDataMatrix<float>(
    std::string file,
    std::string group,
    bool row_major,
    uint32_t buffer_size,
    uint32_t chunk_size,
    uint32_t gzip_level
);
template StoredMatrixWriter<double> createAnnDataMatrix<double>(
    std::string file,
    std::string group,
    bool row_major,
    uint32_t buffer_size,
    uint32_t chunk_size,
    uint32_t gzip_level
);

template void createAnnDataObsVarIfMissing<uint32_t>(
    MatrixLoader<uint32_t> &mat, std::string file, bool row_major, uint32_t gzip_level
);
template void createAnnDataObsVarIfMissing<uint64_t>(
    MatrixLoader<uint64_t> &mat, std::string file, bool row_major, uint32_t gzip_level
);
template void createAnnDataObsVarIfMissing<float>(
    MatrixLoader<float> &mat, std::string file, bool row_major, uint32_t gzip_level
);
template void createAnnDataObsVarIfMissing<double>(
    MatrixLoader<double> &mat, std::string file, bool row_major, uint32_t gzip_level
);
// End explicit template instantiations

std::string getAnnDataMatrixType(std::string file, std::string group) {
    HighFive::SilenceHDF5 s;
    HighFive::File h5file = openH5ForReading(file);
    if (group == "") group = "/";

    HighFive::DataType t;
    if (h5file.getObjectType(group) == HighFive::ObjectType::Group) {
        t = h5file.getGroup(group).getDataSet("data").getDataType();
    } else {
        t = h5file.getDataSet(group).getDataType();
    }

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
        "getAnnDataMatrixType: unrecognized type for group " + group + " in file " + file
    );
}

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

bool isRowOrientedAnnDataMatrix(std::string file, std::string group) {
    HighFive::SilenceHDF5 s;
    HighFive::File h5file = openH5ForReading(file);

    std::vector<uint32_t> dims;
    bool row_major, is_sparse;
    readAnnDataDims(h5file, group, dims, row_major, is_sparse);

    return row_major;
}

template <typename T>
void readMember(HighFive::DataSet &&dataset, std::string name, std::vector<T> &out) {
    auto dims = dataset.getDimensions();
    if (dims.size() != 1) {
        std::ostringstream ss;
        ss << "readMember: dims must be 1, found " << dims.size();
        throw std::runtime_error(ss.str());
    }
    HighFive::CompoundType base_type = dataset.getDataType();
    HighFive::DataType field_type = HighFive::create_datatype<T>();
    bool found = false;
    // Error checking that the type exists
    for (auto &t : base_type.getMembers()) {
        if (t.name == name) {
            if (t.base_type.getClass() == field_type.getClass()) {
                found = true;
                break;
            }
            std::ostringstream ss;
            ss << "Type of member \"" << name << "\" in file (" << t.base_type.string()
               << ") does not match class of requested type (" << field_type.string() << ")";
            throw HighFive::DataTypeException(ss.str());
        }
    }
    if (!found) {
        std::ostringstream ss;
        ss << "Member \"" << name << "\" not found in compound data type";
        throw HighFive::DataTypeException(ss.str());
    }
    // Since HDF5 is pretty good about compatibility, we can ask for a "type" that has just the one
    // field we care about
    HighFive::CompoundType subtype({{name.c_str(), field_type}});

    out.resize(dims[0]);
    // Use a data-converter like the fancy read function (makes string conversion etc. work)
    auto r = HighFive::details::data_converter::get_reader<std::vector<T>>(dims, out, field_type);
    dataset.read_raw(r.getPointer(), subtype);
    // re-arrange results
    r.unserialize(out);
}

} // end namespace BPCells
