#include <Rcpp.h>

#include "matrixIterators/MatrixIterator.h"
#include "matrixIterators/PackedMatrix.h"
#include "arrayIO/vector.h"
#include "arrayIO/binaryfile.h"
#include "arrayIO/hdf5.h"

using namespace Rcpp;
using namespace BPCells;

// PACKED MATRIX

// [[Rcpp::export]]
List write_packed_matrix_cpp(SEXP matrix) {
    XPtr<MatrixLoader<uint32_t>> loader(matrix);
    
    std::vector<uint32_t> val_data, val_idx, row_data, row_starts, row_idx, col_ptr, row_count;

    auto writer = PackedMatrixWriter(
        std::make_unique<ZVecUIntWriter>(val_data), std::make_unique<ZVecUIntWriter>(val_idx), std::make_unique<ZVecUIntWriter>(row_data), 
        std::make_unique<ZVecUIntWriter>(row_starts), std::make_unique<ZVecUIntWriter>(row_idx), std::make_unique<ZVecUIntWriter>(col_ptr), 
        std::make_unique<ZVecUIntWriter>(row_count)
    );
    writer.write(*loader, &Rcpp::checkUserInterrupt);
    loader.release();
    return List::create(
        Named("val_data") = IntegerVector(val_data.begin(), val_data.end()),
        Named("val_idx") = IntegerVector(val_idx.begin(), val_idx.end()),
        Named("row_data") = IntegerVector(row_data.begin(), row_data.end()),
        Named("row_starts") = IntegerVector(row_starts.begin(), row_starts.end()),
        Named("row_idx") = IntegerVector(row_idx.begin(), row_idx.end()),
        Named("col_ptr") = IntegerVector(col_ptr.begin(), col_ptr.end()),
        Named("row_count") = IntegerVector(row_count.begin(), row_count.end())
    );
}

// [[Rcpp::export]]
SEXP iterate_packed_matrix_cpp(const S4 s4) {
    IntegerVector val_data = s4.slot("val_data");
    IntegerVector val_idx = s4.slot("val_idx");
    IntegerVector row_data = s4.slot("row_data");
    IntegerVector row_starts = s4.slot("row_starts");
    IntegerVector row_idx = s4.slot("row_idx");
    IntegerVector col_ptr = s4.slot("col_ptr");
    IntegerVector row_count = s4.slot("row_count");
    return Rcpp::wrap(
        XPtr<MatrixLoader<uint32_t>>(new PackedMatrix(
            std::make_unique<ZVecUIntReader>((const uint32_t *) val_data.cbegin(), val_data.size()),
            std::make_unique<ZVecUIntReader>((const uint32_t *) val_idx.cbegin(), val_idx.size()),
            std::make_unique<ZVecUIntReader>((const uint32_t *) row_data.cbegin(), row_data.size()),
            std::make_unique<ZVecUIntReader>((const uint32_t *) row_starts.cbegin(), row_starts.size()),
            std::make_unique<ZVecUIntReader>((const uint32_t *) row_idx.cbegin(), row_idx.size()),
            std::make_unique<ZVecUIntReader>((const uint32_t *) col_ptr.cbegin(), col_ptr.size()),
            std::make_unique<ZVecUIntReader>((const uint32_t *) row_count.cbegin(), row_count.size())
        ))
    );
}


// [[Rcpp::export]]
List dims_matrix_file_cpp(std::string dir, uint32_t buffer_size) {
    std::string version(loadVersionMatrixDir(dir));
    IntegerVector dims(2);
    bool compressed = false;
    if (version == "v1-unpacked") {
        UnpackedMatrix mat(openUnpackedMatrixDir(dir, buffer_size));
        dims[0] = mat.rows();
        dims[1] = mat.cols();
    } else if (version == "v1-packed") {
        PackedMatrix mat(openPackedMatrixDir(dir, buffer_size));
        dims[0] = mat.rows();
        dims[1] = mat.cols();
        compressed = true;
    } else {
        throw std::runtime_error(std::string("Packed matrix directory has unrecognized version ") + version);
    }
    
    return List::create(
        Named("dims") = dims,
        Named("compressed") = compressed
    );
}

// [[Rcpp::export]]
SEXP iterate_packed_matrix_file_cpp(std::string dir, uint32_t buffer_size) {
    return Rcpp::wrap(
        XPtr<PackedMatrix>(new PackedMatrix(openPackedMatrixDir(dir, buffer_size)))
    );
}

// [[Rcpp::export]]
void write_packed_matrix_file_cpp(SEXP mat, std::string dir, uint32_t buffer_size) {
    XPtr<MatrixLoader<uint32_t>> loader(mat);
    
    PackedMatrixWriter writer(createPackedMatrixDir(dir, buffer_size));

    writer.write(*loader, &Rcpp::checkUserInterrupt);
}


// [[Rcpp::export]]
List dims_matrix_hdf5_cpp(std::string file, std::string group, uint32_t buffer_size) {
    std::string version(loadVersionMatrixH5(file, group));
    IntegerVector dims(2);
    bool compressed = false;
    if (version == "v1-unpacked") {
        UnpackedMatrix mat(openUnpackedMatrixH5(file, group, buffer_size));
        dims[0] = mat.rows();
        dims[1] = mat.cols();
    } else if (version == "v1-packed") {
        PackedMatrix mat(openPackedMatrixH5(file, group, buffer_size));
        dims[0] = mat.rows();
        dims[1] = mat.cols();
        compressed = true;
    } else {
        throw std::runtime_error(std::string("Packed matrix hdf5 has unrecognized version ") + version);
    }
    
    return List::create(
        Named("dims") = dims,
        Named("compressed") = compressed
    );
}

// [[Rcpp::export]]
SEXP iterate_packed_matrix_hdf5_cpp(std::string file, std::string group, uint32_t buffer_size) {
    return Rcpp::wrap(
        XPtr<PackedMatrix>(new PackedMatrix(openPackedMatrixH5(file, group, buffer_size)))
    );
}

// [[Rcpp::export]]
void write_packed_matrix_hdf5_cpp(SEXP mat, std::string file, std::string group, uint32_t buffer_size, uint32_t chunk_size) {
    XPtr<MatrixLoader<uint32_t>> loader(mat);
    
    PackedMatrixWriter writer(createPackedMatrixH5(file, group, buffer_size, chunk_size));

    writer.write(*loader, &Rcpp::checkUserInterrupt);
}


// UNPACKED MATRIX

// [[Rcpp::export]]
List write_unpacked_matrix_cpp(SEXP matrix) {
    XPtr<MatrixLoader<uint32_t>> loader(matrix);
    
    std::vector<uint32_t> val, row, col_ptr, row_count;

    auto writer = UnpackedMatrixWriter(
        std::make_unique<ZVecUIntWriter>(val), std::make_unique<ZVecUIntWriter>(row), 
        std::make_unique<ZVecUIntWriter>(col_ptr), std::make_unique<ZVecUIntWriter>(row_count)
    );
    writer.write(*loader, &Rcpp::checkUserInterrupt);
    loader.release();
    return List::create(
        Named("val") = IntegerVector(val.begin(), val.end()),
        Named("row") = IntegerVector(row.begin(), row.end()),
        Named("col_ptr") = IntegerVector(col_ptr.begin(), col_ptr.end()),
        Named("row_count") = IntegerVector(row_count.begin(), row_count.end())
    );
}

// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_cpp(const S4 s4) {
    IntegerVector val = s4.slot("val");
    IntegerVector row = s4.slot("row");
    IntegerVector col_ptr = s4.slot("col_ptr");
    IntegerVector row_count = s4.slot("row_count");
    return Rcpp::wrap(
        XPtr<MatrixLoader<uint32_t>>(new UnpackedMatrix(
            std::make_unique<ZVecUIntReader>((const uint32_t *) val.cbegin(), val.size()),
            std::make_unique<ZVecUIntReader>((const uint32_t *) row.cbegin(), row.size()),
            std::make_unique<ZVecUIntReader>((const uint32_t *) col_ptr.cbegin(), col_ptr.size()),
            std::make_unique<ZVecUIntReader>((const uint32_t *) row_count.cbegin(), row_count.size())
        ))
    );
}


// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_file_cpp(std::string dir, uint32_t buffer_size) {
    return Rcpp::wrap(
        XPtr<UnpackedMatrix>(new UnpackedMatrix(openUnpackedMatrixDir(dir, buffer_size)))
    );
}

// [[Rcpp::export]]
void write_unpacked_matrix_file_cpp(SEXP mat, std::string dir, uint32_t buffer_size) {
    XPtr<MatrixLoader<uint32_t>> loader(mat);
    
    UnpackedMatrixWriter writer(createUnpackedMatrixDir(dir, buffer_size));

    writer.write(*loader, &Rcpp::checkUserInterrupt);
}


// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_hdf5_cpp(std::string file, std::string group, uint32_t buffer_size) {
    return Rcpp::wrap(
        XPtr<UnpackedMatrix>(new UnpackedMatrix(openUnpackedMatrixH5(file, group, buffer_size)))
    );
}

// [[Rcpp::export]]
void write_unpacked_matrix_hdf5_cpp(SEXP mat, std::string file, std::string group, uint32_t buffer_size, uint32_t chunk_size) {
    XPtr<MatrixLoader<uint32_t>> loader(mat);
    
    UnpackedMatrixWriter writer(createUnpackedMatrixH5(file, group, buffer_size, chunk_size));

    writer.write(*loader, &Rcpp::checkUserInterrupt);
}

// [[Rcpp::export]]
SEXP iterate_matrix_10x_hdf5_cpp(std::string file, uint32_t buffer_size) {
    return Rcpp::wrap(
        XPtr<UnpackedMatrix>(new UnpackedMatrix(open10xFeatureMatrix(file, "matrix", buffer_size)))
    );
}

// [[Rcpp::export]]
List dims_matrix_10x_hdf5_cpp(std::string file, uint32_t buffer_size) {
    IntegerVector dims(2);
    UnpackedMatrix mat(open10xFeatureMatrix(file, "matrix", buffer_size));
    dims[0] = mat.rows();
    dims[1] = mat.cols();

    return List::create(
        Named("dims") = dims    
    );
}