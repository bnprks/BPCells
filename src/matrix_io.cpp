#include <Rcpp.h>

#include "matrixIterators/ImportMatrixHDF5.h"
#include "matrixIterators/MatrixIterator.h"
#include "matrixIterators/StoredMatrix.h"

#include "arrayIO/vector.h"
#include "arrayIO/binaryfile.h"
#include "arrayIO/hdf5.h"

#include "R_array_io.h"

using namespace Rcpp;
using namespace BPCells;

// [[Rcpp::export]]
SEXP iterate_packed_matrix_cpp(S4 s4, StringVector row_names, StringVector col_names) {
    S4ReaderBuilder rb(s4);
    return Rcpp::wrap(XPtr<MatrixLoader<uint32_t>>(new StoredMatrix<uint32_t>(
            StoredMatrix<uint32_t>::openPacked(rb, 1024, 
                std::make_unique<RcppStringReader>(row_names), 
                std::make_unique<RcppStringReader>(col_names))
    )));
}

// [[Rcpp::export]]
SEXP write_packed_matrix_cpp(SEXP matrix) {
    XPtr<MatrixLoader<uint32_t>> loader(matrix);
    MatrixIterator iter(*loader);
    
    ListWriterBuilder wb;
    StoredMatrixWriter::createPacked(wb).write(iter, &Rcpp::checkUserInterrupt);

    return wb.getList();
}

// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_cpp(S4 s4, StringVector row_names, StringVector col_names) {
    S4ReaderBuilder rb(s4);
    return Rcpp::wrap(XPtr<MatrixLoader<uint32_t>>(new StoredMatrix<uint32_t>(
            StoredMatrix<uint32_t>::openUnpacked(rb, 
                std::make_unique<RcppStringReader>(row_names), 
                std::make_unique<RcppStringReader>(col_names))
    )));
}

// [[Rcpp::export]]
SEXP write_unpacked_matrix_cpp(SEXP matrix) {
    XPtr<MatrixLoader<uint32_t>> loader(matrix);
    MatrixIterator iter(*loader);
    
    ListWriterBuilder wb;
    StoredMatrixWriter::createUnpacked(wb).write(iter, &Rcpp::checkUserInterrupt);

    return wb.getList();
}

template<class T>
List dims_matrix(const StoredMatrix<T> &mat) {
    IntegerVector dims(2);

    dims[0] = mat.rows();
    dims[1] = mat.cols();
    
    uint32_t row_name_count = mat.rows() == 0 || mat.rowNames(0) == NULL ? 0 : mat.rows();
    uint32_t col_name_count = mat.cols() == 0 || mat.colNames(0) == NULL ? 0 : mat.cols();
    StringVector row_names(row_name_count);
    StringVector col_names(col_name_count);
    for (uint32_t i = 0; i < row_name_count; i++) {
        if (mat.rowNames(i) == NULL) {
            throw std::runtime_error("Matrix has some rownames, but not equal to row number");
        }
        row_names[i] = mat.rowNames(i);
    }
    for (uint32_t i = 0; i < col_name_count; i++) {
        if (mat.colNames(i) == NULL) throw std::runtime_error("Matrix has some colnames, but not equal to col number");
        col_names[i] = mat.colNames(i);
    }

    return List::create(
        Named("dims") = dims,
        Named("row_names") = row_names,
        Named("col_names") = col_names
    );
}

List dims_matrix_reader_builder(ReaderBuilder &rb) {
    std::string version = rb.readVersion();
    
    if (version == "unpacked-uint-matrix-v1") {
        // StoredMatrix<uint32_t> mat = StoredMatrix<uint32_t>::openUnpacked(rb);
        List l = dims_matrix(StoredMatrix<uint32_t>::openUnpacked(rb));
        l["compressed"] = false;
        return l;
    } else if (version == "packed-uint-matrix-v1") {
        // mat = std::make_unique<StoredMatrix<uint32_t>>(StoredMatrix<uint32_t>::openPacked(rb, 1024));
        List l = dims_matrix(StoredMatrix<uint32_t>::openUnpacked(rb));
        l["compressed"] = true;
        return l;
    } else {
        throw std::runtime_error(std::string("Matrix directory has unrecognized version ") + version);
    }
}

// [[Rcpp::export]]
List dims_matrix_file_cpp(std::string dir, uint32_t buffer_size) {
    FileReaderBuilder rb(dir, buffer_size);
    return dims_matrix_reader_builder(rb);
}


// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_file_cpp(std::string dir, uint32_t buffer_size, StringVector row_names, StringVector col_names) {
    FileReaderBuilder rb(dir, buffer_size);
    return Rcpp::wrap(XPtr<MatrixLoader<uint32_t>>(new StoredMatrix<uint32_t>(
            StoredMatrix<uint32_t>::openUnpacked(rb, 
                std::make_unique<RcppStringReader>(row_names), 
                std::make_unique<RcppStringReader>(col_names))
    )));
}

// [[Rcpp::export]]
void write_unpacked_matrix_file_cpp(SEXP matrix, std::string dir, uint32_t buffer_size) {
    XPtr<MatrixLoader<uint32_t>> loader(matrix);
    MatrixIterator<uint32_t> iter(*loader);
    FileWriterBuilder wb(dir, buffer_size);    
    StoredMatrixWriter::createUnpacked(wb).write(iter, &Rcpp::checkUserInterrupt);
}

// [[Rcpp::export]]
SEXP iterate_packed_matrix_file_cpp(std::string dir, uint32_t buffer_size, StringVector row_names, StringVector col_names) {
    FileReaderBuilder rb(dir, buffer_size);
    return Rcpp::wrap(XPtr<MatrixLoader<uint32_t>>(new StoredMatrix<uint32_t>(
            StoredMatrix<uint32_t>::openPacked(rb, 1024, 
                std::make_unique<RcppStringReader>(row_names), 
                std::make_unique<RcppStringReader>(col_names))
    )));
}

// [[Rcpp::export]]
void write_packed_matrix_file_cpp(SEXP matrix, std::string dir, uint32_t buffer_size) {
    XPtr<MatrixLoader<uint32_t>> loader(matrix);
    MatrixIterator<uint32_t> iter(*loader);
    FileWriterBuilder wb(dir, buffer_size);    
    StoredMatrixWriter::createPacked(wb).write(iter, &Rcpp::checkUserInterrupt);
}


// [[Rcpp::export]]
List dims_matrix_hdf5_cpp(std::string file, std::string group, uint32_t buffer_size) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return dims_matrix_reader_builder(rb);
}


// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_hdf5_cpp(std::string file, std::string group, uint32_t buffer_size, StringVector row_names, StringVector col_names) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return Rcpp::wrap(XPtr<MatrixLoader<uint32_t>>(new StoredMatrix<uint32_t>(
            StoredMatrix<uint32_t>::openUnpacked(rb, 
                std::make_unique<RcppStringReader>(row_names), 
                std::make_unique<RcppStringReader>(col_names))
    )));
}

// [[Rcpp::export]]
void write_unpacked_matrix_hdf5_cpp(SEXP matrix, std::string file, std::string group, uint32_t buffer_size, uint32_t chunk_size) {
    XPtr<MatrixLoader<uint32_t>> loader(matrix);
    MatrixIterator<uint32_t> iter(*loader);

    H5WriterBuilder wb(file, group, buffer_size, chunk_size);
    StoredMatrixWriter::createUnpacked(wb).write(iter, &Rcpp::checkUserInterrupt);
}

// [[Rcpp::export]]
SEXP iterate_packed_matrix_hdf5_cpp(std::string file, std::string group, uint32_t buffer_size, StringVector row_names, StringVector col_names) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return Rcpp::wrap(XPtr<MatrixLoader<uint32_t>>(new StoredMatrix<uint32_t>(
            StoredMatrix<uint32_t>::openPacked(rb, 1024, 
                std::make_unique<RcppStringReader>(row_names), 
                std::make_unique<RcppStringReader>(col_names))
    )));
}

// [[Rcpp::export]]
List dims_matrix_10x_hdf5_cpp(std::string file, std::string group, uint32_t buffer_size) {
    StoredMatrix<uint32_t> mat = open10xFeatureMatrix(file, group, buffer_size);
    return dims_matrix(mat);
}

// [[Rcpp::export]]
SEXP iterate_matrix_10x_hdf5_cpp(std::string file, std::string group, uint32_t buffer_size) {
    return Rcpp::wrap(
        XPtr<StoredMatrix<uint32_t>>(new StoredMatrix(open10xFeatureMatrix(file, group, buffer_size)))
    );
}

// [[Rcpp::export]]
List dims_matrix_anndata_hdf5_cpp(std::string file, std::string group, uint32_t buffer_size) {
    List l = dims_matrix(openAnnDataMatrix(file, group, buffer_size));
    l["transpose"] = isRowOrientedAnnDataMatrix(file, group);
    return l;
}

// [[Rcpp::export]]
SEXP iterate_matrix_anndata_hdf5_cpp(std::string file, std::string group, uint32_t buffer_size) {
    return Rcpp::wrap(
        XPtr<StoredMatrix<double>>(new StoredMatrix(openAnnDataMatrix(file, group, buffer_size)))
    );
}