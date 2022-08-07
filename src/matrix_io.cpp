#include <Rcpp.h>

#include "matrixIterators/ImportMatrixHDF5.h"
#include "matrixIterators/MatrixIterator.h"
#include "matrixIterators/StoredMatrix.h"
#include "matrixIterators/StoredMatrixWriter.h"
#include "matrixIterators/StoredMatrixTransposeWriter.h"

#include "arrayIO/vector.h"
#include "arrayIO/binaryfile.h"
#include "arrayIO/hdf5.h"

#include "R_array_io.h"

using namespace Rcpp;
using namespace BPCells;



template<class T>
List dims_matrix(StoredMatrix<T> &&mat) {
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
        List l = dims_matrix(StoredMatrix<uint32_t>::openUnpacked(rb));
        l["compressed"] = false; l["type"] = "uint32_t";
        return l;
    } else if (version == "packed-uint-matrix-v1") {
        List l = dims_matrix(StoredMatrix<uint32_t>::openPacked(rb));
        l["compressed"] = true; l["type"] = "uint32_t";
        return l;
    } else if (version == "unpacked-float-matrix-v1") {
        List l = dims_matrix(StoredMatrix<float>::openUnpacked(rb));
        l["compressed"] = false; l["type"] = "float";
        return l;
    } else if (version == "packed-float-matrix-v1") {
        List l = dims_matrix(StoredMatrix<float>::openPacked(rb));
        l["compressed"] = true; l["type"] = "float";
        return l;
    } else if (version == "unpacked-double-matrix-v1") {
        List l = dims_matrix(StoredMatrix<double>::openUnpacked(rb));
        l["compressed"] = false; l["type"] = "double";
        return l;
    } else if (version == "packed-double-matrix-v1") {
        List l = dims_matrix(StoredMatrix<double>::openPacked(rb));
        l["compressed"] = true; l["type"] = "double";
        return l;
    } else {
        throw std::runtime_error(std::string("Matrix has unrecognized version ") + version);
    }
}


///////// MATRIX TRANSPOSE FUNCTIONS //////////

template<typename T>
SEXP write_matrix_transpose(SEXP matrix, std::string tmpdir, size_t load_bytes, size_t sort_buffer_bytes) {
    XPtr<MatrixLoader<T>> input(matrix);
    StoredMatrixTransposeWriter<T> transpose(tmpdir.c_str(), load_bytes, sort_buffer_bytes);
    transpose.write(*input);
    return Rcpp::wrap(
        XPtr<StoredMatrix<T>>(new StoredMatrix<T>(transpose.read()))
    );
}

// [[Rcpp::export]]
SEXP write_matrix_transpose_uint32_t_cpp(SEXP matrix, std::string tmpdir, size_t load_bytes, size_t sort_buffer_bytes) {
    return write_matrix_transpose<uint32_t>(matrix, tmpdir, load_bytes, sort_buffer_bytes);
}

// [[Rcpp::export]]
SEXP write_matrix_transpose_float_cpp(SEXP matrix, std::string tmpdir, size_t load_bytes, size_t sort_buffer_bytes) {
    return write_matrix_transpose<float>(matrix, tmpdir, load_bytes, sort_buffer_bytes);
}

// [[Rcpp::export]]
SEXP write_matrix_transpose_double_cpp(SEXP matrix, std::string tmpdir, size_t load_bytes, size_t sort_buffer_bytes) {
    return write_matrix_transpose<double>(matrix, tmpdir, load_bytes, sort_buffer_bytes);
}


///////// MEM MATRIX FUNCTIONS //////////
template<typename T>
SEXP iterate_packed_matrix(ReaderBuilder &rb, const StringVector row_names, const StringVector col_names) {
    return Rcpp::wrap(XPtr<MatrixLoader<T>>(new StoredMatrix<T>(
        StoredMatrix<T>::openPacked(rb, 1024,
            std::make_unique<RcppStringReader>(row_names), 
            std::make_unique<RcppStringReader>(col_names))
    )));
}

template<typename T>
SEXP iterate_unpacked_matrix(ReaderBuilder &rb, const StringVector row_names, const StringVector col_names) {
    return Rcpp::wrap(XPtr<MatrixLoader<T>>(new StoredMatrix<T>(
        StoredMatrix<T>::openUnpacked(rb,
            std::make_unique<RcppStringReader>(row_names), 
            std::make_unique<RcppStringReader>(col_names))
    )));
}

// [[Rcpp::export]]
SEXP iterate_packed_matrix_mem_uint32_t_cpp(S4 s4, const StringVector row_names, const StringVector col_names) {
    S4ReaderBuilder rb(s4);
    return iterate_packed_matrix<uint32_t>(rb, row_names, col_names);
}
// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_mem_uint32_t_cpp(S4 s4, const StringVector row_names, const StringVector col_names) {
    S4ReaderBuilder rb(s4);
    return iterate_unpacked_matrix<uint32_t>(rb, row_names, col_names);
}
// [[Rcpp::export]]
SEXP iterate_packed_matrix_mem_float_cpp(S4 s4, const StringVector row_names, const StringVector col_names) {
    S4ReaderBuilder rb(s4);
    return iterate_packed_matrix<float>(rb, row_names, col_names);
}
// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_mem_float_cpp(S4 s4, const StringVector row_names, const StringVector col_names) {
    S4ReaderBuilder rb(s4);
    return iterate_unpacked_matrix<float>(rb, row_names, col_names);
}
// [[Rcpp::export]]
SEXP iterate_packed_matrix_mem_double_cpp(S4 s4, const StringVector row_names, const StringVector col_names) {
    S4ReaderBuilder rb(s4);
    return iterate_packed_matrix<double>(rb, row_names, col_names);
}
// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_mem_double_cpp(S4 s4, const StringVector row_names, const StringVector col_names) {
    S4ReaderBuilder rb(s4);
    return iterate_unpacked_matrix<double>(rb, row_names, col_names);
}

template<typename T>
void write_packed_matrix(WriterBuilder &wb, SEXP matrix) {
    XPtr<MatrixLoader<T>> loader(matrix);
    loader->restart();
    StoredMatrixWriter<T>::createPacked(wb).write(*loader, &Rcpp::checkUserInterrupt);
}

template<typename T>
void write_unpacked_matrix(WriterBuilder &wb, SEXP matrix) {
    XPtr<MatrixLoader<T>> loader(matrix);
    loader->restart();
    StoredMatrixWriter<T>::createUnpacked(wb).write(*loader, &Rcpp::checkUserInterrupt);
}

// [[Rcpp::export]]
SEXP write_packed_matrix_mem_uint32_t_cpp(SEXP matrix) {
    ListWriterBuilder wb; write_packed_matrix<uint32_t>(wb, matrix); return wb.getList();
}

// [[Rcpp::export]]
SEXP write_unpacked_matrix_mem_uint32_t_cpp(SEXP matrix) {
    ListWriterBuilder wb; write_unpacked_matrix<uint32_t>(wb, matrix); return wb.getList();
}

// [[Rcpp::export]]
SEXP write_packed_matrix_mem_float_cpp(SEXP matrix) {
    ListWriterBuilder wb; write_packed_matrix<float>(wb, matrix); return wb.getList();
}

// [[Rcpp::export]]
SEXP write_unpacked_matrix_mem_float_cpp(SEXP matrix) {
    ListWriterBuilder wb; write_unpacked_matrix<float>(wb, matrix); return wb.getList();
}

// [[Rcpp::export]]
SEXP write_packed_matrix_mem_double_cpp(SEXP matrix) {
    ListWriterBuilder wb; write_packed_matrix<double>(wb, matrix); return wb.getList();
}

// [[Rcpp::export]]
SEXP write_unpacked_matrix_mem_double_cpp(SEXP matrix) {
    ListWriterBuilder wb; write_unpacked_matrix<double>(wb, matrix); return wb.getList();
}

///////// FILE MATRIX FUNCTIONS //////////

// [[Rcpp::export]]
List dims_matrix_file_cpp(std::string dir, uint32_t buffer_size) {
    FileReaderBuilder rb(dir, buffer_size);
    return dims_matrix_reader_builder(rb);
}

// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_file_uint32_t_cpp(std::string dir, uint32_t buffer_size, const StringVector row_names, const StringVector col_names) {
    FileReaderBuilder rb(dir, buffer_size);
    return iterate_unpacked_matrix<uint32_t>(rb, row_names, col_names);
}
// [[Rcpp::export]]
SEXP iterate_packed_matrix_file_uint32_t_cpp(std::string dir, uint32_t buffer_size, const StringVector row_names, const StringVector col_names) {
    FileReaderBuilder rb(dir, buffer_size);
    return iterate_packed_matrix<uint32_t>(rb, row_names, col_names);
}
// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_file_float_cpp(std::string dir, uint32_t buffer_size, const StringVector row_names, const StringVector col_names) {
    FileReaderBuilder rb(dir, buffer_size);
    return iterate_unpacked_matrix<float>(rb, row_names, col_names);
}
// [[Rcpp::export]]
SEXP iterate_packed_matrix_file_float_cpp(std::string dir, uint32_t buffer_size, const StringVector row_names, const StringVector col_names) {
    FileReaderBuilder rb(dir, buffer_size);
    return iterate_packed_matrix<float>(rb, row_names, col_names);
}
// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_file_double_cpp(std::string dir, uint32_t buffer_size, const StringVector row_names, const StringVector col_names) {
    FileReaderBuilder rb(dir, buffer_size);
    return iterate_unpacked_matrix<double>(rb, row_names, col_names);
}
// [[Rcpp::export]]
SEXP iterate_packed_matrix_file_double_cpp(std::string dir, uint32_t buffer_size, const StringVector row_names, const StringVector col_names) {
    FileReaderBuilder rb(dir, buffer_size);
    return iterate_packed_matrix<double>(rb, row_names, col_names);
}
// [[Rcpp::export]]
void write_unpacked_matrix_file_uint32_t_cpp(SEXP matrix, std::string dir, uint32_t buffer_size) {
    FileWriterBuilder wb(dir, buffer_size);    
    write_unpacked_matrix<uint32_t>(wb, matrix);
}
// [[Rcpp::export]]
void write_packed_matrix_file_uint32_t_cpp(SEXP matrix, std::string dir, uint32_t buffer_size) {
    FileWriterBuilder wb(dir, buffer_size);    
    write_packed_matrix<uint32_t>(wb, matrix);
}
// [[Rcpp::export]]
void write_unpacked_matrix_file_float_cpp(SEXP matrix, std::string dir, uint32_t buffer_size) {
    FileWriterBuilder wb(dir, buffer_size);    
    write_unpacked_matrix<float>(wb, matrix);
}
// [[Rcpp::export]]
void write_packed_matrix_file_float_cpp(SEXP matrix, std::string dir, uint32_t buffer_size) {
    FileWriterBuilder wb(dir, buffer_size);    
    write_packed_matrix<float>(wb, matrix);
}
// [[Rcpp::export]]
void write_unpacked_matrix_file_double_cpp(SEXP matrix, std::string dir, uint32_t buffer_size) {
    FileWriterBuilder wb(dir, buffer_size);    
    write_unpacked_matrix<double>(wb, matrix);
}
// [[Rcpp::export]]
void write_packed_matrix_file_double_cpp(SEXP matrix, std::string dir, uint32_t buffer_size) {
    FileWriterBuilder wb(dir, buffer_size);    
    write_packed_matrix<double>(wb, matrix);
}

///////// HDF5 MATRIX FUNCTIONS //////////

// [[Rcpp::export]]
List dims_matrix_hdf5_cpp(std::string file, std::string group, uint32_t buffer_size) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return dims_matrix_reader_builder(rb);
}

// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_hdf5_uint32_t_cpp(std::string file, std::string group, uint32_t buffer_size, const StringVector row_names, const StringVector col_names) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return iterate_unpacked_matrix<uint32_t>(rb, row_names, col_names);
}
// [[Rcpp::export]]
SEXP iterate_packed_matrix_hdf5_uint32_t_cpp(std::string file, std::string group, uint32_t buffer_size, const StringVector row_names, const StringVector col_names) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return iterate_packed_matrix<uint32_t>(rb, row_names, col_names);
}
// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_hdf5_float_cpp(std::string file, std::string group, uint32_t buffer_size, const StringVector row_names, const StringVector col_names) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return iterate_unpacked_matrix<float>(rb, row_names, col_names);
}
// [[Rcpp::export]]
SEXP iterate_packed_matrix_hdf5_float_cpp(std::string file, std::string group, uint32_t buffer_size, const StringVector row_names, const StringVector col_names) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return iterate_packed_matrix<float>(rb, row_names, col_names);
}
// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_hdf5_double_cpp(std::string file, std::string group, uint32_t buffer_size, const StringVector row_names, const StringVector col_names) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return iterate_unpacked_matrix<double>(rb, row_names, col_names);
}
// [[Rcpp::export]]
SEXP iterate_packed_matrix_hdf5_double_cpp(std::string file, std::string group, uint32_t buffer_size, const StringVector row_names, const StringVector col_names) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return iterate_packed_matrix<double>(rb, row_names, col_names);
}

// [[Rcpp::export]]
void write_unpacked_matrix_hdf5_uint32_t_cpp(SEXP matrix, std::string file, std::string group, uint32_t buffer_size, uint32_t chunk_size) {
    H5WriterBuilder wb(file, group, buffer_size, chunk_size);
    write_unpacked_matrix<uint32_t>(wb, matrix);
}
// [[Rcpp::export]]
void write_packed_matrix_hdf5_uint32_t_cpp(SEXP matrix, std::string file, std::string group, uint32_t buffer_size, uint32_t chunk_size) {
    H5WriterBuilder wb(file, group, buffer_size, chunk_size);
    write_packed_matrix<uint32_t>(wb, matrix);
}
// [[Rcpp::export]]
void write_unpacked_matrix_hdf5_float_cpp(SEXP matrix, std::string file, std::string group, uint32_t buffer_size, uint32_t chunk_size) {
    H5WriterBuilder wb(file, group, buffer_size, chunk_size);
    write_unpacked_matrix<float>(wb, matrix);
}
// [[Rcpp::export]]
void write_packed_matrix_hdf5_float_cpp(SEXP matrix, std::string file, std::string group, uint32_t buffer_size, uint32_t chunk_size) {
    H5WriterBuilder wb(file, group, buffer_size, chunk_size);
    write_packed_matrix<float>(wb, matrix);
}
// [[Rcpp::export]]
void write_unpacked_matrix_hdf5_double_cpp(SEXP matrix, std::string file, std::string group, uint32_t buffer_size, uint32_t chunk_size) {
    H5WriterBuilder wb(file, group, buffer_size, chunk_size);
    write_unpacked_matrix<double>(wb, matrix);
}
// [[Rcpp::export]]
void write_packed_matrix_hdf5_double_cpp(SEXP matrix, std::string file, std::string group, uint32_t buffer_size, uint32_t chunk_size) {
    H5WriterBuilder wb(file, group, buffer_size, chunk_size);
    write_packed_matrix<double>(wb, matrix);
}

///////// 3rd PARTY COMPATIBILITY HDF5 MATRIX FUNCTIONS //////////

// [[Rcpp::export]]
List dims_matrix_10x_hdf5_cpp(std::string file, uint32_t buffer_size) {
    return dims_matrix(open10xFeatureMatrix(file, buffer_size));
}

// [[Rcpp::export]]
SEXP iterate_matrix_10x_hdf5_cpp(std::string file, uint32_t buffer_size) {
    return Rcpp::wrap(
        XPtr<StoredMatrix<uint32_t>>(new StoredMatrix(open10xFeatureMatrix(file, buffer_size)))
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
        XPtr<StoredMatrix<float>>(new StoredMatrix(openAnnDataMatrix(file, group, buffer_size)))
    );
}