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
List dims_matrix(StoredMatrix<T> &&mat, bool transpose) {
    IntegerVector dims(2);
    
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

    if (transpose) {
        dims[0] = mat.cols();
        dims[1] = mat.rows();
        return List::create(
            Named("dims") = dims,
            Named("row_names") = col_names,
            Named("col_names") = row_names,
            Named("transpose") = transpose
        );
    } else {
        dims[0] = mat.rows();
        dims[1] = mat.cols();
        return List::create(
            Named("dims") = dims,
            Named("row_names") = row_names,
            Named("col_names") = col_names,
            Named("transpose") = transpose
        );
    }
}

List dims_matrix_reader_builder(ReaderBuilder &rb) {
    std::string version = rb.readVersion();
    auto storage_order_reader = rb.openStringReader("storage_order");
    auto storage_order = storage_order_reader->get(0);
    bool row_major = storage_order == std::string_view("row");

    if (version == "unpacked-uint-matrix-v1") {
        List l = dims_matrix(StoredMatrix<uint32_t>::openUnpacked(rb), row_major);
        l["compressed"] = false; l["type"] = "uint32_t";
        return l;
    } else if (version == "packed-uint-matrix-v1") {
        List l = dims_matrix(StoredMatrix<uint32_t>::openPacked(rb), row_major);
        l["compressed"] = true; l["type"] = "uint32_t";
        return l;
    } else if (version == "unpacked-float-matrix-v1") {
        List l = dims_matrix(StoredMatrix<float>::openUnpacked(rb), row_major);
        l["compressed"] = false; l["type"] = "float";
        return l;
    } else if (version == "packed-float-matrix-v1") {
        List l = dims_matrix(StoredMatrix<float>::openPacked(rb), row_major);
        l["compressed"] = true; l["type"] = "float";
        return l;
    } else if (version == "unpacked-double-matrix-v1") {
        List l = dims_matrix(StoredMatrix<double>::openUnpacked(rb), row_major);
        l["compressed"] = false; l["type"] = "double";
        return l;
    } else if (version == "packed-double-matrix-v1") {
        List l = dims_matrix(StoredMatrix<double>::openPacked(rb), row_major);
        l["compressed"] = true; l["type"] = "double";
        return l;
    } else {
        throw std::runtime_error(std::string("Matrix has unrecognized version ") + version);
    }
}


///////// MATRIX TRANSPOSE FUNCTIONS //////////

template<typename T>
SEXP write_matrix_transpose_dir(SEXP matrix, std::string outdir, std::string tmpdir, uint64_t load_bytes, uint64_t sort_buffer_bytes, bool row_major) {
    XPtr<MatrixLoader<T>> input(matrix);
    FileWriterBuilder wb(outdir);
    StoredMatrixTransposeWriter<T> transpose(wb, tmpdir.c_str(), load_bytes, sort_buffer_bytes, row_major);
    transpose.write(*input);
    FileReaderBuilder rb(outdir);
    return Rcpp::wrap(XPtr<StoredMatrix<T>>(new StoredMatrix<T>(
        StoredMatrix<T>::openPacked(rb)
    )));
}

// [[Rcpp::export]]
SEXP write_matrix_transpose_uint32_t_cpp(SEXP matrix, std::string outdir, std::string tmpdir, uint64_t load_bytes, uint64_t sort_buffer_bytes, bool row_major) {
    return write_matrix_transpose_dir<uint32_t>(matrix, outdir, tmpdir, load_bytes, sort_buffer_bytes, row_major);
}

// [[Rcpp::export]]
SEXP write_matrix_transpose_float_cpp(SEXP matrix, std::string outdir, std::string tmpdir, uint64_t load_bytes, uint64_t sort_buffer_bytes, bool row_major) {
    return write_matrix_transpose_dir<float>(matrix, outdir, tmpdir, load_bytes, sort_buffer_bytes, row_major);
}

// [[Rcpp::export]]
SEXP write_matrix_transpose_double_cpp(SEXP matrix, std::string outdir, std::string tmpdir, uint64_t load_bytes, uint64_t sort_buffer_bytes, bool row_major) {
    return write_matrix_transpose_dir<double>(matrix, outdir, tmpdir, load_bytes, sort_buffer_bytes, row_major);
}


///////// MEM MATRIX FUNCTIONS //////////
template<typename T>
SEXP iterate_packed_matrix(ReaderBuilder &rb, const StringVector row_names, const StringVector col_names, uint32_t row_count) {
    return Rcpp::wrap(XPtr<MatrixLoader<T>>(new StoredMatrix<T>(
        StoredMatrix<T>::openPacked(rb, 1024,
            std::make_unique<RcppStringReader>(row_names), 
            std::make_unique<RcppStringReader>(col_names),
            row_count
        )
    )));
}

template<typename T>
SEXP iterate_unpacked_matrix(ReaderBuilder &rb, const StringVector row_names, const StringVector col_names, uint32_t row_count) {
    return Rcpp::wrap(XPtr<MatrixLoader<T>>(new StoredMatrix<T>(
        StoredMatrix<T>::openUnpacked(rb,
            std::make_unique<RcppStringReader>(row_names), 
            std::make_unique<RcppStringReader>(col_names),
            row_count
        )
    )));
}

// [[Rcpp::export]]
SEXP iterate_packed_matrix_mem_uint32_t_cpp(S4 s4, const StringVector row_names, const StringVector col_names, uint32_t row_count) {
    S4ReaderBuilder rb(s4);
    return iterate_packed_matrix<uint32_t>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_mem_uint32_t_cpp(S4 s4, const StringVector row_names, const StringVector col_names, uint32_t row_count) {
    S4ReaderBuilder rb(s4);
    return iterate_unpacked_matrix<uint32_t>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_packed_matrix_mem_float_cpp(S4 s4, const StringVector row_names, const StringVector col_names, uint32_t row_count) {
    S4ReaderBuilder rb(s4);
    return iterate_packed_matrix<float>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_mem_float_cpp(S4 s4, const StringVector row_names, const StringVector col_names, uint32_t row_count) {
    S4ReaderBuilder rb(s4);
    return iterate_unpacked_matrix<float>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_packed_matrix_mem_double_cpp(S4 s4, const StringVector row_names, const StringVector col_names, uint32_t row_count) {
    S4ReaderBuilder rb(s4);
    return iterate_packed_matrix<double>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_mem_double_cpp(S4 s4, const StringVector row_names, const StringVector col_names, uint32_t row_count) {
    S4ReaderBuilder rb(s4);
    return iterate_unpacked_matrix<double>(rb, row_names, col_names, row_count);
}

template<typename T>
void write_packed_matrix(WriterBuilder &wb, SEXP matrix, bool row_major) {
    XPtr<MatrixLoader<T>> loader(matrix);
    loader->restart();
    StoredMatrixWriter<T>::createPacked(wb, row_major).write(*loader, &Rcpp::checkUserInterrupt);
}

template<typename T>
void write_unpacked_matrix(WriterBuilder &wb, SEXP matrix, bool row_major) {
    XPtr<MatrixLoader<T>> loader(matrix);
    loader->restart();
    StoredMatrixWriter<T>::createUnpacked(wb, row_major).write(*loader, &Rcpp::checkUserInterrupt);
}

// [[Rcpp::export]]
SEXP write_packed_matrix_mem_uint32_t_cpp(SEXP matrix, bool row_major) {
    ListWriterBuilder wb; write_packed_matrix<uint32_t>(wb, matrix, row_major); return wb.getList();
}

// [[Rcpp::export]]
SEXP write_unpacked_matrix_mem_uint32_t_cpp(SEXP matrix, bool row_major) {
    ListWriterBuilder wb; write_unpacked_matrix<uint32_t>(wb, matrix, row_major); return wb.getList();
}

// [[Rcpp::export]]
SEXP write_packed_matrix_mem_float_cpp(SEXP matrix, bool row_major) {
    ListWriterBuilder wb; write_packed_matrix<float>(wb, matrix, row_major); return wb.getList();
}

// [[Rcpp::export]]
SEXP write_unpacked_matrix_mem_float_cpp(SEXP matrix, bool row_major) {
    ListWriterBuilder wb; write_unpacked_matrix<float>(wb, matrix, row_major); return wb.getList();
}

// [[Rcpp::export]]
SEXP write_packed_matrix_mem_double_cpp(SEXP matrix, bool row_major) {
    ListWriterBuilder wb; write_packed_matrix<double>(wb, matrix, row_major); return wb.getList();
}

// [[Rcpp::export]]
SEXP write_unpacked_matrix_mem_double_cpp(SEXP matrix, bool row_major) {
    ListWriterBuilder wb; write_unpacked_matrix<double>(wb, matrix, row_major); return wb.getList();
}

///////// FILE MATRIX FUNCTIONS //////////

// [[Rcpp::export]]
List dims_matrix_file_cpp(std::string dir, uint32_t buffer_size) {
    FileReaderBuilder rb(dir, buffer_size);
    return dims_matrix_reader_builder(rb);
}

// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_file_uint32_t_cpp(std::string dir, uint32_t buffer_size, const StringVector row_names, const StringVector col_names, uint32_t row_count) {
    FileReaderBuilder rb(dir, buffer_size);
    return iterate_unpacked_matrix<uint32_t>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_packed_matrix_file_uint32_t_cpp(std::string dir, uint32_t buffer_size, const StringVector row_names, const StringVector col_names, uint32_t row_count) {
    FileReaderBuilder rb(dir, buffer_size);
    return iterate_packed_matrix<uint32_t>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_file_float_cpp(std::string dir, uint32_t buffer_size, const StringVector row_names, const StringVector col_names, uint32_t row_count) {
    FileReaderBuilder rb(dir, buffer_size);
    return iterate_unpacked_matrix<float>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_packed_matrix_file_float_cpp(std::string dir, uint32_t buffer_size, const StringVector row_names, const StringVector col_names, uint32_t row_count) {
    FileReaderBuilder rb(dir, buffer_size);
    return iterate_packed_matrix<float>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_file_double_cpp(std::string dir, uint32_t buffer_size, const StringVector row_names, const StringVector col_names, uint32_t row_count) {
    FileReaderBuilder rb(dir, buffer_size);
    return iterate_unpacked_matrix<double>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_packed_matrix_file_double_cpp(std::string dir, uint32_t buffer_size, const StringVector row_names, const StringVector col_names, uint32_t row_count) {
    FileReaderBuilder rb(dir, buffer_size);
    return iterate_packed_matrix<double>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
void write_unpacked_matrix_file_uint32_t_cpp(SEXP matrix, std::string dir, uint32_t buffer_size, bool row_major) {
    FileWriterBuilder wb(dir, buffer_size);    
    write_unpacked_matrix<uint32_t>(wb, matrix, row_major);
}
// [[Rcpp::export]]
void write_packed_matrix_file_uint32_t_cpp(SEXP matrix, std::string dir, uint32_t buffer_size, bool row_major) {
    FileWriterBuilder wb(dir, buffer_size);    
    write_packed_matrix<uint32_t>(wb, matrix, row_major);
}
// [[Rcpp::export]]
void write_unpacked_matrix_file_float_cpp(SEXP matrix, std::string dir, uint32_t buffer_size, bool row_major) {
    FileWriterBuilder wb(dir, buffer_size);    
    write_unpacked_matrix<float>(wb, matrix, row_major);
}
// [[Rcpp::export]]
void write_packed_matrix_file_float_cpp(SEXP matrix, std::string dir, uint32_t buffer_size, bool row_major) {
    FileWriterBuilder wb(dir, buffer_size);    
    write_packed_matrix<float>(wb, matrix, row_major);
}
// [[Rcpp::export]]
void write_unpacked_matrix_file_double_cpp(SEXP matrix, std::string dir, uint32_t buffer_size, bool row_major) {
    FileWriterBuilder wb(dir, buffer_size);    
    write_unpacked_matrix<double>(wb, matrix, row_major);
}
// [[Rcpp::export]]
void write_packed_matrix_file_double_cpp(SEXP matrix, std::string dir, uint32_t buffer_size, bool row_major) {
    FileWriterBuilder wb(dir, buffer_size);    
    write_packed_matrix<double>(wb, matrix, row_major);
}

///////// HDF5 MATRIX FUNCTIONS //////////

// [[Rcpp::export]]
List dims_matrix_hdf5_cpp(std::string file, std::string group, uint32_t buffer_size) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return dims_matrix_reader_builder(rb);
}

// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_hdf5_uint32_t_cpp(std::string file, std::string group, uint32_t buffer_size, const StringVector row_names, const StringVector col_names, uint32_t row_count) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return iterate_unpacked_matrix<uint32_t>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_packed_matrix_hdf5_uint32_t_cpp(std::string file, std::string group, uint32_t buffer_size, const StringVector row_names, const StringVector col_names, uint32_t row_count) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return iterate_packed_matrix<uint32_t>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_hdf5_float_cpp(std::string file, std::string group, uint32_t buffer_size, const StringVector row_names, const StringVector col_names, uint32_t row_count) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return iterate_unpacked_matrix<float>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_packed_matrix_hdf5_float_cpp(std::string file, std::string group, uint32_t buffer_size, const StringVector row_names, const StringVector col_names, uint32_t row_count) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return iterate_packed_matrix<float>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_hdf5_double_cpp(std::string file, std::string group, uint32_t buffer_size, const StringVector row_names, const StringVector col_names, uint32_t row_count) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return iterate_unpacked_matrix<double>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_packed_matrix_hdf5_double_cpp(std::string file, std::string group, uint32_t buffer_size, const StringVector row_names, const StringVector col_names, uint32_t row_count) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return iterate_packed_matrix<double>(rb, row_names, col_names, row_count);
}

// [[Rcpp::export]]
void write_unpacked_matrix_hdf5_uint32_t_cpp(SEXP matrix, std::string file, std::string group, uint32_t buffer_size, uint32_t chunk_size, bool row_major) {
    H5WriterBuilder wb(file, group, buffer_size, chunk_size);
    write_unpacked_matrix<uint32_t>(wb, matrix, row_major);
}
// [[Rcpp::export]]
void write_packed_matrix_hdf5_uint32_t_cpp(SEXP matrix, std::string file, std::string group, uint32_t buffer_size, uint32_t chunk_size, bool row_major) {
    H5WriterBuilder wb(file, group, buffer_size, chunk_size);
    write_packed_matrix<uint32_t>(wb, matrix, row_major);
}
// [[Rcpp::export]]
void write_unpacked_matrix_hdf5_float_cpp(SEXP matrix, std::string file, std::string group, uint32_t buffer_size, uint32_t chunk_size, bool row_major) {
    H5WriterBuilder wb(file, group, buffer_size, chunk_size);
    write_unpacked_matrix<float>(wb, matrix, row_major);
}
// [[Rcpp::export]]
void write_packed_matrix_hdf5_float_cpp(SEXP matrix, std::string file, std::string group, uint32_t buffer_size, uint32_t chunk_size, bool row_major) {
    H5WriterBuilder wb(file, group, buffer_size, chunk_size);
    write_packed_matrix<float>(wb, matrix, row_major);
}
// [[Rcpp::export]]
void write_unpacked_matrix_hdf5_double_cpp(SEXP matrix, std::string file, std::string group, uint32_t buffer_size, uint32_t chunk_size, bool row_major) {
    H5WriterBuilder wb(file, group, buffer_size, chunk_size);
    write_unpacked_matrix<double>(wb, matrix, row_major);
}
// [[Rcpp::export]]
void write_packed_matrix_hdf5_double_cpp(SEXP matrix, std::string file, std::string group, uint32_t buffer_size, uint32_t chunk_size, bool row_major) {
    H5WriterBuilder wb(file, group, buffer_size, chunk_size);
    write_packed_matrix<double>(wb, matrix, row_major);
}

///////// 3rd PARTY COMPATIBILITY HDF5 MATRIX FUNCTIONS //////////

// [[Rcpp::export]]
List dims_matrix_10x_hdf5_cpp(std::string file, uint32_t buffer_size) {
    return dims_matrix(open10xFeatureMatrix(file, buffer_size), false);
}

// [[Rcpp::export]]
SEXP iterate_matrix_10x_hdf5_cpp(std::string file, uint32_t buffer_size) {
    return Rcpp::wrap(
        XPtr<StoredMatrix<uint32_t>>(new StoredMatrix(open10xFeatureMatrix(file, buffer_size)))
    );
}

// [[Rcpp::export]]
void write_matrix_10x_hdf5_cpp(
    SEXP matrix, 
    std::string path, 
    StringVector barcodes, 
    StringVector feature_ids, 
    StringVector feature_names, 
    StringVector feature_types, 
    List feature_metadata, 
    uint32_t buffer_size, 
    uint32_t chunk_size
) {
    XPtr<MatrixLoader<uint32_t>> loader(matrix);
    loader->restart();
    std::map<std::string, std::unique_ptr<StringReader>> metadata;
    StringVector metadata_names = feature_metadata.names();
    for (int i = 0; i < feature_metadata.size(); i++) {
        metadata[std::string(metadata_names[i])] = std::make_unique<RcppStringReader>(feature_metadata[i]);
    }
    StoredMatrixWriter<uint32_t> w = create10xFeatureMatrix(
        path, 
        RcppStringReader(barcodes), 
        RcppStringReader(feature_ids), 
        RcppStringReader(feature_names), 
        RcppStringReader(feature_types),
        metadata,
        buffer_size,
        chunk_size 
    );
    w.write(*loader, &Rcpp::checkUserInterrupt);
}

// [[Rcpp::export]]
List dims_matrix_anndata_hdf5_cpp(std::string file, std::string group, uint32_t buffer_size) {
    bool transpose = isRowOrientedAnnDataMatrix(file, group);
    List l = dims_matrix(openAnnDataMatrix(file, group, buffer_size), transpose);
    return l;
}

// [[Rcpp::export]]
SEXP iterate_matrix_anndata_hdf5_cpp(std::string file, std::string group, uint32_t buffer_size) {
    return Rcpp::wrap(
        XPtr<StoredMatrix<float>>(new StoredMatrix(openAnnDataMatrix(file, group, buffer_size)))
    );
}