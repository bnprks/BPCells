#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <Rcpp.h>

#include "matrixIterators/ImportMatrixHDF5.h"
#include "matrixIterators/MatrixIterator.h"
#include "matrixIterators/MatrixMarketImport.h"
#include "matrixIterators/StoredMatrix.h"
#include "matrixIterators/StoredMatrixTransposeWriter.h"
#include "matrixIterators/StoredMatrixWriter.h"

#include "arrayIO/binaryfile.h"
#include "arrayIO/hdf5.h"
#include "arrayIO/vector.h"

#include "R_array_io.h"
#include "R_interrupts.h"
#include "R_xptr_wrapper.h"

using namespace Rcpp;
using namespace BPCells;

template <class T> List dims_matrix(StoredMatrix<T> &&mat, bool transpose) {
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
        if (mat.colNames(i) == NULL)
            throw std::runtime_error("Matrix has some colnames, but not equal to col number");
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

    if (version == "unpacked-uint-matrix-v1" || version == "unpacked-uint-matrix-v2") {
        List l = dims_matrix(StoredMatrix<uint32_t>::openUnpacked(rb), row_major);
        l["compressed"] = false;
        l["type"] = "uint32_t";
        return l;
    } else if (version == "packed-uint-matrix-v1" || version == "packed-uint-matrix-v2") {
        List l = dims_matrix(StoredMatrix<uint32_t>::openPacked(rb), row_major);
        l["compressed"] = true;
        l["type"] = "uint32_t";
        return l;
    } else if (version == "unpacked-float-matrix-v1" || version == "unpacked-float-matrix-v2") {
        List l = dims_matrix(StoredMatrix<float>::openUnpacked(rb), row_major);
        l["compressed"] = false;
        l["type"] = "float";
        return l;
    } else if (version == "packed-float-matrix-v1" || version == "packed-float-matrix-v2") {
        List l = dims_matrix(StoredMatrix<float>::openPacked(rb), row_major);
        l["compressed"] = true;
        l["type"] = "float";
        return l;
    } else if (version == "unpacked-double-matrix-v1" || version == "unpacked-double-matrix-v2") {
        List l = dims_matrix(StoredMatrix<double>::openUnpacked(rb), row_major);
        l["compressed"] = false;
        l["type"] = "double";
        return l;
    } else if (version == "packed-double-matrix-v1" || version == "packed-double-matrix-v2") {
        List l = dims_matrix(StoredMatrix<double>::openPacked(rb), row_major);
        l["compressed"] = true;
        l["type"] = "double";
        return l;
    } else {
        throw std::runtime_error(std::string("Matrix has unrecognized version ") + version);
    }
}

///////// MATRIX TRANSPOSE FUNCTIONS //////////

template <typename T>
SEXP write_matrix_transpose_dir(
    SEXP matrix,
    std::string outdir,
    std::string tmpdir,
    uint64_t load_bytes,
    uint64_t sort_buffer_bytes,
    bool row_major
) {
    FileWriterBuilder wb(outdir);
    StoredMatrixTransposeWriter<T> transpose(
        wb, tmpdir.c_str(), load_bytes, sort_buffer_bytes, row_major
    );
    auto mat = take_unique_xptr<MatrixLoader<T>>(matrix);
    run_with_R_interrupt_check(&MatrixWriter<T>::write, &transpose, std::ref(*mat));

    FileReaderBuilder rb(outdir);
    return make_unique_xptr<StoredMatrix<T>>(StoredMatrix<T>::openPacked(rb));
}

// [[Rcpp::export]]
SEXP write_matrix_transpose_uint32_t_cpp(
    SEXP matrix,
    std::string outdir,
    std::string tmpdir,
    uint64_t load_bytes,
    uint64_t sort_buffer_bytes,
    bool row_major
) {
    return write_matrix_transpose_dir<uint32_t>(
        matrix, outdir, tmpdir, load_bytes, sort_buffer_bytes, row_major
    );
}

// [[Rcpp::export]]
SEXP write_matrix_transpose_float_cpp(
    SEXP matrix,
    std::string outdir,
    std::string tmpdir,
    uint64_t load_bytes,
    uint64_t sort_buffer_bytes,
    bool row_major
) {
    return write_matrix_transpose_dir<float>(
        matrix, outdir, tmpdir, load_bytes, sort_buffer_bytes, row_major
    );
}

// [[Rcpp::export]]
SEXP write_matrix_transpose_double_cpp(
    SEXP matrix,
    std::string outdir,
    std::string tmpdir,
    uint64_t load_bytes,
    uint64_t sort_buffer_bytes,
    bool row_major
) {
    return write_matrix_transpose_dir<double>(
        matrix, outdir, tmpdir, load_bytes, sort_buffer_bytes, row_major
    );
}

///////// MEM MATRIX FUNCTIONS //////////
template <typename T>
SEXP iterate_packed_matrix(
    ReaderBuilder &rb,
    const StringVector row_names,
    const StringVector col_names,
    uint32_t row_count
) {
    return make_unique_xptr<StoredMatrix<T>>(StoredMatrix<T>::openPacked(
        rb,
        1024,
        std::make_unique<RcppStringReader>(row_names),
        std::make_unique<RcppStringReader>(col_names),
        row_count
    ));
}

template <typename T>
SEXP iterate_unpacked_matrix(
    ReaderBuilder &rb,
    const StringVector row_names,
    const StringVector col_names,
    uint32_t row_count
) {
    return make_unique_xptr<StoredMatrix<T>>(StoredMatrix<T>::openUnpacked(
        rb,
        std::make_unique<RcppStringReader>(row_names),
        std::make_unique<RcppStringReader>(col_names),
        row_count
    ));
}

// [[Rcpp::export]]
SEXP iterate_packed_matrix_mem_uint32_t_cpp(
    S4 s4, const StringVector row_names, const StringVector col_names, uint32_t row_count
) {
    S4ReaderBuilder rb(s4);
    return iterate_packed_matrix<uint32_t>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_mem_uint32_t_cpp(
    S4 s4, const StringVector row_names, const StringVector col_names, uint32_t row_count
) {
    S4ReaderBuilder rb(s4);
    return iterate_unpacked_matrix<uint32_t>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_packed_matrix_mem_float_cpp(
    S4 s4, const StringVector row_names, const StringVector col_names, uint32_t row_count
) {
    S4ReaderBuilder rb(s4);
    return iterate_packed_matrix<float>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_mem_float_cpp(
    S4 s4, const StringVector row_names, const StringVector col_names, uint32_t row_count
) {
    S4ReaderBuilder rb(s4);
    return iterate_unpacked_matrix<float>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_packed_matrix_mem_double_cpp(
    S4 s4, const StringVector row_names, const StringVector col_names, uint32_t row_count
) {
    S4ReaderBuilder rb(s4);
    return iterate_packed_matrix<double>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_mem_double_cpp(
    S4 s4, const StringVector row_names, const StringVector col_names, uint32_t row_count
) {
    S4ReaderBuilder rb(s4);
    return iterate_unpacked_matrix<double>(rb, row_names, col_names, row_count);
}

template <typename T> void write_packed_matrix(WriterBuilder &wb, SEXP matrix, bool row_major) {
    auto loader = take_unique_xptr<MatrixLoader<T>>(matrix);
    loader->restart();
    run_with_R_interrupt_check(
        &StoredMatrixWriter<T>::write,
        StoredMatrixWriter<T>::createPacked(wb, row_major),
        std::ref(*loader)
    );
}

template <typename T> void write_unpacked_matrix(WriterBuilder &wb, SEXP matrix, bool row_major) {
    auto loader = take_unique_xptr<MatrixLoader<T>>(matrix);
    loader->restart();
    run_with_R_interrupt_check(
        &StoredMatrixWriter<T>::write,
        StoredMatrixWriter<T>::createUnpacked(wb, row_major),
        std::ref(*loader)
    );
}

// [[Rcpp::export]]
SEXP write_packed_matrix_mem_uint32_t_cpp(SEXP matrix, bool row_major) {
    ListWriterBuilder wb;
    write_packed_matrix<uint32_t>(wb, matrix, row_major);
    return wb.getList();
}

// [[Rcpp::export]]
SEXP write_unpacked_matrix_mem_uint32_t_cpp(SEXP matrix, bool row_major) {
    ListWriterBuilder wb;
    write_unpacked_matrix<uint32_t>(wb, matrix, row_major);
    return wb.getList();
}

// [[Rcpp::export]]
SEXP write_packed_matrix_mem_float_cpp(SEXP matrix, bool row_major) {
    ListWriterBuilder wb;
    write_packed_matrix<float>(wb, matrix, row_major);
    return wb.getList();
}

// [[Rcpp::export]]
SEXP write_unpacked_matrix_mem_float_cpp(SEXP matrix, bool row_major) {
    ListWriterBuilder wb;
    write_unpacked_matrix<float>(wb, matrix, row_major);
    return wb.getList();
}

// [[Rcpp::export]]
SEXP write_packed_matrix_mem_double_cpp(SEXP matrix, bool row_major) {
    ListWriterBuilder wb;
    write_packed_matrix<double>(wb, matrix, row_major);
    return wb.getList();
}

// [[Rcpp::export]]
SEXP write_unpacked_matrix_mem_double_cpp(SEXP matrix, bool row_major) {
    ListWriterBuilder wb;
    write_unpacked_matrix<double>(wb, matrix, row_major);
    return wb.getList();
}

///////// FILE MATRIX FUNCTIONS //////////

// [[Rcpp::export]]
List dims_matrix_file_cpp(std::string dir, uint32_t buffer_size) {
    FileReaderBuilder rb(dir, buffer_size);
    return dims_matrix_reader_builder(rb);
}

// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_file_uint32_t_cpp(
    std::string dir,
    uint32_t buffer_size,
    const StringVector row_names,
    const StringVector col_names,
    uint32_t row_count
) {
    FileReaderBuilder rb(dir, buffer_size);
    return iterate_unpacked_matrix<uint32_t>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_packed_matrix_file_uint32_t_cpp(
    std::string dir,
    uint32_t buffer_size,
    const StringVector row_names,
    const StringVector col_names,
    uint32_t row_count
) {
    FileReaderBuilder rb(dir, buffer_size);
    return iterate_packed_matrix<uint32_t>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_file_float_cpp(
    std::string dir,
    uint32_t buffer_size,
    const StringVector row_names,
    const StringVector col_names,
    uint32_t row_count
) {
    FileReaderBuilder rb(dir, buffer_size);
    return iterate_unpacked_matrix<float>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_packed_matrix_file_float_cpp(
    std::string dir,
    uint32_t buffer_size,
    const StringVector row_names,
    const StringVector col_names,
    uint32_t row_count
) {
    FileReaderBuilder rb(dir, buffer_size);
    return iterate_packed_matrix<float>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_file_double_cpp(
    std::string dir,
    uint32_t buffer_size,
    const StringVector row_names,
    const StringVector col_names,
    uint32_t row_count
) {
    FileReaderBuilder rb(dir, buffer_size);
    return iterate_unpacked_matrix<double>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_packed_matrix_file_double_cpp(
    std::string dir,
    uint32_t buffer_size,
    const StringVector row_names,
    const StringVector col_names,
    uint32_t row_count
) {
    FileReaderBuilder rb(dir, buffer_size);
    return iterate_packed_matrix<double>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
void write_unpacked_matrix_file_uint32_t_cpp(
    SEXP matrix, std::string dir, uint32_t buffer_size, bool allow_overwrite, bool row_major
) {
    FileWriterBuilder wb(dir, buffer_size, allow_overwrite);
    write_unpacked_matrix<uint32_t>(wb, matrix, row_major);
}
// [[Rcpp::export]]
void write_packed_matrix_file_uint32_t_cpp(
    SEXP matrix, std::string dir, uint32_t buffer_size, bool allow_overwrite, bool row_major
) {
    FileWriterBuilder wb(dir, buffer_size, allow_overwrite);
    write_packed_matrix<uint32_t>(wb, matrix, row_major);
}
// [[Rcpp::export]]
void write_unpacked_matrix_file_float_cpp(
    SEXP matrix, std::string dir, uint32_t buffer_size, bool allow_overwrite, bool row_major
) {
    FileWriterBuilder wb(dir, buffer_size, allow_overwrite);
    write_unpacked_matrix<float>(wb, matrix, row_major);
}
// [[Rcpp::export]]
void write_packed_matrix_file_float_cpp(
    SEXP matrix, std::string dir, uint32_t buffer_size, bool allow_overwrite, bool row_major
) {
    FileWriterBuilder wb(dir, buffer_size, allow_overwrite);
    write_packed_matrix<float>(wb, matrix, row_major);
}
// [[Rcpp::export]]
void write_unpacked_matrix_file_double_cpp(
    SEXP matrix, std::string dir, uint32_t buffer_size, bool allow_overwrite, bool row_major
) {
    FileWriterBuilder wb(dir, buffer_size, allow_overwrite);
    write_unpacked_matrix<double>(wb, matrix, row_major);
}
// [[Rcpp::export]]
void write_packed_matrix_file_double_cpp(
    SEXP matrix, std::string dir, uint32_t buffer_size, bool allow_overwrite, bool row_major
) {
    FileWriterBuilder wb(dir, buffer_size, allow_overwrite);
    write_packed_matrix<double>(wb, matrix, row_major);
}

///////// HDF5 MATRIX FUNCTIONS //////////

// [[Rcpp::export]]
List dims_matrix_hdf5_cpp(std::string file, std::string group, uint32_t buffer_size) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return dims_matrix_reader_builder(rb);
}

// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_hdf5_uint32_t_cpp(
    std::string file,
    std::string group,
    uint32_t buffer_size,
    const StringVector row_names,
    const StringVector col_names,
    uint32_t row_count
) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return iterate_unpacked_matrix<uint32_t>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_packed_matrix_hdf5_uint32_t_cpp(
    std::string file,
    std::string group,
    uint32_t buffer_size,
    const StringVector row_names,
    const StringVector col_names,
    uint32_t row_count
) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return iterate_packed_matrix<uint32_t>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_hdf5_float_cpp(
    std::string file,
    std::string group,
    uint32_t buffer_size,
    const StringVector row_names,
    const StringVector col_names,
    uint32_t row_count
) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return iterate_unpacked_matrix<float>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_packed_matrix_hdf5_float_cpp(
    std::string file,
    std::string group,
    uint32_t buffer_size,
    const StringVector row_names,
    const StringVector col_names,
    uint32_t row_count
) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return iterate_packed_matrix<float>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_unpacked_matrix_hdf5_double_cpp(
    std::string file,
    std::string group,
    uint32_t buffer_size,
    const StringVector row_names,
    const StringVector col_names,
    uint32_t row_count
) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return iterate_unpacked_matrix<double>(rb, row_names, col_names, row_count);
}
// [[Rcpp::export]]
SEXP iterate_packed_matrix_hdf5_double_cpp(
    std::string file,
    std::string group,
    uint32_t buffer_size,
    const StringVector row_names,
    const StringVector col_names,
    uint32_t row_count
) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return iterate_packed_matrix<double>(rb, row_names, col_names, row_count);
}

// [[Rcpp::export]]
void write_unpacked_matrix_hdf5_uint32_t_cpp(
    SEXP matrix,
    std::string file,
    std::string group,
    uint32_t buffer_size,
    uint32_t chunk_size,
    bool allow_overwrite,
    bool row_major,
    uint32_t gzip_level
) {
    H5WriterBuilder wb(file, group, buffer_size, chunk_size, allow_overwrite, gzip_level);
    write_unpacked_matrix<uint32_t>(wb, matrix, row_major);
}
// [[Rcpp::export]]
void write_packed_matrix_hdf5_uint32_t_cpp(
    SEXP matrix,
    std::string file,
    std::string group,
    uint32_t buffer_size,
    uint32_t chunk_size,
    bool allow_overwrite,
    bool row_major,
    uint32_t gzip_level
) {
    H5WriterBuilder wb(file, group, buffer_size, chunk_size, allow_overwrite, gzip_level);
    write_packed_matrix<uint32_t>(wb, matrix, row_major);
}
// [[Rcpp::export]]
void write_unpacked_matrix_hdf5_float_cpp(
    SEXP matrix,
    std::string file,
    std::string group,
    uint32_t buffer_size,
    uint32_t chunk_size,
    bool allow_overwrite,
    bool row_major,
    uint32_t gzip_level
) {
    H5WriterBuilder wb(file, group, buffer_size, chunk_size, allow_overwrite, gzip_level);
    write_unpacked_matrix<float>(wb, matrix, row_major);
}
// [[Rcpp::export]]
void write_packed_matrix_hdf5_float_cpp(
    SEXP matrix,
    std::string file,
    std::string group,
    uint32_t buffer_size,
    uint32_t chunk_size,
    bool allow_overwrite,
    bool row_major,
    uint32_t gzip_level
) {
    H5WriterBuilder wb(file, group, buffer_size, chunk_size, allow_overwrite, gzip_level);
    write_packed_matrix<float>(wb, matrix, row_major);
}
// [[Rcpp::export]]
void write_unpacked_matrix_hdf5_double_cpp(
    SEXP matrix,
    std::string file,
    std::string group,
    uint32_t buffer_size,
    uint32_t chunk_size,
    bool allow_overwrite,
    bool row_major,
    uint32_t gzip_level
) {
    H5WriterBuilder wb(file, group, buffer_size, chunk_size, allow_overwrite, gzip_level);
    write_unpacked_matrix<double>(wb, matrix, row_major);
}
// [[Rcpp::export]]
void write_packed_matrix_hdf5_double_cpp(
    SEXP matrix,
    std::string file,
    std::string group,
    uint32_t buffer_size,
    uint32_t chunk_size,
    bool allow_overwrite,
    bool row_major,
    uint32_t gzip_level
) {
    H5WriterBuilder wb(file, group, buffer_size, chunk_size, allow_overwrite, gzip_level);
    write_packed_matrix<double>(wb, matrix, row_major);
}

///////// 3rd PARTY COMPATIBILITY HDF5 MATRIX FUNCTIONS //////////

// 10x matrix //////////

// [[Rcpp::export]]
List dims_matrix_10x_hdf5_cpp(std::string file, std::string group, uint32_t buffer_size) {
    std::string type = get10xMatrixType(file, group);
    List l;
    if (type == "uint32_t") {
        l = dims_matrix(open10xFeatureMatrix<uint32_t>(file, group, buffer_size), false);
    } else if (type == "uint64_t") {
        l = dims_matrix(open10xFeatureMatrix<uint64_t>(file, group, buffer_size), false);
    } else if (type == "float") {
        l = dims_matrix(open10xFeatureMatrix<float>(file, group, buffer_size), false);
    } else if (type == "double") {
        l = dims_matrix(open10xFeatureMatrix<double>(file, group, buffer_size), false);
    } else {
        throw std::runtime_error("dims_matrix_10x_hdf5_cpp: Unrecognized matrix type " + type);
    }
    l["type"] = type;
    return l;
}

// [[Rcpp::export]]
SEXP iterate_matrix_10x_hdf5_cpp(
    std::string file,
    std::string group,
    uint32_t buffer_size,
    const StringVector row_names,
    const StringVector col_names
) {
    std::string type = get10xMatrixType(file, group);
    if (type == "uint32_t") {
        return make_unique_xptr<StoredMatrix<uint32_t>>(open10xFeatureMatrix<uint32_t>(
            file,
            group,
            buffer_size,
            std::make_unique<RcppStringReader>(row_names),
            std::make_unique<RcppStringReader>(col_names)
        ));
    } else if (type == "uint64_t") {
        return make_unique_xptr<StoredMatrix<uint64_t>>(open10xFeatureMatrix<uint64_t>(
            file,
            group,
            buffer_size,
            std::make_unique<RcppStringReader>(row_names),
            std::make_unique<RcppStringReader>(col_names)
        ));
    } else if (type == "float") {
        return make_unique_xptr<StoredMatrix<float>>(open10xFeatureMatrix<float>(
            file,
            group,
            buffer_size,
            std::make_unique<RcppStringReader>(row_names),
            std::make_unique<RcppStringReader>(col_names)
        ));
    } else if (type == "double") {
        return make_unique_xptr<StoredMatrix<double>>(open10xFeatureMatrix<double>(
            file,
            group,
            buffer_size,
            std::make_unique<RcppStringReader>(row_names),
            std::make_unique<RcppStringReader>(col_names)
        ));
    } else {
        throw std::runtime_error("iterate_matrix_10x_hdf5_cpp: Unrecognized matrix type " + type);
    }
}

template <typename T>
void write_matrix_10x_hdf5_base(
    SEXP matrix,
    std::string path,
    StringVector barcodes,
    StringVector feature_ids,
    StringVector feature_names,
    StringVector feature_types,
    List feature_metadata,
    uint32_t buffer_size,
    uint32_t chunk_size,
    uint32_t gzip_level
) {
    auto loader = take_unique_xptr<MatrixLoader<T>>(matrix);
    loader->restart();
    std::map<std::string, std::unique_ptr<StringReader>> metadata;
    StringVector metadata_names = feature_metadata.names();
    for (int i = 0; i < feature_metadata.size(); i++) {
        metadata[std::string(metadata_names[i])] =
            std::make_unique<RcppStringReader>(Rcpp::as<StringVector>(feature_metadata[i]));
    }
    StoredMatrixWriter<T> w = create10xFeatureMatrix<T>(
        path,
        RcppStringReader(barcodes),
        RcppStringReader(feature_ids),
        RcppStringReader(feature_names),
        RcppStringReader(feature_types),
        metadata,
        buffer_size,
        chunk_size,
        gzip_level
    );
    run_with_R_interrupt_check(&StoredMatrixWriter<T>::write, &w, std::ref(*loader));
}

// [[Rcpp::export]]
void write_matrix_10x_hdf5_cpp(
    SEXP matrix,
    std::string path,
    std::string type,
    StringVector barcodes,
    StringVector feature_ids,
    StringVector feature_names,
    StringVector feature_types,
    List feature_metadata,
    uint32_t buffer_size,
    uint32_t chunk_size,
    uint32_t gzip_level
) {
    if (type == "uint32_t") {
        write_matrix_10x_hdf5_base<uint32_t>(
            matrix,
            path,
            barcodes,
            feature_ids,
            feature_names,
            feature_types,
            feature_metadata,
            buffer_size,
            chunk_size,
            gzip_level
        );
    } else if (type == "uint64_t") {
        write_matrix_10x_hdf5_base<uint64_t>(
            matrix,
            path,
            barcodes,
            feature_ids,
            feature_names,
            feature_types,
            feature_metadata,
            buffer_size,
            chunk_size,
            gzip_level
        );
    } else if (type == "float") {
        write_matrix_10x_hdf5_base<float>(
            matrix,
            path,
            barcodes,
            feature_ids,
            feature_names,
            feature_types,
            feature_metadata,
            buffer_size,
            chunk_size,
            gzip_level
        );
    } else if (type == "double") {
        write_matrix_10x_hdf5_base<double>(
            matrix,
            path,
            barcodes,
            feature_ids,
            feature_names,
            feature_types,
            feature_metadata,
            buffer_size,
            chunk_size,
            gzip_level
        );
    } else {
        throw std::runtime_error("write_matrix_10x_hdf5_cpp: unsupported type " + type);
    }
}

// AnnData //////////

// [[Rcpp::export]]
List dims_matrix_anndata_hdf5_cpp(std::string file, std::string group, uint32_t buffer_size) {
    bool transpose = isRowOrientedAnnDataMatrix(file, group);
    std::string type = getAnnDataMatrixType(file, group);
    List l;
    if (type == "uint32_t") {
        l = dims_matrix(openAnnDataMatrix<uint32_t>(file, group, buffer_size), transpose);
    } else if (type == "float") {
        l = dims_matrix(openAnnDataMatrix<float>(file, group, buffer_size), transpose);
    } else if (type == "double") {
        l = dims_matrix(openAnnDataMatrix<double>(file, group, buffer_size), transpose);
    } else {
        throw std::runtime_error("dims_matrix_anndata_hdf5_cpp: Unrecognized matrix type " + type);
    }
    l["type"] = type;
    return l;
}

// [[Rcpp::export]]
SEXP iterate_matrix_anndata_hdf5_cpp(
    std::string file,
    std::string group,
    std::string type,
    uint32_t buffer_size,
    const StringVector row_names,
    const StringVector col_names
) {
    std::string file_type = getAnnDataMatrixType(file, group);
    if (type != file_type) {
        std::stringstream ss;
        ss << "iterate_matrix_anndata_hdf5_cpp: Found type '" << file_type << "' expected '" << type
           << "'";
        throw std::runtime_error(ss.str());
    }
    if (type == "uint32_t") {
        return make_unique_xptr<StoredMatrix<uint32_t>>(openAnnDataMatrix<uint32_t>(
            file,
            group,
            buffer_size,
            std::make_unique<RcppStringReader>(row_names),
            std::make_unique<RcppStringReader>(col_names)
        ));
    } else if (type == "float") {
        return make_unique_xptr<StoredMatrix<float>>(openAnnDataMatrix<float>(
            file,
            group,
            buffer_size,
            std::make_unique<RcppStringReader>(row_names),
            std::make_unique<RcppStringReader>(col_names)
        ));
    } else if (type == "double") {
        return make_unique_xptr<StoredMatrix<double>>(openAnnDataMatrix<double>(
            file,
            group,
            buffer_size,
            std::make_unique<RcppStringReader>(row_names),
            std::make_unique<RcppStringReader>(col_names)
        ));
    } else {
        throw std::runtime_error("iterate_matrix_anndata_hdf5_cpp: Unsupported type " + type);
    }
}

template <typename T>
void write_matrix_anndata_hdf5_base(
    SEXP matrix,
    std::string file,
    std::string group,
    bool row_major,
    uint32_t buffer_size,
    uint32_t chunk_size,
    uint32_t gzip_level
) {
    auto loader = take_unique_xptr<MatrixLoader<T>>(matrix);
    loader->restart();

    StoredMatrixWriter<T> w =
        createAnnDataMatrix<T>(file, group, row_major, buffer_size, chunk_size, gzip_level);
    run_with_R_interrupt_check(&StoredMatrixWriter<T>::write, &w, std::ref(*loader));
    createAnnDataObsVarIfMissing(*loader, file, row_major, gzip_level);
}

// [[Rcpp::export]]
void write_matrix_anndata_hdf5_cpp(
    SEXP matrix,
    std::string file,
    std::string group,
    std::string type,
    bool row_major,
    uint32_t buffer_size,
    uint32_t chunk_size,
    uint32_t gzip_level
) {
    if (type == "uint32_t") {
        write_matrix_anndata_hdf5_base<uint32_t>(
            matrix, file, group, row_major, buffer_size, chunk_size, gzip_level
        );
    } else if (type == "float") {
        write_matrix_anndata_hdf5_base<float>(
            matrix, file, group, row_major, buffer_size, chunk_size, gzip_level
        );
    } else if (type == "double") {
        write_matrix_anndata_hdf5_base<double>(
            matrix, file, group, row_major, buffer_size, chunk_size, gzip_level
        );
    } else {
        throw std::runtime_error("write_matrix_anndata_hdf5_cpp: unsupported type " + type);
    }
}

// [[Rcpp::export]]
std::vector<std::string>
read_hdf5_string_cpp(std::string path, std::string group, uint32_t buffer_size) {
    H5ReaderBuilder rb(path, "/", buffer_size);
    auto reader = rb.openStringReader(group);
    std::vector<std::string> res;
    for (uint32_t i = 0; i < reader->size(); i++) {
        res.push_back(reader->get(i));
    }
    return res;
}

// [[Rcpp::export]]
bool hdf5_group_exists_cpp(std::string path, std::string group) {
    H5ReaderBuilder rb(path, "/", 1);
    return rb.getGroup().exist(group);
}

// [[Rcpp::export]]
std::vector<std::string> hdf5_group_objnames_cpp(std::string path, std::string group) {
    H5ReaderBuilder rb(path, group, 1);
    return rb.getGroup().listObjectNames();
}

// MTX format

// [[Rcpp::export]]
void import_matrix_market_cpp(
    std::string mtx_path,
    std::vector<std::string> row_names,
    std::vector<std::string> col_names,
    std::string outdir,
    std::string tmpdir,
    uint64_t load_bytes,
    uint64_t sort_buffer_bytes,
    bool row_major
) {
    FileWriterBuilder wb(outdir);

    run_with_R_interrupt_check(
        importMtx,
        mtx_path,
        std::move(row_names),
        std::move(col_names),
        std::ref(wb),
        tmpdir.c_str(),
        load_bytes,
        sort_buffer_bytes,
        row_major
    );
}