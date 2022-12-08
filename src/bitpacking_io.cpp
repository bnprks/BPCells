#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <Rcpp.h>

#include "arrayIO/binaryfile.h"
#include "arrayIO/bp128.h"
#include "arrayIO/hdf5.h"
#include "arrayIO/vector.h"

using namespace Rcpp;
using namespace BPCells;

// [[Rcpp::export]]
IntegerVector read_integer_vector(SEXP input) {
    XPtr<UIntReader> reader(input);
    IntegerVector output(reader->size());
    for (size_t i = 0; i < reader->size(); i++) {
        output[i] = reader->read_one();
    }
    return output;
}

// [[Rcpp::export]]
SEXP open_file_reader(std::string path) {
    return Rcpp::wrap(
        XPtr<UIntReader>(new UIntReader(std::make_unique<FileUIntReader>(path.c_str()), 1024))
    );
}

// [[Rcpp::export]]
SEXP open_bp128_d1z(SEXP data, SEXP idx, SEXP starts, uint32_t count) {
    XPtr<UIntReader> data_ptr(data), idx_ptr(idx), starts_ptr(starts);
    return Rcpp::wrap(XPtr<UIntReader>(new UIntReader(
        std::make_unique<BP128_D1Z_UIntReader>(
            std::move(*data_ptr), std::move(*idx_ptr), std::move(*starts_ptr), count
        ),
        1024,
        1024
    )));
}

// [[Rcpp::export]]
SEXP open_bp128_d1(SEXP data, SEXP idx, SEXP starts, uint32_t count) {
    XPtr<UIntReader> data_ptr(data), idx_ptr(idx), starts_ptr(starts);
    return Rcpp::wrap(XPtr<UIntReader>(new UIntReader(
        std::make_unique<BP128_D1_UIntReader>(
            std::move(*data_ptr), std::move(*idx_ptr), std::move(*starts_ptr), count
        ),
        1024,
        1024
    )));
}

// [[Rcpp::export]]
SEXP open_bp128_for(SEXP data, SEXP idx, uint32_t count) {
    XPtr<UIntReader> data_ptr(data), idx_ptr(idx);
    return Rcpp::wrap(XPtr<UIntReader>(new UIntReader(
        std::make_unique<BP128_FOR_UIntReader>(std::move(*data_ptr), std::move(*idx_ptr), count),
        1024,
        1024
    )));
}
