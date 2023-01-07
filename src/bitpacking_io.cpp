#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <Rcpp.h>

#include "arrayIO/binaryfile.h"
#include "arrayIO/bp128.h"
#include "arrayIO/hdf5.h"
#include "arrayIO/vector.h"

#include "lib/sleef_wrapper.h"

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

// [[Rcpp::export]]
std::string simd_vec_version() {
    switch (_SIMDBP128_MODE_) {
        case _SIMDBP128_X86_FULL: return "x86_full_simd";
        case _SIMDBP128_X86: return "x86_simd";
        case _SIMDBP128_ARM_NEON: return "arm_neon";
        case _SIMDBP128_C_FALLBACK: return "c_fallback";
    }
    return "Unknown";
}

// [[Rcpp::export]]
std::string simd_sleef_version() {
    switch (_SIMDBP128_MODE_) {
        case _BPCELLS_SLEEF_FALLBACK: return "FALLBACK";
        case _BPCELLS_SLEEF_SSE2: return "SSE2";
        case _BPCELLS_SLEEF_AVX: return "AVX";
        case _BPCELLS_SLEEF_AVX2: return "AVX2";
        case _BPCELLS_SLEEF_NEON: return "NEON";
    }
    return "Unknown";
}