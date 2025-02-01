// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <Rcpp.h>

#include "bpcells-cpp/arrayIO/binaryfile.h"
#include "bpcells-cpp/arrayIO/bp128.h"
#include "bpcells-cpp/arrayIO/hdf5.h"
#include "bpcells-cpp/arrayIO/vector.h"

#include "bpcells-cpp/simd/math.h"
#include "bpcells-cpp/simd/current_target.h"
#include "bpcells-cpp/simd/bp128.h"

#include <hwy/aligned_allocator.h>

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
SEXP open_bp128_d1z(SEXP data, SEXP idx, SEXP idx_offsets, SEXP starts, uint32_t count) {
    XPtr<UIntReader> data_ptr(data), idx_ptr(idx), starts_ptr(starts);
    XPtr<ULongReader> idx_offsets_ptr(idx_offsets);
    return Rcpp::wrap(XPtr<UIntReader>(new UIntReader(
        std::make_unique<BP128_D1Z_UIntReader>(
            std::move(*data_ptr), std::move(*idx_ptr), std::move(*idx_offsets_ptr), std::move(*starts_ptr), count
        ),
        1024,
        1024
    )));
}

// [[Rcpp::export]]
SEXP open_bp128_d1(SEXP data, SEXP idx, SEXP idx_offsets, SEXP starts, uint32_t count) {
    XPtr<UIntReader> data_ptr(data), idx_ptr(idx), starts_ptr(starts);
    XPtr<ULongReader> idx_offsets_ptr(idx_offsets);
    return Rcpp::wrap(XPtr<UIntReader>(new UIntReader(
        std::make_unique<BP128_D1_UIntReader>(
            std::move(*data_ptr), std::move(*idx_ptr), std::move(*idx_offsets_ptr), std::move(*starts_ptr), count
        ),
        1024,
        1024
    )));
}

// [[Rcpp::export]]
SEXP open_bp128_for(SEXP data, SEXP idx, SEXP idx_offsets, uint32_t count) {
    XPtr<UIntReader> data_ptr(data), idx_ptr(idx);
    XPtr<ULongReader> idx_offsets_ptr(idx_offsets);
    return Rcpp::wrap(XPtr<UIntReader>(new UIntReader(
        std::make_unique<BP128_FOR_UIntReader>(std::move(*data_ptr), std::move(*idx_ptr), std::move(*idx_offsets_ptr), count),
        1024,
        1024
    )));
}

// [[Rcpp::export]]
std::string simd_version() {
    return simd::current_target();
}

// [[Rcpp::export]]
std::vector<std::string> available_simd_versions() {
    return simd::supported_targets();
}

// [[Rcpp::export]]
void set_simd_version(std::string version) {
    simd::set_target(version);
}

// [[Rcpp::export]]
std::string simd_version_bp128() {
    return simd::bp128::current_target();
}

// BENCHMARKING FUNCTIONS
// These are mainly for performance benchmarking, as they only go from/to memory,
// and require pre-allocated memory buffers in order to use.
// Also, they don't report how much of the output memory buffer actually gets used which
// is inconvenient for real-world use

void write_unbuffered(const uint32_t *input, size_t n, BulkNumWriter<uint32_t> &to, size_t chunk_size=2048) {
    size_t i = 0;
    std::vector<uint32_t> buf(chunk_size);
    while (i < n) {
        size_t write_size = std::min(chunk_size, n-i);
        memcpy(buf.data(), input + i, sizeof(uint32_t) * write_size);
        i += to.write(buf.data(), write_size);
    }
    to.finalize();
}

void read_unbuffered(BulkNumReader<uint32_t> &from, uint32_t *output, size_t capacity) {
    size_t i = 0;
    while (true) {
        uint64_t loaded = from.load(output + i, capacity-i);
        if (loaded == 0) break;
        if (i == capacity) throw std::runtime_error("read_unbuffered(): Ran out of capacity in write destination");
        i += loaded;
    }
}

// [[Rcpp::export]]
void write_bp128(IntegerVector input, IntegerVector out_data, IntegerVector out_idx, NumericVector out_idx_offsets) {
    BP128UIntWriter writer(
        NumWriter<uint32_t>(std::make_unique<PointerNumWriter<uint32_t>>((uint32_t *) &out_data[0], out_data.size()),1024),
        NumWriter<uint32_t>(std::make_unique<PointerNumWriter<uint32_t>>((uint32_t *) &out_idx[0], out_idx.size()),1024),
        NumWriter<uint64_t>(std::make_unique<PointerNumWriter<uint64_t>>((uint64_t *) &out_idx_offsets[0], out_idx_offsets.size()),1024)
    );
    write_unbuffered((uint32_t *) &input[0], input.size(), writer);
}

// [[Rcpp::export]]
void write_bp128_for(IntegerVector input, IntegerVector out_data, IntegerVector out_idx, NumericVector out_idx_offsets) {
    BP128_FOR_UIntWriter writer(
        NumWriter<uint32_t>(std::make_unique<PointerNumWriter<uint32_t>>((uint32_t *) &out_data[0], out_data.size()),1024),
        NumWriter<uint32_t>(std::make_unique<PointerNumWriter<uint32_t>>((uint32_t *) &out_idx[0], out_idx.size()),1024),
        NumWriter<uint64_t>(std::make_unique<PointerNumWriter<uint64_t>>((uint64_t *) &out_idx_offsets[0], out_idx_offsets.size()),1024)
    );
    write_unbuffered((uint32_t *) &input[0], input.size(), writer);
}

// [[Rcpp::export]]
void write_bp128_d1(IntegerVector input, IntegerVector out_data, IntegerVector out_idx, NumericVector out_idx_offsets, IntegerVector out_starts) {
    BP128_D1_UIntWriter writer(
        NumWriter<uint32_t>(std::make_unique<PointerNumWriter<uint32_t>>((uint32_t *) &out_data[0], out_data.size()),1024),
        NumWriter<uint32_t>(std::make_unique<PointerNumWriter<uint32_t>>((uint32_t *) &out_idx[0], out_idx.size()),1024),
        NumWriter<uint64_t>(std::make_unique<PointerNumWriter<uint64_t>>((uint64_t *) &out_idx_offsets[0], out_idx_offsets.size()),1024),
        NumWriter<uint32_t>(std::make_unique<PointerNumWriter<uint32_t>>((uint32_t *) &out_starts[0], out_starts.size()), 1024)
    );
    write_unbuffered((uint32_t *) &input[0], input.size(), writer);
}

// [[Rcpp::export]]
void write_bp128_d1z(IntegerVector input, IntegerVector out_data, IntegerVector out_idx, NumericVector out_idx_offsets, IntegerVector out_starts) {
    BP128_D1Z_UIntWriter writer(
        NumWriter<uint32_t>(std::make_unique<PointerNumWriter<uint32_t>>((uint32_t *) &out_data[0], out_data.size()),1024),
        NumWriter<uint32_t>(std::make_unique<PointerNumWriter<uint32_t>>((uint32_t *) &out_idx[0], out_idx.size()),1024),
        NumWriter<uint64_t>(std::make_unique<PointerNumWriter<uint64_t>>((uint64_t *) &out_idx_offsets[0], out_idx_offsets.size()),1024),
        NumWriter<uint32_t>(std::make_unique<PointerNumWriter<uint32_t>>((uint32_t *) &out_starts[0], out_starts.size()), 1024)
    );
    write_unbuffered((uint32_t *) &input[0], input.size(), writer);
}

// [[Rcpp::export]]
void read_bp128(IntegerVector input_data, IntegerVector input_idx, NumericVector input_idx_offsets, IntegerVector out, double count) {
    BP128UIntReader reader(
        NumReader<uint32_t>(std::make_unique<VecNumReader<uint32_t>>((uint32_t *) &input_data[0], input_data.size()),1024),
        NumReader<uint32_t>(std::make_unique<VecNumReader<uint32_t>>((uint32_t *) &input_idx[0], input_idx.size()),1024),
        NumReader<uint64_t>(std::make_unique<VecNumReader<uint64_t>>((uint64_t *) &input_idx_offsets[0], input_idx_offsets.size()),1024),
        (uint64_t) count
    );

    read_unbuffered(reader, (uint32_t *) &out[0], out.size());
}

// [[Rcpp::export]]
void read_bp128_for(IntegerVector input_data, IntegerVector input_idx, NumericVector input_idx_offsets, IntegerVector out, double count) {
    BP128_FOR_UIntReader reader(
        NumReader<uint32_t>(std::make_unique<VecNumReader<uint32_t>>((uint32_t *) &input_data[0], input_data.size()),1024),
        NumReader<uint32_t>(std::make_unique<VecNumReader<uint32_t>>((uint32_t *) &input_idx[0], input_idx.size()),1024),
        NumReader<uint64_t>(std::make_unique<VecNumReader<uint64_t>>((uint64_t *) &input_idx_offsets[0], input_idx_offsets.size()),1024),
        (uint64_t) count
    );

    read_unbuffered(reader, (uint32_t *) &out[0], out.size());
}

// [[Rcpp::export]]
void read_bp128_d1(IntegerVector input_data, IntegerVector input_idx, NumericVector input_idx_offsets, IntegerVector input_starts, IntegerVector out, double count) {
    BP128_D1_UIntReader reader(
        NumReader<uint32_t>(std::make_unique<VecNumReader<uint32_t>>((uint32_t *) &input_data[0], input_data.size()),1024),
        NumReader<uint32_t>(std::make_unique<VecNumReader<uint32_t>>((uint32_t *) &input_idx[0], input_idx.size()),1024),
        NumReader<uint64_t>(std::make_unique<VecNumReader<uint64_t>>((uint64_t *) &input_idx_offsets[0], input_idx_offsets.size()),1024),
        NumReader<uint32_t>(std::make_unique<VecNumReader<uint32_t>>((uint32_t *) &input_starts[0], input_starts.size()),1024),
        (uint64_t) count
    );

    read_unbuffered(reader, (uint32_t *) &out[0], out.size());
}

// [[Rcpp::export]]
void read_bp128_d1z(IntegerVector input_data, IntegerVector input_idx, NumericVector input_idx_offsets, IntegerVector input_starts, IntegerVector out, double count) {
    BP128_D1Z_UIntReader reader(
        NumReader<uint32_t>(std::make_unique<VecNumReader<uint32_t>>((uint32_t *) &input_data[0], input_data.size()),1024),
        NumReader<uint32_t>(std::make_unique<VecNumReader<uint32_t>>((uint32_t *) &input_idx[0], input_idx.size()),1024),
        NumReader<uint64_t>(std::make_unique<VecNumReader<uint64_t>>((uint64_t *) &input_idx_offsets[0], input_idx_offsets.size()),1024),
        NumReader<uint32_t>(std::make_unique<VecNumReader<uint32_t>>((uint32_t *) &input_starts[0], input_starts.size()),1024),
        (uint64_t) count
    );

    read_unbuffered(reader, (uint32_t *) &out[0], out.size());
}

// Performance note: This function does need to have a copy to avoid modifying the input `end` buffer
// [[Rcpp::export]]
void write_bp128_end(IntegerVector end, IntegerVector start, IntegerVector out_data, IntegerVector out_idx, NumericVector out_idx_offsets) {
    // Compress the result of end - start, mimicking the algorithm used within BPCells to store fragments
    // (i.e. do simd::sub or simd::add on each load step)

    BP128UIntWriter writer(
        NumWriter<uint32_t>(std::make_unique<PointerNumWriter<uint32_t>>((uint32_t *) &out_data[0], out_data.size()),1024),
        NumWriter<uint32_t>(std::make_unique<PointerNumWriter<uint32_t>>((uint32_t *) &out_idx[0], out_idx.size()),1024),
        NumWriter<uint64_t>(std::make_unique<PointerNumWriter<uint64_t>>((uint64_t *) &out_idx_offsets[0], out_idx_offsets.size()),1024)
    );

    uint32_t *start_data = (uint32_t *) &start[0];
    uint32_t *end_data = (uint32_t *) &end[0];
    size_t capacity = end.size();

    constexpr uint64_t chunk_size = 2048;
    uint32_t buf[chunk_size];
    uint64_t pos = 0;
    while (pos < capacity) {
        size_t write_size = std::min(chunk_size, capacity-pos);
        memcpy(buf, end_data + pos, sizeof(uint32_t) * write_size);

        if (write_size   == 0) break;
        simd::sub(buf, start_data + pos, write_size);

        pos += writer.write(buf, write_size);
    }
    writer.finalize();
}

// [[Rcpp::export]]
void read_bp128_end(IntegerVector input_data, IntegerVector input_idx, NumericVector input_idx_offsets, IntegerVector start, IntegerVector out, double count) {
    // Decompress the result of input + start, mimicking the algorithm used within BPCells to store fragments
    // (i.e. do simd::sub or simd::add on each load step)
    BP128UIntReader reader(
        NumReader<uint32_t>(std::make_unique<VecNumReader<uint32_t>>((uint32_t *) &input_data[0], input_data.size()),1024),
        NumReader<uint32_t>(std::make_unique<VecNumReader<uint32_t>>((uint32_t *) &input_idx[0], input_idx.size()),1024),
        NumReader<uint64_t>(std::make_unique<VecNumReader<uint64_t>>((uint64_t *) &input_idx_offsets[0], input_idx_offsets.size()),1024),
        (uint64_t) count
    );

    uint32_t *start_data = (uint32_t *) &start[0];
    uint32_t *output = (uint32_t *) &out[0];
    size_t capacity = out.size();
    uint64_t chunk_size = 2048;

    size_t i = 0;
    while (true) {
        uint64_t loaded = reader.load(output + i, std::min<uint64_t>(chunk_size, capacity-i));
        if (loaded == 0) break;
        simd::add(output + i, start_data + i, loaded);
        if (i == capacity) throw std::runtime_error("read_unbuffered(): Ran out of capacity in write destination");
        i += loaded;
    }
}