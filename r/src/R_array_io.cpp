// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "R_array_io.h"

#include "bpcells-cpp/arrayIO/array_interfaces.h"
#include "bpcells-cpp/arrayIO/vector.h"
#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <Rcpp.h>

using namespace Rcpp;
using namespace BPCells;


// [[Rcpp::export]]
NumericVector convert_ulong_to_numeric(const NumericVector &ulong_vec) {
    NumericVector numeric_vec(ulong_vec.size());

    for (int64_t i = 0; i < ulong_vec.size(); i++) {
        uint64_t val = *((uint64_t *) &ulong_vec[i]);
        numeric_vec[i] = (double) val;
    }
    return numeric_vec;
}

// These are interfaces for making R objects play nice with the arrayIO interfaces
// used for storing matrices & fragments
// Note that we do some bit reinterpretation: floats represented as ints
// and longs represented as doubles

RcppStringReader::RcppStringReader(const StringVector &data) : data(data) {}
const char *RcppStringReader::get(uint64_t idx) {
    if ((int64_t)idx < data.size()) return data[idx];
    return NULL;
}
uint64_t RcppStringReader::size() { return data.size(); }

S4ReaderBuilder::S4ReaderBuilder(S4 s4, uint32_t load_size) : s4(s4), load_size(load_size) {}
UIntReader S4ReaderBuilder::openUIntReader(std::string name) {
    IntegerVector v = s4.slot(name);
    return UIntReader(std::make_unique<VecUIntReader>((uint32_t *)&v[0], v.size()), load_size);
}
ULongReader S4ReaderBuilder::openULongReader(std::string name) {
    NumericVector v = s4.slot(name);
    return DoubleReader(
        std::make_unique<VecNumReader<double>>((double *)&v[0], v.size()), load_size
    ).convert<uint64_t>();
}
FloatReader S4ReaderBuilder::openFloatReader(std::string name) {
    IntegerVector v = s4.slot(name);
    return FloatReader(std::make_unique<VecNumReader<float>>((float *)&v[0], v.size()), load_size);
}
DoubleReader S4ReaderBuilder::openDoubleReader(std::string name) {
    NumericVector v = s4.slot(name);
    return DoubleReader(
        std::make_unique<VecNumReader<double>>((double *)&v[0], v.size()), load_size
    );
}
std::unique_ptr<StringReader> S4ReaderBuilder::openStringReader(std::string name) {
    return std::make_unique<RcppStringReader>(Rcpp::as<StringVector>(s4.slot(name)));
}
std::string S4ReaderBuilder::readVersion() { return s4.slot("version"); }

ListWriterBuilder::ListWriterBuilder(uint32_t chunk_size)
    : VecReaderWriterBuilder(chunk_size) {}

List ListWriterBuilder::getList() {
    List l;
    // Clear out the maps while constructing the list
    for (auto it = int_vecs.begin(); it != int_vecs.end(); it = int_vecs.erase(it)) {
        l[it->first] =
            IntegerVector((int *)it->second.data(), (int *)it->second.data() + it->second.size());
    }

    // Floats
    for (auto it = float_vecs.begin(); it != float_vecs.end(); it = float_vecs.erase(it)) {
        l[it->first] =
            IntegerVector((int *)it->second.data(), (int *)it->second.data() + it->second.size());
    }

    // Longs
    for (auto it = long_vecs.begin(); it != long_vecs.end(); it = long_vecs.erase(it)) {
        std::vector<double> vec_converted(it->second.size());
        for (uint64_t i = 0; i < vec_converted.size(); i++) {
            vec_converted[i] = (double) it->second[i];
        }
        l[it->first] = NumericVector(
            (double *)vec_converted.data(), (double *)vec_converted.data() + vec_converted.size()
        );
    }

    // Doubles
    for (auto it = double_vecs.begin(); it != double_vecs.end(); it = double_vecs.erase(it)) {
        l[it->first] = NumericVector(
            (double *)it->second.data(), (double *)it->second.data() + it->second.size()
        );
    }

    for (auto it = string_vecs.begin(); it != string_vecs.end(); it = string_vecs.erase(it)) {
        l[it->first] = StringVector::import(it->second.begin(), it->second.end());
    }

    for (auto const &v : string_vecs) {
        l[v.first] = StringVector::import(v.second.begin(), v.second.end());
    }
    l["version"] = version;
    return l;
}
