#pragma once

#include "arrayIO/array_interfaces.h"
#include "arrayIO/vector.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace BPCells;

// These are interfaces for making R objects play nice with the arrayIO interfaces
// used for storing matrices & fragments
// Note that we do some bit reinterpretation: floats represented as ints
// and longs represented as doubles

class RcppStringReader : public StringReader {
  private:
    const StringVector data;

  public:
    inline RcppStringReader(const StringVector &data) : data(data) {}
    inline const char *get(uint32_t idx) const override {
        if (idx < data.size()) return data[idx];
        return NULL;
    }
    inline uint32_t size() const override { return data.size(); }
};

class S4ReaderBuilder : public ReaderBuilder {
  private:
    S4 s4;
    uint32_t load_size;

  public:
    inline S4ReaderBuilder(S4 s4, uint32_t load_size = 1024) : s4(s4), load_size(load_size) {}
    inline UIntReader openUIntReader(std::string name) override {
        IntegerVector v = s4.slot(name);
        return UIntReader(std::make_unique<VecUIntReader>((uint32_t *)&v[0], v.size()), load_size);
    }
    inline ULongReader openULongReader(std::string name) override {
        throw std::logic_error("S4Reader doesn't implement ULong type");
    }
    inline FloatReader openFloatReader(std::string name) override {
        IntegerVector v = s4.slot(name);
        return FloatReader(
            std::make_unique<VecNumReader<float>>((float *)&v[0], v.size()), load_size
        );
    }
    inline DoubleReader openDoubleReader(std::string name) override {
        NumericVector v = s4.slot(name);
        return DoubleReader(
            std::make_unique<VecNumReader<double>>((double *)&v[0], v.size()), load_size
        );
    }
    inline std::unique_ptr<StringReader> openStringReader(std::string name) override {
        return std::make_unique<RcppStringReader>(s4.slot(name));
    }
    inline std::string readVersion() override { return s4.slot("version"); }
};

class ListWriterBuilder final : public VecReaderWriterBuilder {
  public:
    ListWriterBuilder(uint32_t chunk_size = 1024) : VecReaderWriterBuilder(chunk_size) {}

    List getList() {
        List l;
        // Clear out the maps while constructing the list
        for (auto it = int_vecs.begin(); it != int_vecs.end(); it = int_vecs.erase(it)) {
            l[it->first] = IntegerVector(
                (int *)it->second.data(), (int *)it->second.data() + it->second.size()
            );
        }

        // Floats
        for (auto it = float_vecs.begin(); it != float_vecs.end(); it = float_vecs.erase(it)) {
            l[it->first] = IntegerVector(
                (int *)it->second.data(), (int *)it->second.data() + it->second.size()
            );
        }

        // Longs
        for (auto it = long_vecs.begin(); it != long_vecs.end(); it = long_vecs.erase(it)) {
            l[it->first] = NumericVector(
                (double *)it->second.data(), (double *)it->second.data() + it->second.size()
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
};