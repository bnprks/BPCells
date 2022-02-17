#pragma once

#include <Rcpp.h>
#include "arrayIO/array_interfaces.h"

using namespace Rcpp;
using namespace BPCells;

// These are interfaces for making R objects play nice with the arrayIO interfaces
// used for storing matrices & fragments

class RcppStringReader : public StringReader {
private:
    StringVector data;
public:
    inline RcppStringReader(StringVector data) : data(data) {}
    inline const char* get(uint32_t idx) const override {
        if (idx < data.size()) return data[idx];
        return NULL;
    }
    inline uint32_t size() const override {return data.size();} 
};

class S4ReaderBuilder : public ReaderBuilder {
private:
    S4 s4;
    uint32_t load_size;
public:
    inline S4ReaderBuilder(S4 s4, uint32_t load_size=1024) : s4(s4), load_size(load_size) {}
    inline UIntReader openUIntReader(std::string name) override {
        IntegerVector v = s4.slot(name);
        return UIntReader(
            std::make_unique<VecUIntReader>((uint32_t *) &v[0], v.size()),
            load_size
        );
    }
    inline std::unique_ptr<StringReader> openStringReader(std::string name) override {
        return std::make_unique<RcppStringReader>(s4.slot(name));
    }
    inline std::string readVersion() override {
        return s4.slot("version");
    }
};

class ListWriterBuilder final : public WriterBuilder {
protected:
    std::map<std::string, std::vector<uint32_t>> int_vecs;
    std::map<std::string, std::vector<std::string>> string_vecs;
    std::string version;
    uint32_t chunk_size;
public:
    ListWriterBuilder(uint32_t chunk_size = 1024) : chunk_size(chunk_size) {}

    UIntWriter createUIntWriter(std::string name) override {
        int_vecs[name] = std::vector<uint32_t>();
        return UIntWriter(
            std::make_unique<VecUIntWriter>(int_vecs.at(name)),
            chunk_size
        );
    }
    std::unique_ptr<StringWriter> createStringWriter(std::string name) override {
        string_vecs.emplace(name, std::vector<std::string>());
        return std::make_unique<VecStringWriter>(string_vecs.at(name));
    }
    void writeVersion(std::string v) override {version = v;}
    
    List getList() {
        List l;
        // Clear out the maps while constructing the list
        for(auto it = int_vecs.begin(); it != int_vecs.end(); it = int_vecs.erase(it)) {
            l[it->first] = IntegerVector((int *) it->second.data(), (int *) it->second.data() + it->second.size());
        }

        for(auto it = string_vecs.begin(); it != string_vecs.end(); it = string_vecs.erase(it)) {
            l[it->first] = StringVector::import(it->second.begin(), it->second.end());
        }

        for (auto const &v: string_vecs) {
            l[v.first] = StringVector::import(v.second.begin(), v.second.end());
        }
        l["version"] = version;
        return l;
    }
};