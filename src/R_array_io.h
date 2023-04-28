#pragma once

#include "arrayIO/array_interfaces.h"
#include "arrayIO/vector.h"
#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <Rcpp.h>

// These are interfaces for making R objects play nice with the arrayIO interfaces
// used for storing matrices & fragments
// Note that we do some bit reinterpretation: floats represented as ints
// and longs represented as doubles

class RcppStringReader : public BPCells::StringReader {
  private:
    const Rcpp::StringVector data;

  public:
    RcppStringReader(const Rcpp::StringVector &data);
    const char *get(uint64_t idx) const override;
    uint64_t size() const override;
};

class S4ReaderBuilder : public BPCells::ReaderBuilder {
  private:
    Rcpp::S4 s4;
    uint32_t load_size;

  public:
    S4ReaderBuilder(Rcpp::S4 s4, uint32_t load_size = 1024);
    BPCells::UIntReader openUIntReader(std::string name) override;
    BPCells::ULongReader openULongReader(std::string name) override;
    BPCells::FloatReader openFloatReader(std::string name) override;
    BPCells::DoubleReader openDoubleReader(std::string name) override;
    std::unique_ptr<BPCells::StringReader> openStringReader(std::string name) override;
    std::string readVersion() override;
};

class ListWriterBuilder final : public BPCells::VecReaderWriterBuilder {
  public:
    ListWriterBuilder(uint32_t chunk_size = 1024);

    Rcpp::List getList();
};