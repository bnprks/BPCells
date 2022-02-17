#pragma once

#include <filesystem>

#include "array_interfaces.h"

#include "../lib/highfive/H5DataSet.hpp"
#include "../lib/highfive/H5DataSpace.hpp"
#include "../lib/highfive/H5File.hpp"
#include "../lib/highfive/H5Utility.hpp"

namespace BPCells {

class H5UIntWriter : public UIntBulkWriter {
private:
    HighFive::DataSet dataset;
    HighFive::DataType datatype = HighFive::create_datatype<uint32_t>();

    static HighFive::DataSet createH5DataSet(HighFive::Group group, std::string group_path, uint32_t chunk_size);
public:
    H5UIntWriter(const HighFive::Group &group, std::string path, uint32_t chunk_size = 1024);
    
    uint32_t write(uint32_t *in, uint32_t count) override;
};


class H5UIntReader : public UIntBulkReader {
private:
    HighFive::DataSet dataset;
    size_t pos = 0;
    HighFive::DataType datatype = HighFive::create_datatype<uint32_t>();
public:
    H5UIntReader(const HighFive::Group &group, std::string path);
    
    // Return total number of integers in the reader
    uint32_t size() const override;

    // Change the next load to start at index pos
    void seek(uint32_t pos) override;

    // Copy up to `count` integers into `out`, returning the actual number copied.
    // Will always load >0 unless there is no more input
    uint32_t load(uint32_t *out, uint32_t count) override;
};

class H5StringReader : public StringReader {
private:
    std::vector<std::string> data;
public:
    H5StringReader(const HighFive::Group &group, std::string path);
    const char* get(uint32_t idx) const override;
    uint32_t size() const override;
};

class H5StringWriter : public StringWriter {
private:
    HighFive::Group group;
    std::string path;
public:
    H5StringWriter(const HighFive::Group &group, std::string path);
    void write(const StringReader &reader) override;
};

class H5WriterBuilder final : public WriterBuilder {
protected:
    HighFive::Group group;
    uint32_t buffer_size;
    uint32_t chunk_size;
public:
    H5WriterBuilder(std::string file, std::string group, uint32_t buffer_size = 8192, uint32_t chunk_size = 1024);
    UIntWriter createUIntWriter(std::string name) override;
    std::unique_ptr<StringWriter> createStringWriter(std::string name) override;
    void writeVersion(std::string version) override;
};

class H5ReaderBuilder final : public ReaderBuilder {
    HighFive::Group group;
    uint32_t buffer_size;
    uint32_t read_size;
public:
    H5ReaderBuilder(std::string file, std::string group, uint32_t buffer_size, uint32_t read_size=1024);
    UIntReader openUIntReader(std::string name) override;
    std::unique_ptr<StringReader> openStringReader(std::string name) override;
    std::string readVersion() override;
};

} // end namespace BPCells