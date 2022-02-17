#pragma once

#include <stdint.h>
#include <map>
#include <vector>
#include <string>

#include "array_interfaces.h"

namespace BPCells {
    

class VecUIntWriter final : public UIntBulkWriter {
private:
    std::vector<uint32_t> &vec;
public:
    VecUIntWriter(std::vector<uint32_t> &vec);
    
    uint32_t write(uint32_t *in, uint32_t count) override;
};

class VecUIntReader : public UIntBulkReader {
private:
    const uint32_t *vec;
    uint32_t capacity;
    uint32_t pos = 0;
public:
    VecUIntReader(const uint32_t *vec, std::size_t capacity);

    // Return total number of integers in the reader
    uint32_t size() const override;

    // Change the next load to start at index pos
    void seek(uint32_t pos) override;

    // Copy up to `count` integers into `out`, returning the actual number copied.
    // Will always load >0 unless there is no more input
    uint32_t load(uint32_t *out, uint32_t count) override;
};

class VecStringWriter : public StringWriter {
private:
    std::vector<std::string> &data;
public:
    VecStringWriter(std::vector<std::string> &data);
    void write(const StringReader &reader) override;
};

class VecReaderWriterBuilder final : public WriterBuilder, public ReaderBuilder {
protected:
    std::map<std::string, std::vector<uint32_t>> int_vecs;
    std::map<std::string, std::vector<std::string>> string_vecs;
    std::string version;
    uint32_t chunk_size;
public:
    VecReaderWriterBuilder(uint32_t chunk_size = 1024);
    UIntWriter createUIntWriter(std::string name) override;
    std::unique_ptr<StringWriter> createStringWriter(std::string name) override;
    void writeVersion(std::string version) override;
    UIntReader openUIntReader(std::string name) override;
    std::unique_ptr<StringReader> openStringReader(std::string name) override;
    std::string readVersion() override;

    std::vector<uint32_t>& getIntVec(std::string name);
    std::vector<std::string>& getStringVec(std::string name);
};

} // end namespace BPCells