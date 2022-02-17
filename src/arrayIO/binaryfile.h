#pragma once

#include <cstring>
#include <fstream>
#include <vector>
#include <filesystem>

#include "array_interfaces.h"

namespace BPCells {

class FileUIntWriter final : public UIntBulkWriter {
private:
    std::ofstream file;
public:
    FileUIntWriter(const char* path);
    uint32_t write(uint32_t *in, uint32_t count) override;
};


class FileUIntReader final : public UIntBulkReader {
protected:
    std::ifstream file;
    uint32_t total_size;
    bool byte_swap;
public:
    FileUIntReader(const char* path);
    
    // Return total number of integers in the reader
    uint32_t size() const override;

    // Change the next load to start at index pos
    void seek(uint32_t pos) override;

    // Copy up to `count` integers into `out`, returning the actual number copied.
    // Will always load >0 unless there is no more input
    uint32_t load(uint32_t *out, uint32_t count) override;
};

std::vector<std::string> readLines(std::filesystem::path path);

class FileStringReader final : public StringReader {
private:
    std::vector<std::string> data;
public:
    FileStringReader(std::filesystem::path path);
    const char* get(uint32_t idx) const override;
    uint32_t size() const override;
};

class FileStringWriter final : public StringWriter {
private:
    std::filesystem::path path;
public:
    FileStringWriter(std::filesystem::path path);
    void write(const StringReader &reader) override;
};

class FileWriterBuilder final : public WriterBuilder {
protected:
    std::filesystem::path dir;
    uint32_t buffer_size;
public:
    FileWriterBuilder(std::string dir, uint32_t buffer_size);
    UIntWriter createUIntWriter(std::string name) override;
    std::unique_ptr<StringWriter> createStringWriter(std::string name) override;
    void writeVersion(std::string version) override;
};

class FileReaderBuilder final : public ReaderBuilder {
    std::filesystem::path dir;
    uint32_t buffer_size;
    uint32_t read_size;
public:
    FileReaderBuilder(std::string dir, uint32_t buffer_size, uint32_t read_size=1024);
    UIntReader openUIntReader(std::string name) override;
    std::unique_ptr<StringReader> openStringReader(std::string name) override;
    std::string readVersion() override;
};


} // end namespace BPCells