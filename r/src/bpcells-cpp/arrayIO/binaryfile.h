// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <cstring>
#include <fstream>
#include <vector>

#include "../utils/filesystem_compat.h"


#include "array_interfaces.h"

#ifdef _MSC_VER
#include <stdlib.h>
#define bswap_32(x) _byteswap_ulong(x)
#else
#define bswap_32(x) __builtin_bswap32(x)
#endif

namespace BPCells {

template <class T> static std::pair<uint32_t, uint32_t> file_header_magic_number() {
    // UINT32v1 on little endian system
    if constexpr (std::is_same_v<T, uint32_t>) return {0x544e4955, 0x31763233};

    // UINT64v1 on little endian system
    if constexpr (std::is_same_v<T, uint64_t>) return {0x544e4955, 0x31763436};

    // DOUBLEv1 on little endian system
    if constexpr (std::is_same_v<T, double>) return {0x42554F44, 0x3176454C};

    // FLOATSv1 on little endian system
    if constexpr (std::is_same_v<T, float>) return {0x414f4c46, 0x31765354};
}

template <class T> class FileNumWriter final : public BulkNumWriter<T> {
  private:
    std::ofstream file;

  public:
    FileNumWriter(const char *path) {
        // Make sure we get exceptions when things fail
        file.exceptions(std::ofstream::failbit | std::ofstream::badbit);

        // Turn off I/O buffering (Removed because I can't get it to match reasonable performance
        // when I do the buffering manually) file.rdbuf()->pubsetbuf(NULL, 0);

        file.open(path, std::ios_base::binary);
        std::pair<uint32_t, uint32_t> header = file_header_magic_number<T>();
        file.write((char *)&header.first, 4);
        file.write((char *)&header.second, 4);
    }
    uint64_t write(T *in, uint64_t count) override {
        file.write((char *)in, count * sizeof(T));
        return count;
    }
};

using FileUIntWriter = FileNumWriter<uint32_t>;

template <class T> class FileNumReader final : public BulkNumReader<T> {
  protected:
    std::ifstream file;
    uint64_t total_size;

  public:
    FileNumReader(const char *path) {
        file.open(path, std::ios_base::binary);
#if _WIN32
        // Windows-specific: Try to increase maximum open file limit if it is causing
        // the file open operation to fail
        if (!file && _getmaxstdio() == 512) {
            if (_getmaxstdio() == 512) {
                int new_max = 8192;
                while (_setmaxstdio(new_max) == -1 && new_max > 512) {
                    new_max /= 2;
                }
                file.open(path, std::ios_base::binary);
            }
        }
#endif
        if (!file) {
            throw std::runtime_error(
                std::string("Error opening file: ") + strerror(errno) + ": " + path
            );
        }

        uint32_t header[2];
        file.read((char *)header, 8);
        std::pair<uint32_t, uint32_t> magic_number = file_header_magic_number<T>();
        if (header[0] == magic_number.first && header[1] == magic_number.second) {
            // Don't do anything here now, since we removed bytswapping given
            // a lack of accessible testing on big-endian architectures
        } else if (bswap_32(header[0]) == magic_number.first &&
                bswap_32(header[1]) == magic_number.second ) {
            throw std::invalid_argument(
                std::string("Support for big-endian architectures not yet implemented")
            );
        } else {
            throw std::invalid_argument(
                std::string(
                    "File header doesn't match magic number (UINT32v1 or byteswapped TNIU1v23): "
                ) +
                path
            );
        }

        // Detect the file size & cache it
        uint64_t cur = file.tellg();
        file.seekg(0, file.end);
        total_size = (file.tellg() / sizeof(T)) - 8 / sizeof(T);
        file.seekg(cur);
    }

    // Return total number of numbers in the reader
    uint64_t size() const override { return total_size; }

    // Change the next load to start at index pos
    void seek(uint64_t pos) override { file.seekg(8 + pos * sizeof(T)); }

    // Copy up to `count` integers into `out`, returning the actual number copied.
    // Will always load >0 unless there is no more input
    uint64_t load(T *out, uint64_t count) override {
        file.read((char *)out, sizeof(T) * count);
        uint64_t read_count = file.gcount() / sizeof(T);
        return read_count;
    }
};

using FileUIntReader = FileNumReader<uint32_t>;

std::vector<std::string> readLines(std_fs::path path);

class FileStringReader final : public StringReader {
  private:
    bool data_ready = false;
    std_fs::path path;
    std::vector<std::string> data;

    inline void ensureDataReady();
  public:
    FileStringReader(std_fs::path path);
    const char *get(uint64_t idx) override;
    uint64_t size() override;
};

class FileStringWriter final : public StringWriter {
  private:
    std_fs::path path;

  public:
    FileStringWriter(std_fs::path path);
    void write(StringReader &reader) override;
};

class FileWriterBuilder final : public WriterBuilder {
  protected:
    std_fs::path dir;
    uint64_t buffer_size;

  public:
    FileWriterBuilder() = default;
    FileWriterBuilder(std::string dir, uint64_t buffer_size = 8192, bool allow_exists = false);
    FileWriterBuilder &operator=(FileWriterBuilder &&other) = default;
    UIntWriter createUIntWriter(std::string name) override;
    ULongWriter createULongWriter(std::string name) override;
    FloatWriter createFloatWriter(std::string name) override;
    DoubleWriter createDoubleWriter(std::string name) override;
    std::unique_ptr<StringWriter> createStringWriter(std::string name) override;
    void writeVersion(std::string version) override;
    void deleteWriter(std::string name) override;
};

class FileReaderBuilder final : public ReaderBuilder {
    std_fs::path dir;
    uint64_t buffer_size;
    uint64_t read_size;

  public:
    FileReaderBuilder() = default;
    FileReaderBuilder(std::string dir, uint64_t buffer_size = 8192, uint64_t read_size = 1024);
    FileReaderBuilder &operator=(FileReaderBuilder &&other) = default;
    UIntReader openUIntReader(std::string name) override;
    ULongReader openULongReader(std::string name) override;
    FloatReader openFloatReader(std::string name) override;
    DoubleReader openDoubleReader(std::string name) override;
    std::unique_ptr<StringReader> openStringReader(std::string name) override;
    std::string readVersion() override;
};

} // end namespace BPCells