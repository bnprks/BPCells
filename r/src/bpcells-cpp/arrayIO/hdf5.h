// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include "array_interfaces.h"
#include "hdf5_threadsafe.h"

namespace BPCells {

template <class T> class H5NumWriter : public BulkNumWriter<T> {
  private:
    H5DataSet1D<T> dataset;

  public:
    H5NumWriter(
        H5Group &group,
        std::string path,
        uint64_t chunk_size = 1024,
        uint32_t gzip_level = 0
    )
        : dataset(group.createDataSet1D<T>(path, 0, chunk_size, gzip_level)) {}

    uint64_t write(T *in, uint64_t count) override {
        dataset.append(in, count);
        return count;
    }
};

using H5UIntWriter = H5NumWriter<uint32_t>;

template <class T> class H5NumReader : public BulkNumReader<T> {
  private:
    H5DataSet1D<T> dataset;
    size_t pos = 0;

  public:
    H5NumReader(const H5Group &group, const std::string &path)
        : dataset(group.openDataSet1D<T>(path)) {}

    // Return total number of integers in the reader
    uint64_t size() const override {
        return dataset.getDimension();
    }

    // Change the next load to start at index pos
    void seek(uint64_t new_pos) override { pos = new_pos; }

    // Copy up to `count` integers into `out`, returning the actual number copied.
    // Will always load >0 unless there is no more input
    uint64_t load(T *out, uint64_t count) override {
        dataset.load(pos, out, count);
        pos += count;
        return count;
    }
};

using H5UIntReader = H5NumReader<uint32_t>;

class H5StringReader : public StringReader {
  private:
    bool data_ready = false;
    H5DataSet1D<std::string> dataset;
    std::vector<std::string> data;

    inline void ensureDataReady();
  public:
    H5StringReader(const H5Group &group, const std::string &path);
    const char *get(uint64_t idx) override;
    uint64_t size() override;
};

class H5StringWriter : public StringWriter {
  private:
    H5Group group;
    std::string path;
    const uint32_t gzip_level;

  public:
    H5StringWriter(const H5Group &group, const std::string &path, uint32_t gzip_level = 0);
    void write(StringReader &reader) override;
};

class H5WriterBuilder final : public WriterBuilder {
  protected:
    H5Group group;
    uint64_t buffer_size;
    uint64_t chunk_size;
    const uint32_t gzip_level;

  public:
    H5WriterBuilder(
        std::string file,
        std::string group,
        uint64_t buffer_size = 8192,
        uint64_t chunk_size = 1024,
        bool allow_exists = false,
        uint32_t gzip_level = 0
    );
    IntWriter createIntWriter(std::string name);
    LongWriter createLongWriter(std::string name);
    UIntWriter createUIntWriter(std::string name) override;
    ULongWriter createULongWriter(std::string name) override;
    FloatWriter createFloatWriter(std::string name) override;
    DoubleWriter createDoubleWriter(std::string name) override;
    template <class T> NumWriter<T> create(std::string name) {
        if constexpr (std::is_same_v<T, int32_t>) {
            return createIntWriter(name);
        }
        if constexpr (std::is_same_v<T, int64_t>) {
            return createLongWriter(name);
        }
        if constexpr (std::is_same_v<T, uint32_t>) {
            return createUIntWriter(name);
        }
        if constexpr (std::is_same_v<T, uint64_t>) {
            return createULongWriter(name);
        }
        if constexpr (std::is_same_v<T, float>) {
            return createFloatWriter(name);
        }
        if constexpr (std::is_same_v<T, double>) {
            return createDoubleWriter(name);
        }
    }

    std::unique_ptr<StringWriter> createStringWriter(std::string name) override;
    void writeVersion(std::string version) override;
    void deleteWriter(std::string name) override;
    H5Group &getGroup();
};

class H5ReaderBuilder final : public ReaderBuilder {
    H5Group group;
    uint64_t buffer_size;
    uint64_t read_size;

  public:
    H5ReaderBuilder(
        std::string file, std::string group, uint64_t buffer_size, uint64_t read_size = 1024
    );
    IntReader openIntReader(std::string name);
    LongReader openLongReader(std::string name);
    UIntReader openUIntReader(std::string name) override;
    ULongReader openULongReader(std::string name) override;
    FloatReader openFloatReader(std::string name) override;
    DoubleReader openDoubleReader(std::string name) override;
    template <class T> NumReader<T> open(std::string name) {
        if constexpr (std::is_same_v<T, int32_t>) {
            return openIntReader(name);
        }
        if constexpr (std::is_same_v<T, int64_t>) {
            return openLongReader(name);
        }
        if constexpr (std::is_same_v<T, uint32_t>) {
            return openUIntReader(name);
        }
        if constexpr (std::is_same_v<T, uint64_t>) {
            return openULongReader(name);
        }
        if constexpr (std::is_same_v<T, float>) {
            return openFloatReader(name);
        }
        if constexpr (std::is_same_v<T, double>) {
            return openDoubleReader(name);
        }
    }
    std::unique_ptr<StringReader> openStringReader(std::string name) override;
    std::string readVersion() override;
    H5Group &getGroup();
};

} // end namespace BPCells