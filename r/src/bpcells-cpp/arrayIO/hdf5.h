// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include "array_interfaces.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Utility.hpp>

namespace BPCells {

template <class T> class H5NumWriter : public BulkNumWriter<T> {
  private:
    HighFive::DataSet dataset;
    HighFive::DataType datatype = HighFive::create_datatype<T>();

    static HighFive::DataSet createH5DataSet(
        HighFive::Group group, std::string group_path, uint64_t chunk_size, uint32_t gzip_level
    ) {
        HighFive::SilenceHDF5 s;
        // Create a dataspace with initial shape and max shape
        HighFive::DataSpace dataspace({0}, {HighFive::DataSpace::UNLIMITED});

        // Use chunking
        HighFive::DataSetCreateProps props;
        props.add(HighFive::Chunking(std::vector<hsize_t>{chunk_size}));
        if (gzip_level > 0) {
            props.add(HighFive::Shuffle());
            props.add(HighFive::Deflate(gzip_level));
        }

        // At one point I considered using more aggressive chunk caching, but I
        // don't think it's necessary anymore
        // HighFive::DataSetAccessProps a_props;
        // a_props.add(HighFive::Caching(521, 50<<20));// 50MB cache for overkill

        if (group.exist(group_path)) {
            group.unlink(group_path);
        }

        // Create the dataset
        return group.createDataSet<T>(group_path, dataspace, props);
    }

  public:
    H5NumWriter(
        const HighFive::Group &group,
        std::string path,
        uint64_t chunk_size = 1024,
        uint32_t gzip_level = 0
    )
        : dataset(createH5DataSet(group, path, chunk_size, gzip_level)) {}

    uint64_t write(T *in, uint64_t count) override {
        uint64_t cur_size = dataset.getDimensions()[0];
        dataset.resize({cur_size + count});
        dataset.select({cur_size}, {count}).write_raw(in, datatype);
        return count;
    }
};

using H5UIntWriter = H5NumWriter<uint32_t>;

template <class T> class H5NumReader : public BulkNumReader<T> {
  private:
    HighFive::DataSet dataset;
    size_t pos = 0;
    HighFive::DataType datatype = HighFive::create_datatype<T>();

  public:
    H5NumReader(const HighFive::Group &group, std::string path) : dataset(group.getDataSet(path)) {}

    // Return total number of integers in the reader
    uint64_t size() const override { return dataset.getDimensions()[0]; }

    // Change the next load to start at index pos
    void seek(uint64_t new_pos) override { pos = new_pos; }

    // Copy up to `count` integers into `out`, returning the actual number copied.
    // Will always load >0 unless there is no more input
    uint64_t load(T *out, uint64_t count) override {
        dataset.select({pos}, {count}).read_raw(out, datatype);
        pos += count;
        return count;
    }
};

using H5UIntReader = H5NumReader<uint32_t>;

class H5StringReader : public StringReader {
  private:
    bool data_ready = false;
    HighFive::DataSet dataset;
    std::vector<std::string> data;

    inline void ensureDataReady();
  public:
    H5StringReader(const HighFive::Group &group, std::string path);
    const char *get(uint64_t idx) override;
    uint64_t size() override;
};

class H5StringWriter : public StringWriter {
  private:
    HighFive::Group group;
    std::string path;
    const uint32_t gzip_level;

  public:
    H5StringWriter(const HighFive::Group &group, std::string path, uint32_t gzip_level = 0);
    void write(StringReader &reader) override;
};

class H5WriterBuilder final : public WriterBuilder {
  protected:
    HighFive::Group group;
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
    HighFive::Group &getGroup();
};

// Try to open a file for read-write, then fall back to read only if needed.
// If we first open a file ReadOnly, it prevents future opening with ReadWrite
// (bad if we want to read + write the same file).
// This retry makes it possible to still open a file if it's read-only though.
HighFive::File openH5ForReading(const std::string &path);

class H5ReaderBuilder final : public ReaderBuilder {
    HighFive::Group group;
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
    HighFive::Group &getGroup();
};

} // end namespace BPCells