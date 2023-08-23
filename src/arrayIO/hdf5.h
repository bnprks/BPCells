#pragma once


#include "array_interfaces.h"

#include "../lib/highfive/H5DataSet.hpp"
#include "../lib/highfive/H5DataSpace.hpp"
#include "../lib/highfive/H5File.hpp"
#include "../lib/highfive/H5Utility.hpp"

namespace BPCells {

template <class T> class H5NumWriter : public BulkNumWriter<T> {
  private:
    HighFive::DataSet dataset;
    HighFive::DataType datatype = HighFive::create_datatype<T>();

    static HighFive::DataSet
    createH5DataSet(HighFive::Group group, std::string group_path, uint64_t chunk_size) {
        HighFive::SilenceHDF5 s;
        // Create a dataspace with initial shape and max shape
        HighFive::DataSpace dataspace({0}, {HighFive::DataSpace::UNLIMITED});

        // Use chunking
        HighFive::DataSetCreateProps props;
        props.add(HighFive::Chunking(std::vector<hsize_t>{chunk_size}));

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
    H5NumWriter(const HighFive::Group &group, std::string path, uint64_t chunk_size = 1024)
        : dataset(createH5DataSet(group, path, chunk_size)) {}

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
        dataset.select({pos}, {count}).read(out, datatype);
        pos += count;
        return count;
    }
};

using H5UIntReader = H5NumReader<uint32_t>;

class H5StringReader : public StringReader {
  private:
    std::vector<std::string> data;

  public:
    H5StringReader(const HighFive::Group &group, std::string path);
    const char *get(uint64_t idx) const override;
    uint64_t size() const override;
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
    uint64_t buffer_size;
    uint64_t chunk_size;

  public:
    H5WriterBuilder(
        std::string file, std::string group, uint64_t buffer_size = 8192, uint64_t chunk_size = 1024, bool allow_exists = false
    );
    UIntWriter createUIntWriter(std::string name) override;
    ULongWriter createULongWriter(std::string name) override;
    FloatWriter createFloatWriter(std::string name) override;
    DoubleWriter createDoubleWriter(std::string name) override;
    std::unique_ptr<StringWriter> createStringWriter(std::string name) override;
    void writeVersion(std::string version) override;
    void deleteWriter(std::string name) override;
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
    UIntReader openUIntReader(std::string name) override;
    ULongReader openULongReader(std::string name) override;
    FloatReader openFloatReader(std::string name) override;
    DoubleReader openDoubleReader(std::string name) override;
    std::unique_ptr<StringReader> openStringReader(std::string name) override;
    std::string readVersion() override;
    HighFive::Group &getGroup();
};

} // end namespace BPCells