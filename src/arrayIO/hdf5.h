#pragma once

#include <filesystem>

#include "array_interfaces.h"
#include "frags.h"
#include "../matrixIterators/PackedMatrix.h"
#include "../matrixIterators/UnpackedMatrix.h"
#include "../fragmentIterators/UnpackedFragments3.h"
#include "../fragmentIterators/PackedFragments3.h"

#include "../lib/highfive/H5DataSet.hpp"
#include "../lib/highfive/H5DataSpace.hpp"
#include "../lib/highfive/H5File.hpp"
#include "../lib/highfive/H5Utility.hpp"

namespace BPCells {

class MovableDataSet : public HighFive::DataSet {
public:
    MovableDataSet() : HighFive::DataSet(H5I_INVALID_HID) {}
    MovableDataSet(MovableDataSet&& other) = default;
    MovableDataSet(const MovableDataSet&) = delete;
    MovableDataSet(HighFive::DataSet&& other) : HighFive::DataSet(std::move(other)) {}
    MovableDataSet& operator= (MovableDataSet&& other) = default;
    MovableDataSet& operator= (const MovableDataSet&) = delete;
};

class MovableGroup : public HighFive::Group {
public:
    MovableGroup() : HighFive::Group(H5I_INVALID_HID) {}
    MovableGroup(MovableGroup&& other) = default;
    MovableGroup(const MovableGroup&) = delete;
    MovableGroup(HighFive::Group&& other) : HighFive::Group(std::move(other)) {}
    MovableGroup& operator= (MovableGroup&& other) = default;
    MovableGroup& operator= (const MovableGroup&) = delete;
};

class MovableDataType : public HighFive::DataType {
public:
    MovableDataType() : HighFive::DataType(H5I_INVALID_HID) {}
    MovableDataType(MovableDataType&& other) = default;
    MovableDataType(const MovableDataType&) = delete;
    MovableDataType(const HighFive::DataType& other); 
    MovableDataType(HighFive::DataType&& other);
    MovableDataType& operator= (MovableDataType&& other) = default;
    MovableDataType& operator= (const MovableDataType&) = delete;
};

MovableGroup constructH5Group(std::string file_path, std::string group_path);
MovableGroup openH5Group(std::string file_path, std::string group_path);

MovableDataSet constructH5DataSet(std::string file_path, std::string group_path, uint32_t chunk_size);
MovableDataSet constructH5DataSet(HighFive::Group group, std::string group_path, uint32_t chunk_size);

MovableDataSet openH5DataSet(std::string file_path, std::string group_path);
MovableDataSet openH5DataSet(HighFive::Group group, std::string path);

UnpackedMatrix open10xFeatureMatrix(std::string file_path, std::string group_path = "matrix", uint32_t buffer_size = 8192);
std::string loadVersionH5(std::string file_path, std::string group_path);
UnpackedMatrix openUnpackedMatrixH5(std::string file_path, std::string group_path, uint32_t buffer_size = 8192);
UnpackedMatrixWriter createUnpackedMatrixH5(std::string file_path, std::string group_path, uint32_t buffer_size = 8192, uint32_t chunk_size = 1024);
PackedMatrix openPackedMatrixH5(std::string file_path, std::string group_path, uint32_t buffer_size = 8192);
PackedMatrixWriter createPackedMatrixH5(std::string file_path, std::string group_path, uint32_t buffer_size = 8192, uint32_t chunk_size = 1024);

UnpackedFragments3 openUnpackedFragmentsH5(std::string file_path, std::string group_path, uint32_t buffer_size);
UnpackedFragmentsWriter3 createUnpackedFragmentsH5(std::string file_path, std::string group_path, uint32_t buffer_size, uint32_t chunk_size = 1024);

PackedFragments3 openPackedFragmentsH5(std::string file_path, std::string group_path, uint32_t buffer_size);
PackedFragmentsWriter3 createPackedFragmentsH5(std::string file_path, std::string group_path, uint32_t buffer_size, uint32_t chunk_size = 1024);


class ZH5UIntWriter : public UIntWriter {
private:
    MovableDataSet dataset;
    MovableDataType datatype = HighFive::create_datatype<uint32_t>();
    std::vector<uint32_t> buf;
    void flush();
    void _ensureCapacity(size_t capacity) override;
public:
    ~ZH5UIntWriter();
    ZH5UIntWriter() = default;
    ZH5UIntWriter(std::string file_path, std::string group_path, uint32_t buffer_size = 8192, uint32_t chunk_size = 1024);
    ZH5UIntWriter(const HighFive::Group &group, std::string path, uint32_t buffer_size = 8192, uint32_t chunk_size = 1024);
    ZH5UIntWriter(const ZH5UIntWriter&) = delete;
    ZH5UIntWriter& operator= (const ZH5UIntWriter&) = delete;
    ZH5UIntWriter& operator= (ZH5UIntWriter&&) = default;

    void next() override;
};


class ZH5UIntReader : public UIntReader {
private:
    MovableDataSet dataset;
    std::vector<uint32_t> buf;
    size_t pos = 0;
    MovableDataType datatype = HighFive::create_datatype<uint32_t>();
    void _ensureCapacity(size_t capacity) override;
    size_t readPos(uint32_t* out, uint32_t capacity);
public:
    ZH5UIntReader() = default;
    ZH5UIntReader(std::string file_path, std::string group_path, uint32_t buffer_size = 8192);
    ZH5UIntReader(const HighFive::Group &group, std::string path, uint32_t buffer_size = 8192);
    ZH5UIntReader(const ZH5UIntReader&) = delete;
    ZH5UIntReader& operator= (const ZH5UIntReader&) = delete;
    ZH5UIntReader& operator= (ZH5UIntReader&&) = default;

    bool next() override;
    bool seek(size_t pos) override;
};

class H5StringReader : public StringReader {
private:
    std::vector<std::string> data;
public:
    H5StringReader(const MovableDataSet &dataset);
    const char* get(uint32_t idx) const override;
    uint32_t size() const override;
};

class H5StringWriter : public StringWriter {
private:
    std::string file, group;
public:
    H5StringWriter(std::string file, std::string group);
    void write(const StringReader &reader) override;
};


class H5UIntWriter {
private:
    MovableDataSet dataset;
    MovableDataType datatype = HighFive::create_datatype<uint32_t>();
public:
    H5UIntWriter() = default;
    H5UIntWriter(std::string file_path, std::string group_path, uint32_t chunk_size = 250000);
    H5UIntWriter(HighFive::Group group, std::string path, uint32_t chunk_size = 250000);
    H5UIntWriter(const H5UIntWriter&) = delete;
    H5UIntWriter& operator= (const H5UIntWriter&) = delete;
    H5UIntWriter& operator= (H5UIntWriter&&) = default;

    // Append integers to output stream; Throws exception on failure
    void write(const uint32_t *buffer, uint32_t count);
    void finalize();
};

class H5UIntReader {
private:
    MovableDataSet dataset;
    size_t pos = 0;
    MovableDataType datatype = HighFive::create_datatype<uint32_t>();
public:
    H5UIntReader() = default;
    H5UIntReader(std::string file_path, std::string group_path);
    H5UIntReader(HighFive::Group group, std::string path);
    H5UIntReader(const H5UIntReader&) = delete;
    H5UIntReader& operator= (const H5UIntReader&) = delete;
    H5UIntReader& operator= (H5UIntReader&&) = default;

    // Read integers up to count into buffer and return the number actually read. Throw an
    // exception on failures. Must repeatedly return 0 after finishing the stream
    // (and not return 0 before then)
    uint32_t read(uint32_t *buffer, uint32_t count);

    // If this reader is seekable, return the total number of integers stored
    uint32_t size();
    // Seek to the given index in the array
    void seek(const size_t pos);
};

class H5FragmentsSaver {
public:
    using UIntWriter = H5UIntWriter;
    H5FragmentsSaver() = delete;
    H5FragmentsSaver(std::string file_path, std::string group_path, uint32_t chunk_size = 250000);
    H5FragmentsSaver(H5FragmentsSaver&&) = default;
    H5FragmentsSaver(const H5FragmentsSaver&) = delete;
    H5FragmentsSaver& operator= (const H5FragmentsSaver&) = delete;

    UnpackedFrags<UIntWriter> chrWriterUnpacked(uint32_t chr_id);
    PackedFrags<UIntWriter> chrWriterPacked(uint32_t chr_id);

    void writeCellNames(std::vector<std::string> cell_names);
    void writeChrNames(std::vector<std::string> chr_names);
private:
    static MovableGroup openH5Group(std::string file_path, std::string group_path);
    MovableGroup group;
    uint32_t chunk_size;
};

class H5FragmentsLoader {
public:
    using UIntReader = H5UIntReader;
    H5FragmentsLoader() = delete;
    H5FragmentsLoader(std::string file_path, std::string group_path, bool is_packed=false);
    H5FragmentsLoader(H5FragmentsLoader&&) = default;
    H5FragmentsLoader(const H5FragmentsLoader&) = delete;
    H5FragmentsLoader& operator= (const H5FragmentsLoader&) = delete;

    UnpackedFrags<UIntReader> chrReaderUnpacked(uint32_t chr_id);
    PackedFrags<UIntReader> chrReaderPacked(uint32_t chr_id);

    std::vector<std::string> readCellNames();
    std::vector<std::string> readChrNames();

private:
    MovableGroup group;
};

} // end namespace BPCells