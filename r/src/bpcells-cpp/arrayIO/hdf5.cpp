// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "hdf5.h"
#include "../utils/filesystem_compat.h"
#include <mutex>

namespace BPCells {

std::recursive_mutex &bpcells_hdf5_global_lock() {
    static std::recursive_mutex mut;
    return mut;
}

H5StringReader::H5StringReader(H5Group &group, const std::string &path)
    : dataset(group.openDataSet1D(path)) {}

inline void H5StringReader::ensureDataReady() {
    if (data_ready) return;
    dataset.load(0, data, dataset.getDimension());
    data_ready = true;
}
const char *H5StringReader::get(uint64_t idx) {
    ensureDataReady();
    if (idx < data.size()) return data[idx].c_str();
    return NULL;
}
uint64_t H5StringReader::size() {
    return dataset.getDimension();
}

H5StringWriter::H5StringWriter(const H5Group &group, const std::string &path, uint32_t gzip_level)
    : group(group)
    , path(path)
    , gzip_level(gzip_level) {}

void H5StringWriter::write(StringReader &reader) {
    std::vector<std::string> data;
    uint64_t i = 0;
    while (true) {
        const char *s = reader.get(i);
        if (s == NULL) break;
        data.push_back(s);
        i++;
    }

    group.createDataSet1D<std::string>(path, data.size(), 0, gzip_level).write(data);
}

H5WriterBuilder::H5WriterBuilder(
    std::string file,
    std::string group,
    uint64_t buffer_size,
    uint64_t chunk_size,
    bool allow_exists,
    uint32_t gzip_level
)
    : group(file, group, allow_exists ? H5Group::OpenMode::WriteOrCreate : H5Group::OpenMode::Create)
    , buffer_size(buffer_size)
    , chunk_size(chunk_size)
    , gzip_level(gzip_level) {}

IntWriter H5WriterBuilder::createIntWriter(std::string name) {
    return IntWriter(
        std::make_unique<H5NumWriter<int32_t>>(group, name, chunk_size, gzip_level), buffer_size
    );
}

LongWriter H5WriterBuilder::createLongWriter(std::string name) {
    return LongWriter(
        std::make_unique<H5NumWriter<int64_t>>(group, name, chunk_size, gzip_level), buffer_size
    );
}

UIntWriter H5WriterBuilder::createUIntWriter(std::string name) {
    return UIntWriter(
        std::make_unique<H5NumWriter<uint32_t>>(group, name, chunk_size, gzip_level), buffer_size
    );
}

ULongWriter H5WriterBuilder::createULongWriter(std::string name) {
    return ULongWriter(
        std::make_unique<H5NumWriter<uint64_t>>(group, name, chunk_size, gzip_level), buffer_size
    );
}

FloatWriter H5WriterBuilder::createFloatWriter(std::string name) {
    return FloatWriter(
        std::make_unique<H5NumWriter<float>>(group, name, chunk_size, gzip_level), buffer_size
    );
}

DoubleWriter H5WriterBuilder::createDoubleWriter(std::string name) {
    return DoubleWriter(
        std::make_unique<H5NumWriter<double>>(group, name, chunk_size, gzip_level), buffer_size
    );
}

std::unique_ptr<StringWriter> H5WriterBuilder::createStringWriter(std::string name) {
    return std::make_unique<H5StringWriter>(group, name, gzip_level);
}

void H5WriterBuilder::writeVersion(std::string version) {
    group.setAttribute("version", version);
}

void H5WriterBuilder::deleteWriter(std::string name) {
    throw std::logic_error("deleteWriter: HDF5 files don't support deletion");
}

H5Group &H5WriterBuilder::getGroup() { return group; }


H5ReaderBuilder::H5ReaderBuilder(
    std::string file, std::string group, uint64_t buffer_size, uint64_t read_size
)
    : group(file, group, H5Group::OpenMode::Read)
    , buffer_size(buffer_size)
    , read_size(read_size) {}

IntReader H5ReaderBuilder::openIntReader(std::string name) {
    return IntReader(std::make_unique<H5NumReader<int32_t>>(group, name), buffer_size, read_size);
}

LongReader H5ReaderBuilder::openLongReader(std::string name) {
    return LongReader(
        std::make_unique<H5NumReader<int64_t>>(group, name), buffer_size, read_size
    );
}

UIntReader H5ReaderBuilder::openUIntReader(std::string name) {
    return UIntReader(std::make_unique<H5NumReader<uint32_t>>(group, name), buffer_size, read_size);
}

ULongReader H5ReaderBuilder::openULongReader(std::string name) {
    return ULongReader(
        std::make_unique<H5NumReader<uint64_t>>(group, name), buffer_size, read_size
    );
}

FloatReader H5ReaderBuilder::openFloatReader(std::string name) {
    return FloatReader(std::make_unique<H5NumReader<float>>(group, name), buffer_size, read_size);
}

DoubleReader H5ReaderBuilder::openDoubleReader(std::string name) {
    return DoubleReader(std::make_unique<H5NumReader<double>>(group, name), buffer_size, read_size);
}

std::unique_ptr<StringReader> H5ReaderBuilder::openStringReader(std::string name) {
    return std::make_unique<H5StringReader>(group, name);
}

std::string H5ReaderBuilder::readVersion() {
    return group.getAttribute("version");
}

H5Group &H5ReaderBuilder::getGroup() { return group; }

} // end namespace BPCells
