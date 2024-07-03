// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "hdf5.h"
#include "../utils/filesystem_compat.h"

namespace BPCells {

H5StringReader::H5StringReader(const HighFive::Group &group, std::string path) : dataset(group.getDataSet(path)) {}

inline void H5StringReader::ensureDataReady() {
    if (data_ready) return;
    HighFive::DataType type = dataset.getDataType();
    if (type.isVariableStr()) {
        // Workaround for HighFive bug: don't try reading an empty string vector
        if (dataset.getDimensions()[0] > 0) {
            dataset.read(data);
        } else {
            data.resize(0);
        }
    } else {
        uint64_t bytes = type.getSize();
        uint64_t elements = dataset.getDimensions()[0];
        std::vector<char> char_data(bytes * elements);
        dataset.read_raw(char_data.data(), type);
        data.resize(elements);
        for (uint64_t i = 0; i < elements; i++) {
            data[i] = std::string(char_data.data() + bytes * i, char_data.data() + bytes * (i + 1));
        }
    }
    data_ready = true;
}
const char *H5StringReader::get(uint64_t idx) {
    ensureDataReady();
    if (idx < data.size()) return data[idx].c_str();
    return NULL;
}
uint64_t H5StringReader::size() { 
    return dataset.getDimensions()[0]; 
}

H5StringWriter::H5StringWriter(const HighFive::Group &group, std::string path, uint32_t gzip_level)
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
    HighFive::SilenceHDF5 s;

    HighFive::DataSetCreateProps props;
    if (gzip_level > 0) {
        props.add(HighFive::Deflate(gzip_level));
    }

    if (group.exist(path)) {
        group.unlink(path);
    }
    HighFive::DataSet ds =
        group.createDataSet<std::string>(path, HighFive::DataSpace::From(data), props);
    // Safety check to avoid an ASan complaint in the R test "AnnData and 10x row/col rename works"
    if (data.size() > 0) {
        ds.write(data);
    }
}

HighFive::Group createH5Group(std::string file_path, std::string group_path, bool allow_exists) {
    HighFive::SilenceHDF5 s;
    if (group_path == "") group_path = "/";

    std_fs::path path(file_path);
    if (path.has_parent_path() && !std_fs::exists(path.parent_path())) {
        std_fs::create_directories(path.parent_path());
    }

    HighFive::File file(file_path, HighFive::File::OpenOrCreate);
    try {
        HighFive::Group ret(file.getGroup(group_path));
        if (!allow_exists && ret.getNumberObjects() != 0) {
            throw std::runtime_error("Requested hdf5 group is not empty");
        }
        return ret;
    } catch (const HighFive::GroupException &e) {
        return file.createGroup(group_path);
    }
}

H5WriterBuilder::H5WriterBuilder(
    std::string file,
    std::string group,
    uint64_t buffer_size,
    uint64_t chunk_size,
    bool allow_exists,
    uint32_t gzip_level
)
    : group(createH5Group(file, group, allow_exists))
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
    if (group.hasAttribute("version")) {
        group.getAttribute("version").write(version);
    } else {
        group.createAttribute<std::string>("version", HighFive::DataSpace::From(version))
            .write(version);
    }
}

void H5WriterBuilder::deleteWriter(std::string name) {
    throw std::logic_error("deleteWriter: HDF5 files don't support deletion");
}

HighFive::Group &H5WriterBuilder::getGroup() { return group; }

// Try to open a file for read-write, then fall back to read only if needed.
// If we first open a file ReadOnly, it prevents future opening with ReadWrite
// (bad if we want to read + write the same file).
// This retry makes it possible to still open a file if it's read-only though.
HighFive::File openH5ForReading(const std::string &path) {
    try {
        HighFive::SilenceHDF5 s;
        return HighFive::File(path, HighFive::File::ReadWrite);
    } catch (const HighFive::FileException &f) {
        return HighFive::File(path, HighFive::File::ReadOnly);
    }
}

H5ReaderBuilder::H5ReaderBuilder(
    std::string file, std::string group, uint64_t buffer_size, uint64_t read_size
)
    : group(openH5ForReading(file).getGroup(group == "" ? std::string("/") : group))
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
    std::string version;
    group.getAttribute("version").read(version);
    return version;
}

HighFive::Group &H5ReaderBuilder::getGroup() { return group; }

} // end namespace BPCells
