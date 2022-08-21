#include "hdf5.h"

namespace BPCells {


H5StringReader::H5StringReader(const HighFive::Group &group, std::string path) {
    HighFive::SilenceHDF5 s;

    HighFive::DataSet d(group.getDataSet(path));
    HighFive::DataType type = d.getDataType();
    if (type.isVariableStr()) {
        d.read(data);
    } else {
        uint32_t bytes = type.getSize();
        uint32_t elements = d.getDimensions()[0];
        std::vector<char> char_data(bytes * elements);
        d.read(char_data.data(), type);
        data.resize(elements);
        for (uint32_t i = 0; i < elements; i++) {
            data[i] = std::string(char_data.data() + bytes*i, char_data.data() + bytes*(i+1));
        }
    }
}
const char* H5StringReader::get(uint32_t idx) const {
    if (idx < data.size()) return data[idx].c_str();
    return NULL;
}
uint32_t H5StringReader::size() const {return data.size();}

H5StringWriter::H5StringWriter(const HighFive::Group &group, std::string path) :
    group(group), path(path) {}

void H5StringWriter::write(const StringReader &reader) {
    std::vector<std::string> data;
    uint32_t i = 0;
    while (true) {
        const char* s = reader.get(i);
        if (s == NULL) break;
        data.push_back(s);
        i++;
    }
    HighFive::SilenceHDF5 s;
    HighFive::DataSet ds = group.createDataSet<std::string>(path, HighFive::DataSpace::From(data));
    ds.write(data);
}

HighFive::Group createH5Group(std::string file_path, std::string group_path) {
    HighFive::SilenceHDF5 s;
    if (group_path == "") group_path = "/";
    
    std::filesystem::path path(file_path);
    if (path.has_parent_path() && !std::filesystem::exists(path.parent_path())) {
        std::filesystem::create_directories(path.parent_path());
    }

    HighFive::File file(file_path, HighFive::File::OpenOrCreate);
    try {
        HighFive::Group ret(file.getGroup(group_path));
        if (ret.getNumberObjects() != 0) {
            throw std::runtime_error("Requested hdf5 group is not empty");
        }
        return ret;
    } catch (const HighFive::GroupException &e) {
        return file.createGroup(group_path);
    }
}

H5WriterBuilder::H5WriterBuilder(std::string file, std::string group, uint32_t buffer_size, uint32_t chunk_size) :
    group(createH5Group(file, group)), buffer_size(buffer_size), chunk_size(chunk_size) {}

UIntWriter H5WriterBuilder::createUIntWriter(std::string name) {
    return UIntWriter(std::make_unique<H5NumWriter<uint32_t>>(group, name, chunk_size), buffer_size);
}

ULongWriter H5WriterBuilder::createULongWriter(std::string name) {
    return ULongWriter(std::make_unique<H5NumWriter<uint64_t>>(group, name, chunk_size), buffer_size);
}

FloatWriter H5WriterBuilder::createFloatWriter(std::string name) {
    return FloatWriter(std::make_unique<H5NumWriter<float>>(group, name, chunk_size), buffer_size);
}

DoubleWriter H5WriterBuilder::createDoubleWriter(std::string name) {
    return DoubleWriter(std::make_unique<H5NumWriter<double>>(group, name, chunk_size), buffer_size);
}

std::unique_ptr<StringWriter> H5WriterBuilder::createStringWriter(std::string name) {
    return std::make_unique<H5StringWriter>(group, name);
}

void H5WriterBuilder::writeVersion(std::string version) {
    group.createAttribute<std::string>("version", HighFive::DataSpace::From(version)).write(version);
}

void H5WriterBuilder::deleteWriter(std::string name) {
    throw std::logic_error("deleteWriter: HDF5 files don't support deletion");
}

H5ReaderBuilder::H5ReaderBuilder(std::string file, std::string group, uint32_t buffer_size, uint32_t read_size) :
    group(HighFive::File(
            file,
            HighFive::File::ReadWrite
        ).getGroup(group == "" ? std::string("/") : group)),
    buffer_size(buffer_size), read_size(read_size) {}

UIntReader H5ReaderBuilder::openUIntReader(std::string name) {
    return UIntReader(std::make_unique<H5NumReader<uint32_t>>(group, name), buffer_size, read_size);
}

ULongReader H5ReaderBuilder::openULongReader(std::string name) {
    return ULongReader(std::make_unique<H5NumReader<uint64_t>>(group, name), buffer_size, read_size);
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

HighFive::Group& H5ReaderBuilder::getGroup() {return group;}

} // end namespace BPCells