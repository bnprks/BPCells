#include "hdf5.h"

namespace BPCells {


H5UIntWriter::H5UIntWriter(const HighFive::Group &group, std::string path, uint32_t chunk_size) : 
    dataset(createH5DataSet(group, path, chunk_size)) {}

HighFive::DataSet H5UIntWriter::createH5DataSet(HighFive::Group group, std::string group_path, uint32_t chunk_size) {
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

    // Create the dataset
    return group.createDataSet<uint32_t>(group_path, dataspace, props);
}

uint32_t H5UIntWriter::write(uint32_t *in, uint32_t count) {
    uint32_t cur_size = dataset.getDimensions()[0];
    dataset.resize({cur_size + count});
    dataset.select({cur_size}, {count}).write_raw(in, datatype);
    return count;
}

H5UIntReader::H5UIntReader(const HighFive::Group &group, std::string path) :
    dataset(group.getDataSet(path)) {}

uint32_t H5UIntReader::size() const {return dataset.getDimensions()[0];}

void H5UIntReader::seek(uint32_t new_pos) {pos = new_pos;}

uint32_t H5UIntReader::load(uint32_t *out, uint32_t count) {
    dataset.select({pos}, {count}).read(out, datatype);
    pos += count;
    return count;
}

H5StringReader::H5StringReader(const HighFive::Group &group, std::string path) {
    group.getDataSet(path).read(data);
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
    return UIntWriter(
        std::make_unique<H5UIntWriter>(group, name, chunk_size),
        buffer_size
    );
}
std::unique_ptr<StringWriter> H5WriterBuilder::createStringWriter(std::string name) {
    return std::make_unique<H5StringWriter>(group, name);
}

void H5WriterBuilder::writeVersion(std::string version) {
    group.createAttribute<std::string>("version", HighFive::DataSpace::From(version)).write(version);
}

H5ReaderBuilder::H5ReaderBuilder(std::string file, std::string group, uint32_t buffer_size, uint32_t read_size) :
    group(HighFive::File(
            file,
            HighFive::File::ReadWrite
        ).getGroup(group == "" ? std::string("/") : group)),
    buffer_size(buffer_size), read_size(read_size) {}

UIntReader H5ReaderBuilder::openUIntReader(std::string name) {
    return UIntReader(
        std::make_unique<H5UIntReader>(group, name),
        buffer_size,
        read_size
    );
}

std::unique_ptr<StringReader> H5ReaderBuilder::openStringReader(std::string name) {
    return std::make_unique<H5StringReader>(group, name);
}

std::string H5ReaderBuilder::readVersion() {
    std::string version;
    group.getAttribute("version").read(version);
    return version;
}

}; // end namespace BPCells