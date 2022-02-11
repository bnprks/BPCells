#include "hdf5.h"
#include <hdf5.h>


namespace BPCells {



MovableDataType::MovableDataType(const HighFive::DataType& other) : HighFive::DataType(other) {}
MovableDataType::MovableDataType(HighFive::DataType&& other) : HighFive::DataType(std::move(other)) {}

ZH5UIntWriter::ZH5UIntWriter(std::string file_path, std::string group_path, uint32_t buffer_size, uint32_t chunk_size) : 
    dataset(constructH5DataSet(file_path, group_path, chunk_size)), buf(buffer_size/sizeof(uint32_t)) {

    write_buffer.data = buf.data();
    write_buffer.size = 0;
}
ZH5UIntWriter::ZH5UIntWriter(const HighFive::Group &group, std::string path, uint32_t buffer_size, uint32_t chunk_size) : 
    dataset(constructH5DataSet(group, path, chunk_size)), buf(buffer_size/sizeof(uint32_t)) {

    write_buffer.data = buf.data();
    write_buffer.size = 0;
}

ZH5UIntWriter::~ZH5UIntWriter() {
    flush(); 
}

void ZH5UIntWriter::flush() {
    if (write_buffer.data != NULL) {
        uint32_t write_size = data() + capacity() - buf.data();

        uint32_t cur_size = dataset.getDimensions()[0];
        dataset.resize({cur_size + write_size});
        dataset.select({cur_size}, {write_size}).write_raw(buf.data(), datatype);
    }
    write_buffer.data = buf.data();
    write_buffer.size = buf.size();
}

void ZH5UIntWriter::_ensureCapacity(size_t new_capacity) {
    if (buf.data() + buf.size() - data() > new_capacity) {
        //We have enough space before the end of the buffer, so just
        //expand our advertized capacity
        write_buffer.size = buf.data() + buf.size() - data();
    } else if (new_capacity > buf.size()) {
        throw std::runtime_error("Requested write capacity larger than buffer size");
    } else {
        write_buffer.size = 0;
        flush();
    }
}

// This is basically identical to the FileUIntWriter implementation
void ZH5UIntWriter::next() {
    if (data() + capacity() < buf.data() + buf.size()) {
        // We have more buffer left, so don't output yet,
        // just have the client write more data into the buffer
        write_buffer.data += capacity();
        write_buffer.size = buf.data() + buf.size() - write_buffer.data;
    } else {
        flush();
    }
}

ZH5UIntReader::ZH5UIntReader(std::string file_path, std::string group_path, uint32_t buffer_size) :
    dataset(openH5DataSet(file_path, group_path)), buf(buffer_size/sizeof(uint32_t)) {

    total_size = dataset.getDimensions()[0];
}
ZH5UIntReader::ZH5UIntReader(const HighFive::Group &group, std::string path, uint32_t buffer_size) :
    dataset(group.getDataSet(path)), buf(buffer_size/sizeof(uint32_t)) {

    total_size = dataset.getDimensions()[0];
}

void ZH5UIntReader::_ensureCapacity(size_t new_capacity) {
    if (data() != buf.data()) {
        std::memmove(buf.data(), data(), capacity() * sizeof(uint32_t));
    }
    read_buffer.data = buf.data();
    read_buffer.size += readPos(buf.data() + capacity(), buf.size() - capacity());
    if (capacity() < new_capacity)
        throw std::runtime_error("Not enough remaining data to ensure read capacity");
};

size_t ZH5UIntReader::readPos(uint32_t *out, uint32_t capacity) {
    uint32_t read_count = std::min(capacity, (uint32_t) (total_size - pos));
    dataset.select({pos}, {read_count}).read(out, datatype);
    pos += read_count;
    return read_count;
}

bool ZH5UIntReader::next() {
    read_buffer.size = readPos(buf.data(), buf.size());
    read_buffer.data = buf.data();
    return read_buffer.size > 0;
}
bool ZH5UIntReader::seek(size_t new_pos) {
    // Number of successfully read items in the buffer
    size_t read_count = data() - capacity() - buf.data();
    // File position corresponding to the start of the buffer
    auto buf_pos = pos - read_count;
    if (buf_pos <= new_pos && buf_pos + read_count > new_pos) {
        // Don't seek, just let the next read come from the existing buffer
        read_buffer.data = buf.data() + pos - buf_pos;
        read_buffer.size = read_count - (pos - buf_pos);
        return true;
    } else {
        pos = new_pos;
        return next();
    }
}

H5StringReader::H5StringReader(const MovableDataSet &dataset) {
    dataset.read(data);
}
const char* H5StringReader::get(uint32_t idx) const {
    if (idx < data.size()) return data[idx].c_str();
    return NULL;
}
uint32_t H5StringReader::size() const {return data.size();}

H5StringWriter::H5StringWriter(std::string file, std::string group) : file(file), group(group) {}
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
    HighFive::File f(file, HighFive::File::ReadWrite);
    HighFive::DataSet ds = f.createDataSet<std::string>(group, HighFive::DataSpace::From(data));
    ds.write(data);
}

std::string loadVersionH5(std::string file_path, std::string group_path) {
    HighFive::SilenceHDF5 s;
    if (group_path == "") group_path = "/";
    
    // Load with ReadWrite permissions so that if we open an HDF5 then later want to write
    // back to it, we aren't locked out from having originally opened the file ReadOnly
    HighFive::Group group(HighFive::File(file_path, HighFive::File::ReadWrite).getGroup(group_path));

    std::string version;
    group.getAttribute("version").read(version);
    return version;
}

UnpackedMatrix open10xFeatureMatrix(std::string file_path, std::string group_path, uint32_t buffer_size) {
    HighFive::SilenceHDF5 s;
    if (group_path == "") group_path = "/";
    
    // Load with ReadWrite permissions so that if we open an HDF5 then later want to write
    // back to it, we aren't locked out from having originally opened the file ReadOnly
    HighFive::Group group(HighFive::File(file_path, HighFive::File::ReadWrite).getGroup(group_path));

    return UnpackedMatrix(
        std::make_unique<ZH5UIntReader>(group, "indices", buffer_size), 
        std::make_unique<ZH5UIntReader>(group, "data", buffer_size), 
        std::make_unique<ZH5UIntReader>(group, "indptr", buffer_size), 
        std::make_unique<ZH5UIntReader>(group, "shape", buffer_size)
    );
}

UnpackedMatrix openUnpackedMatrixH5(std::string file_path, std::string group_path, uint32_t buffer_size) {
    HighFive::SilenceHDF5 s;
    if (group_path == "") group_path = "/";
    
    // Load with ReadWrite permissions so that if we open an HDF5 then later want to write
    // back to it, we aren't locked out from having originally opened the file ReadOnly
    HighFive::Group group(HighFive::File(file_path, HighFive::File::ReadWrite).getGroup(group_path));

    std::string version;
    group.getAttribute("version").read(version);
    if (version != "v1-unpacked-matrix")
        throw std::runtime_error("HDF5 group does not have correct version attribute (v1-unpacked-matrix)");

    return UnpackedMatrix(
        std::make_unique<ZH5UIntReader>(group, "val", buffer_size), 
        std::make_unique<ZH5UIntReader>(group, "row", buffer_size), 
        std::make_unique<ZH5UIntReader>(group, "col_ptr", buffer_size), 
        std::make_unique<ZH5UIntReader>(group, "row_count", buffer_size)
    );
}

UnpackedMatrixWriter createUnpackedMatrixH5(std::string file_path, std::string group_path, uint32_t buffer_size, uint32_t chunk_size) {
    HighFive::SilenceHDF5 s;
    if (group_path == "") group_path = "/";
    
    std::filesystem::path path(file_path);
    if (path.has_parent_path() && !std::filesystem::exists(path.parent_path())) {
        std::filesystem::create_directories(path.parent_path());
    }
    
    MovableGroup group = constructH5Group(file_path, group_path);

    std::string version("v1-unpacked-matrix");
    group.createAttribute<std::string>("version", HighFive::DataSpace::From(version)).write(version);

    return UnpackedMatrixWriter(
        std::make_unique<ZH5UIntWriter>(group, "val", buffer_size), 
        std::make_unique<ZH5UIntWriter>(group, "row", buffer_size), 
        std::make_unique<ZH5UIntWriter>(group, "col_ptr", buffer_size), 
        std::make_unique<ZH5UIntWriter>(group, "row_count", buffer_size)
    );
}


PackedMatrix openPackedMatrixH5(std::string file_path, std::string group_path, uint32_t buffer_size) {
    HighFive::SilenceHDF5 s;
    if (group_path == "") group_path = "/";
    
    // Load with ReadWrite permissions so that if we open an HDF5 then later want to write
    // back to it, we aren't locked out from having originally opened the file ReadOnly
    HighFive::Group group(HighFive::File(file_path, HighFive::File::ReadWrite).getGroup(group_path));

    std::string version;
    group.getAttribute("version").read(version);
    if (version != "v1-packed-matrix")
        throw std::runtime_error("HDF5 group does not have correct version attribute (v1-packed-matrix)");

    return PackedMatrix(
        std::make_unique<ZH5UIntReader>(group, "val_data", buffer_size), 
        std::make_unique<ZH5UIntReader>(group, "val_idx", buffer_size), 
        std::make_unique<ZH5UIntReader>(group, "row_data", buffer_size), 
        std::make_unique<ZH5UIntReader>(group, "row_starts", buffer_size), 
        std::make_unique<ZH5UIntReader>(group, "row_idx", buffer_size), 
        std::make_unique<ZH5UIntReader>(group, "col_ptr", buffer_size), 
        std::make_unique<ZH5UIntReader>(group, "row_count", buffer_size)
    );
}

PackedMatrixWriter createPackedMatrixH5(std::string file_path, std::string group_path, uint32_t buffer_size, uint32_t chunk_size) {
    HighFive::SilenceHDF5 s;
    if (group_path == "") group_path = "/";
    
    std::filesystem::path path(file_path);
    if (path.has_parent_path() && !std::filesystem::exists(path.parent_path())) {
        std::filesystem::create_directories(path.parent_path());
    }
    
    MovableGroup group = constructH5Group(file_path, group_path);
    
    std::string version("v1-packed-matrix");
    group.createAttribute<std::string>("version", HighFive::DataSpace::From(version)).write(version);

    return PackedMatrixWriter(
        std::make_unique<ZH5UIntWriter>(group, "val_data", buffer_size), 
        std::make_unique<ZH5UIntWriter>(group, "val_idx", buffer_size), 
        std::make_unique<ZH5UIntWriter>(group, "row_data", buffer_size), 
        std::make_unique<ZH5UIntWriter>(group, "row_starts", buffer_size), 
        std::make_unique<ZH5UIntWriter>(group, "row_idx", buffer_size), 
        std::make_unique<ZH5UIntWriter>(group, "col_ptr", buffer_size), 
        std::make_unique<ZH5UIntWriter>(group, "row_count", buffer_size)
    );
}


UnpackedFragments3 openUnpackedFragmentsH5(std::string file_path, std::string group_path, uint32_t buffer_size) {
    HighFive::SilenceHDF5 s;
    if (group_path == "") group_path = "/";
    
    // Load with ReadWrite permissions so that if we open an HDF5 then later want to write
    // back to it, we aren't locked out from having originally opened the file ReadOnly
    HighFive::Group group(HighFive::File(file_path, HighFive::File::ReadWrite).getGroup(group_path));

    std::string version;
    group.getAttribute("version").read(version);
    if (version != "v1-unpacked-fragments")
        throw std::runtime_error("HDF5 group does not have correct version attribute (v1-unpacked-fragments)");

    return UnpackedFragments3(
        std::make_unique<ZH5UIntReader>(group, "cell", buffer_size), 
        std::make_unique<ZH5UIntReader>(group, "start", buffer_size), 
        std::make_unique<ZH5UIntReader>(group, "end", buffer_size), 
        std::make_unique<ZH5UIntReader>(group, "end_max", buffer_size), 
        std::make_unique<ZH5UIntReader>(group, "chr_ptr", buffer_size), 
        std::make_unique<H5StringReader>(openH5DataSet(group, "chr_names")), 
        std::make_unique<H5StringReader>(openH5DataSet(group, "cell_names"))
    );
}

UnpackedFragmentsWriter3 createUnpackedFragmentsH5(std::string file_path, std::string group_path, uint32_t buffer_size, uint32_t chunk_size) {
    HighFive::SilenceHDF5 s;
    if (group_path == "") group_path = "/";
    
    std::filesystem::path path(file_path);
    if (path.has_parent_path() && !std::filesystem::exists(path.parent_path())) {
        std::filesystem::create_directories(path.parent_path());
    }
    
    MovableGroup group = constructH5Group(file_path, group_path);

    std::string version("v1-unpacked-fragments");
    group.createAttribute<std::string>("version", HighFive::DataSpace::From(version)).write(version);

    return UnpackedFragmentsWriter3(
        std::make_unique<ZH5UIntWriter>(group, "cell", buffer_size), 
        std::make_unique<ZH5UIntWriter>(group, "start", buffer_size), 
        std::make_unique<ZH5UIntWriter>(group, "end", buffer_size), 
        std::make_unique<ZH5UIntWriter>(group, "end_max", buffer_size), 
        std::make_unique<ZH5UIntWriter>(group, "chr_ptr", buffer_size), 
        std::make_unique<H5StringWriter>(file_path, group.getPath() + "/chr_names"), 
        std::make_unique<H5StringWriter>(file_path, group.getPath() + "/cell_names")
    );
}

PackedFragments3 openPackedFragmentsH5(std::string file_path, std::string group_path, uint32_t buffer_size) {
    HighFive::SilenceHDF5 s;
    if (group_path == "") group_path = "/";
    
    // Load with ReadWrite permissions so that if we open an HDF5 then later want to write
    // back to it, we aren't locked out from having originally opened the file ReadOnly
    HighFive::Group group(HighFive::File(file_path, HighFive::File::ReadWrite).getGroup(group_path));

    std::string version;
    group.getAttribute("version").read(version);
    if (version != "v1-packed-fragments")
        throw std::runtime_error("HDF5 group does not have correct version attribute (v1-packed-fragments)");

    return PackedFragments3(
        std::make_unique<ZH5UIntReader>(group, "cell_data", buffer_size),
        std::make_unique<ZH5UIntReader>(group, "cell_idx", buffer_size),
        std::make_unique<ZH5UIntReader>(group, "start_data", buffer_size),
        std::make_unique<ZH5UIntReader>(group, "start_idx", buffer_size),
        std::make_unique<ZH5UIntReader>(group, "start_starts", buffer_size),
        std::make_unique<ZH5UIntReader>(group, "end_data", buffer_size),
        std::make_unique<ZH5UIntReader>(group, "end_idx", buffer_size),
        std::make_unique<ZH5UIntReader>(group, "end_max", buffer_size),
        std::make_unique<ZH5UIntReader>(group, "chr_ptr", buffer_size),

        std::make_unique<H5StringReader>(openH5DataSet(group, "chr_names")), 
        std::make_unique<H5StringReader>(openH5DataSet(group, "cell_names"))
    );
}

PackedFragmentsWriter3 createPackedFragmentsH5(std::string file_path, std::string group_path, uint32_t buffer_size, uint32_t chunk_size) {
    HighFive::SilenceHDF5 s;
    if (group_path == "") group_path = "/";
    
    std::filesystem::path path(file_path);
    if (path.has_parent_path() && !std::filesystem::exists(path.parent_path())) {
        std::filesystem::create_directories(path.parent_path());
    }
    
    MovableGroup group = constructH5Group(file_path, group_path);

    std::string version("v1-packed-fragments");
    group.createAttribute<std::string>("version", HighFive::DataSpace::From(version)).write(version);

    return PackedFragmentsWriter3(
        std::make_unique<ZH5UIntWriter>(group, "cell_data", buffer_size),
        std::make_unique<ZH5UIntWriter>(group, "cell_idx", buffer_size),
        std::make_unique<ZH5UIntWriter>(group, "start_data", buffer_size),
        std::make_unique<ZH5UIntWriter>(group, "start_idx", buffer_size),
        std::make_unique<ZH5UIntWriter>(group, "start_starts", buffer_size),
        std::make_unique<ZH5UIntWriter>(group, "end_data", buffer_size),
        std::make_unique<ZH5UIntWriter>(group, "end_idx", buffer_size),
        std::make_unique<ZH5UIntWriter>(group, "end_max", buffer_size),
        std::make_unique<ZH5UIntWriter>(group, "chr_ptr", buffer_size),

        std::make_unique<H5StringWriter>(file_path, group.getPath() + "/chr_names"), 
        std::make_unique<H5StringWriter>(file_path, group.getPath() + "/cell_names")
    );
}



H5UIntWriter::H5UIntWriter(std::string file_path, std::string group_path, uint32_t chunk_size) :
        dataset(constructH5DataSet(file_path, group_path, chunk_size)) {}

H5UIntWriter::H5UIntWriter(HighFive::Group group, std::string path, uint32_t chunk_size) :
        dataset(constructH5DataSet(group, path, chunk_size)) {}

// Append integers to output stream; Throws exception on failure
void H5UIntWriter::write(const uint32_t *buffer, uint32_t count) {
    uint32_t cur_size = dataset.getDimensions()[0];
    dataset.resize({cur_size + count});
    dataset.select({cur_size}, {count}).write_raw(buffer, datatype);
};





MovableDataSet constructH5DataSet(std::string file_path, std::string group_path, uint32_t chunk_size) {
    HighFive::SilenceHDF5 s;
    HighFive::File file(file_path, HighFive::File::ReadWrite | HighFive::File::OpenOrCreate);

    // Create a dataspace with initial shape and max shape
    HighFive::DataSpace dataspace({0}, {HighFive::DataSpace::UNLIMITED});

    // Use chunking
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(std::vector<hsize_t>{chunk_size}));
    HighFive::DataSetAccessProps a_props;
    a_props.add(HighFive::Caching(521, 50<<20));// 50MB cache for overkill
    
    // Create the dataset
    return MovableDataSet(file.createDataSet<uint32_t>(group_path, dataspace, props, a_props));
}


MovableDataSet constructH5DataSet(HighFive::Group group, std::string group_path, uint32_t chunk_size) {
    HighFive::SilenceHDF5 s;
    // Create a dataspace with initial shape and max shape
    HighFive::DataSpace dataspace({0}, {HighFive::DataSpace::UNLIMITED});

    // Use chunking
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(std::vector<hsize_t>{chunk_size}));
    HighFive::DataSetAccessProps a_props;
    a_props.add(HighFive::Caching(521, 50<<20));// 50MB cache for overkill
    // Create the dataset
    return MovableDataSet(group.createDataSet<uint32_t>(group_path, dataspace, props, a_props));
}

void H5UIntWriter::finalize() {dataset.getFile().flush();} 


H5UIntReader::H5UIntReader(std::string file_path, std::string group_path) :
    dataset(openH5DataSet(file_path, group_path)) {}
H5UIntReader::H5UIntReader(HighFive::Group group, std::string path) :
    dataset(group.getDataSet(path)) {}

MovableDataSet openH5DataSet(std::string file_path, std::string group_path) {
    HighFive::SilenceHDF5 s;
    HighFive::DataSetAccessProps a_props;
    a_props.add(HighFive::Caching(521, 50<<20));// 50MB cache for overkill

    // Load with ReadWrite permissions so that if we open an HDF5 then later want to write
    // back to it, we aren't locked out from having originally opened the file ReadOnly
    return MovableDataSet(HighFive::File(file_path, HighFive::File::ReadWrite).getDataSet(group_path, a_props));
}

MovableDataSet openH5DataSet(HighFive::Group group, std::string path) {
    HighFive::SilenceHDF5 s;
    HighFive::DataSetAccessProps a_props;
    a_props.add(HighFive::Caching(521, 50<<20));// 50MB cache for overkill
    return MovableDataSet(group.getDataSet(path, a_props));
}
// Read integers up to count into buffer and return the number actually read. Throw an
// exception on failures. Must repeatedly return 0 after finishing the stream
// (and not return 0 before then)
uint32_t H5UIntReader::read(uint32_t *buffer, uint32_t count) {
    if (pos + count > dataset.getDimensions()[0]) {
        count = dataset.getDimensions()[0] - pos;
    }
    dataset.select({pos}, {count}).read(buffer, datatype);
    pos += count;
    return count;
}

// If this reader is seekable, return the total number of integers stored
uint32_t H5UIntReader::size() {return dataset.getDimensions()[0];};
// Seek to the given index in the array
void H5UIntReader::seek(const size_t pos) {this->pos = pos;};

H5FragmentsSaver::H5FragmentsSaver(std::string file_path, std::string group_path, uint32_t chunk_size) :
    group(constructH5Group(file_path, group_path)),
    chunk_size(chunk_size) {
    
    std::string version("v1");
    group.createAttribute<std::string>("version", HighFive::DataSpace::From(version)).write(version);
}

MovableGroup constructH5Group(std::string file_path, std::string group_path) {
    HighFive::SilenceHDF5 s;
    if (group_path == "") group_path = "/";
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

UnpackedFrags<H5UIntWriter> H5FragmentsSaver::chrWriterUnpacked(uint32_t chr_id) {
    HighFive::Group chr = group.createGroup(std::string("chromosomes/") + std::to_string(chr_id));
    return UnpackedFrags<H5UIntWriter> {
        H5UIntWriter(chr, "start", chunk_size),
        H5UIntWriter(chr, "end", chunk_size),
        H5UIntWriter(chr, "cell", chunk_size),
        H5UIntWriter(chr, "end_max", chunk_size)
    };
}

PackedFrags<H5UIntWriter> H5FragmentsSaver::chrWriterPacked(uint32_t chr_id) {
    HighFive::Group chr = group.createGroup(std::string("chromosomes/") + std::to_string(chr_id));
    return PackedFrags<H5UIntWriter> {
        H5UIntWriter(chr, "start_data", chunk_size),
        H5UIntWriter(chr, "start_idx", chunk_size),
        H5UIntWriter(chr, "start_starts", chunk_size),
        H5UIntWriter(chr, "end_data", chunk_size),
        H5UIntWriter(chr, "end_idx", chunk_size),
        H5UIntWriter(chr, "end_max", chunk_size),
        H5UIntWriter(chr, "cell_data", chunk_size),
        H5UIntWriter(chr, "cell_idx", chunk_size),
        H5UIntWriter(chr, "count", chunk_size)
    };
}

void H5FragmentsSaver::writeCellNames(std::vector<std::string> cell_names) {
    HighFive::SilenceHDF5 s;
    group.createDataSet<std::string>(
        "cell_names", HighFive::DataSpace::From(cell_names)
    ).write(cell_names);
}
void H5FragmentsSaver::writeChrNames(std::vector<std::string> chr_names) {
    HighFive::SilenceHDF5 s;
    group.createDataSet<std::string>(
        "chr_names", HighFive::DataSpace::From(chr_names)
    ).write(chr_names);
}

H5FragmentsLoader::H5FragmentsLoader(std::string file_path, std::string group_path, bool is_packed) :
    group(openH5Group(file_path, group_path)) {
    HighFive::SilenceHDF5 s;
    
    std::string version;
    group.getAttribute("version").read(version);
    if (version != "v1")
        throw std::runtime_error("HDF5 group does not have correct version attribute");

    // Check that the whole group structure exists
    if (!group.exist("cell_names"))
        throw std::runtime_error("Missing group: cell_names");
    if (!group.exist("chr_names"))
        throw std::runtime_error("Missing group: chr_names");
    
    size_t chr_count = group.getDataSet("chr_names").getSpace().getElementCount();
    
    for (uint32_t i = 0; i < chr_count; i++) {
        HighFive::Group chr = group.getGroup(std::string("chromosomes/") + std::to_string(i));
        if (is_packed) {
            const char *expected_groups[] = {"start_data", "start_idx", "start_starts", "end_data", "end_idx", "end_max", "cell_data", "cell_idx", "count"};
            for (auto g : expected_groups) {
                if (!chr.exist(g))
                    throw std::runtime_error(std::string("Missing group: ") + chr.getPath() + "/" + g);
            }
        } else {
            const char *expected_groups[] = {"start", "end", "cell", "end_max"};
            for (auto g : expected_groups) {
                if (!chr.exist(g))
                    throw std::runtime_error(std::string("Missing group: ") + chr.getPath() + "/" + g);
            }
        }
    }
}

MovableGroup openH5Group(std::string file_path, std::string group_path) {
    HighFive::SilenceHDF5 s;
    if (group_path == "") group_path = "/";
    return HighFive::File(file_path, HighFive::File::ReadWrite).getGroup(group_path);
}

UnpackedFrags<H5UIntReader> H5FragmentsLoader::chrReaderUnpacked(uint32_t chr_id) {
    HighFive::Group chr = group.getGroup(std::string("chromosomes/") + std::to_string(chr_id));
    return UnpackedFrags<H5UIntReader> {
        H5UIntReader(chr, "start"),
        H5UIntReader(chr, "end"),
        H5UIntReader(chr, "cell"),
        H5UIntReader(chr, "end_max")
    };
}

PackedFrags<H5UIntReader> H5FragmentsLoader::chrReaderPacked(uint32_t chr_id) {
    HighFive::Group chr = group.getGroup(std::string("chromosomes/") + std::to_string(chr_id));
    return PackedFrags<H5UIntReader> {
        H5UIntReader(chr, "start_data"),
        H5UIntReader(chr, "start_idx"),
        H5UIntReader(chr, "start_starts"),
        H5UIntReader(chr, "end_data"),
        H5UIntReader(chr, "end_idx"),
        H5UIntReader(chr, "end_max"),
        H5UIntReader(chr, "cell_data"),
        H5UIntReader(chr, "cell_idx"),
        H5UIntReader(chr, "count")
    };
}

std::vector<std::string> H5FragmentsLoader::readCellNames() {
    HighFive::SilenceHDF5 s;
    std::vector<std::string> res;
    group.getDataSet("cell_names").read(res);
    return res;
}
std::vector<std::string> H5FragmentsLoader::readChrNames() {
    HighFive::SilenceHDF5 s;
    std::vector<std::string> res;
    HighFive::DataSet d = group.getDataSet("chr_names");
    d.read(res);
    return res;
}

}; // end namespace BPCells