#include "binaryfile.h"

namespace BPCells {

ZFileUIntWriter::ZFileUIntWriter(const char* path, uint32_t buffer_size) : buf(buffer_size/sizeof(uint32_t)) {
    // Make sure we get exceptions when things fail
    file.exceptions(std::ofstream::failbit);

    
    // Turn off I/O buffering (Removed because I can't get it to match reasonable performance when I do the buffering manually)
    // file.rdbuf()->pubsetbuf(NULL, 0); 
    
    file.open(path, std::ios_base::binary);
    file.write("UINT32v1", 8);

    write_buffer.data = buf.data();
    write_buffer.size = 0;
}

ZFileUIntWriter::~ZFileUIntWriter() {
    flush();
}

void ZFileUIntWriter::flush() {
    if (write_buffer.data != NULL) {
        uint32_t write_size = data() + capacity() - buf.data();
        for (size_t i = 0; i < write_size; i++) {
            buf[i] = htonl(buf[i]);
        }
        file.write((char *) buf.data(), write_size*sizeof(uint32_t));
    }
    
    write_buffer.data = buf.data();
    write_buffer.size = buf.size();
}
void ZFileUIntWriter::_ensureCapacity(size_t new_capacity) {
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

void ZFileUIntWriter::next() {
    if (data() + capacity() < buf.data() + buf.size()) {
        // We have more buffer left, so don't output yet,
        // just have the client write more data into the buffer
        write_buffer.data += capacity();
        write_buffer.size = buf.data() + buf.size() - write_buffer.data;
    } else {
        flush();
    }
}

ZFileUIntReader::ZFileUIntReader(const char* path, uint32_t buffer_size) : buf(buffer_size/sizeof(uint32_t)) {
    // Turn off I/O buffering (Removed because I can't get it to match reasonable performance when I do the buffering manually)
    // file.rdbuf()->pubsetbuf(NULL, 0); 

    file.open(path, std::ios_base::binary);
    if (!file) {
        throw std::runtime_error(std::string("Error opening file: ") + path);
    }
    
    char header[8];
    file.read(header, 8);

    if (strncmp(header, "UINT32v1", 8) != 0) {
        throw std::invalid_argument(std::string("File header doesn't match magic number (UINT32v1): ") + path);
    }

    // Detect the file size & cache it
    uint32_t cur = file.tellg();
    file.seekg(0, file.end);
    total_size = (file.tellg() / sizeof(uint32_t)) - 2;
    file.seekg(cur);
}

void ZFileUIntReader::_ensureCapacity(size_t new_capacity) {
    if (data() != buf.data()) {
        std::memmove(buf.data(), data(), capacity()*sizeof(uint32_t));
    }
    read_buffer.data = buf.data();
    read_buffer.size += readPos(buf.data() + capacity(), buf.size() - capacity());
    if (capacity() < new_capacity)
        throw std::runtime_error("Not enough remaining data to ensure read capacity");
};

size_t ZFileUIntReader::readPos(uint32_t *out, uint32_t capacity) {
    file.read((char *) out, sizeof(uint32_t)*capacity);
    uint32_t read_size = file.gcount() / sizeof(uint32_t);
    for (uint32_t i = 0; i < read_size; i++) {
        out[i] = ntohl(out[i]);
    }
    return read_size;
}

bool ZFileUIntReader::next() {
    read_buffer.size = readPos(buf.data(), buf.size());
    read_buffer.data = buf.data();
    return read_buffer.size > 0;
}

bool ZFileUIntReader::seek(size_t pos) {
    size_t read_count = data() - capacity() - buf.data();
    size_t buf_pos = (size_t) file.tellg() - read_count;
    if (buf_pos < pos && buf_pos + read_count > pos+128) {
        // Don't seek, just let the next read come from the existing buffer
        read_buffer.size = read_buffer.size - read_count;
        read_buffer.data = buf.data() + read_count;
        return true;
    } else {
        file.seekg(8 + pos * sizeof(uint32_t));
        return next();
    }
}

std::string loadVersionMatrixDir(std::string dir) {
    std::filesystem::path path(dir);

    if (!std::filesystem::exists(path)) {
        throw std::runtime_error(std::string("Missing directory: ") + dir);
    }
    std::vector<std::string> version = readLines(path / "version");
    if (version.size() != 1)
        throw std::runtime_error("Matrix version file does not have exactly one line");
    
    return version[0];
}

UnpackedMatrix openUnpackedMatrixDir(std::string dir, uint32_t buffer_size) {
    std::filesystem::path path(dir);

    if (!std::filesystem::exists(path)) {
        throw std::runtime_error(std::string("Missing directory: ") + dir);
    }

    // Check version info tag
    std::vector<std::string> version = readLines(path / "version");
    if (version.size() != 1 || version[0] != "v1-unpacked")
        throw std::runtime_error("version file has incompatible version (not v1-unpacked)");

    return UnpackedMatrix(
        std::make_unique<ZFileUIntReader>((path / "val").c_str(), buffer_size), 
        std::make_unique<ZFileUIntReader>((path / "row").c_str(), buffer_size), 
        std::make_unique<ZFileUIntReader>((path / "col_ptr").c_str(), buffer_size), 
        std::make_unique<ZFileUIntReader>((path / "row_count").c_str(), buffer_size)
    );
}

UnpackedMatrixWriter createUnpackedMatrixDir(std::string dir, uint32_t buffer_size) {
    std::filesystem::path path(dir);
    if (std::filesystem::exists(path)) {
        throw std::runtime_error(std::string("Path already exists: ") + dir);
    }

    std::filesystem::create_directories(path);

    std::ofstream f((path / "version").c_str());
    f << "v1-unpacked\n";

    return UnpackedMatrixWriter(
        std::make_unique<ZFileUIntWriter>((path / "val").c_str(), buffer_size), 
        std::make_unique<ZFileUIntWriter>((path / "row").c_str(), buffer_size), 
        std::make_unique<ZFileUIntWriter>((path / "col_ptr").c_str(), buffer_size), 
        std::make_unique<ZFileUIntWriter>((path / "row_count").c_str(), buffer_size)
    );
}

PackedMatrix openPackedMatrixDir(std::string dir, uint32_t buffer_size) {
    std::filesystem::path path(dir);

    if (!std::filesystem::exists(path)) {
        throw std::runtime_error(std::string("Missing directory: ") + dir);
    }

    // Check version info tag
    std::vector<std::string> version = readLines(path / "version");
    if (version.size() != 1 || version[0] != "v1-packed")
        throw std::runtime_error("version file has incompatible version (not v1-packed)");

    return PackedMatrix(
        std::make_unique<ZFileUIntReader>((path / "val_data").c_str(), buffer_size), 
        std::make_unique<ZFileUIntReader>((path / "val_idx").c_str(), buffer_size), 
        std::make_unique<ZFileUIntReader>((path / "row_data").c_str(), buffer_size), 
        std::make_unique<ZFileUIntReader>((path / "row_starts").c_str(), buffer_size), 
        std::make_unique<ZFileUIntReader>((path / "row_idx").c_str(), buffer_size), 
        std::make_unique<ZFileUIntReader>((path / "col_ptr").c_str(), buffer_size), 
        std::make_unique<ZFileUIntReader>((path / "row_count").c_str(), buffer_size)
    );
}

PackedMatrixWriter createPackedMatrixDir(std::string dir, uint32_t buffer_size) {
    std::filesystem::path path(dir);
    if (std::filesystem::exists(path)) {
        throw std::runtime_error(std::string("Path already exists: ") + dir);
    }

    std::filesystem::create_directories(path);

    std::ofstream f((path / "version").c_str());
    f << "v1-packed\n";

    return PackedMatrixWriter(
        std::make_unique<ZFileUIntWriter>((path / "val_data").c_str(), buffer_size), 
        std::make_unique<ZFileUIntWriter>((path / "val_idx").c_str(), buffer_size), 
        std::make_unique<ZFileUIntWriter>((path / "row_data").c_str(), buffer_size), 
        std::make_unique<ZFileUIntWriter>((path / "row_starts").c_str(), buffer_size), 
        std::make_unique<ZFileUIntWriter>((path / "row_idx").c_str(), buffer_size), 
        std::make_unique<ZFileUIntWriter>((path / "col_ptr").c_str(), buffer_size), 
        std::make_unique<ZFileUIntWriter>((path / "row_count").c_str(), buffer_size)
    );
}




FileUIntWriter::FileUIntWriter(const char* path) {
    file.exceptions(std::ofstream::failbit);
    file.open(path, std::ios_base::binary);
    file.write("UINT32v1", 8);
}

void FileUIntWriter::write(const uint32_t *buffer, uint32_t count) {    
    size_t i;    
    uint32_t x;
    for (i = 0; i + 128 <= count;) {
        for (size_t j = 0; j < 128; j++) {
            order_buf[j] = htonl(buffer[i++]);
        }
        file.write((char *) &order_buf, sizeof(uint32_t)*128);
    }
    for(; i < count; i++) {
        x = htonl(buffer[i]);
        file.write((char *) &x, sizeof(x));
    }
}

void FileUIntWriter::finalize() {file.flush();}


FileUIntReader::FileUIntReader(const char* path) {
    file.open(path, std::ios_base::binary);
    if (!file) {
        throw std::runtime_error(std::string("Error opening file: ") + path);
    }
    
    char header[8];
    file.read(header, 8);

    if (strncmp(header, "UINT32v1", 8) != 0) {
        throw std::invalid_argument(std::string("File header doesn't match magic number (UINT32v1): ") + path);
    }   
}


uint32_t FileUIntReader::read(uint32_t *buffer, uint32_t count) {
    file.read((char *) buffer, sizeof(uint32_t)*count);
    const uint32_t read_count = file.gcount() / sizeof(uint32_t);

    for (uint32_t i = 0; i < read_count; i++) {
        buffer[i] = ntohl(buffer[i]);
    }
    return read_count;
}

uint32_t FileUIntReader::size() {
    uint32_t cur = file.tellg();
    file.seekg(0, file.end);
    uint32_t size = (file.tellg() / sizeof(uint32_t)) - 2;
    file.seekg(cur);
    return size;
}

void FileUIntReader::seek(const size_t pos) {
    file.seekg(8 + pos * sizeof(uint32_t));
}



FileFragmentsSaver::FileFragmentsSaver(std::string dir) : dir(dir) {
    if (std::filesystem::exists(this->dir)) {
        throw std::runtime_error(std::string("Path already exists: ") + dir);
    }

    std::filesystem::create_directories(this->dir);
    std::filesystem::create_directory(this->dir / "chromosomes");

    std::ofstream f((this->dir / "version").c_str());
    f << "v1\n";
}

UnpackedFrags<FileUIntWriter> FileFragmentsSaver::chrWriterUnpacked(uint32_t chr_id) {
    auto chr_dir = this->dir / "chromosomes" / std::to_string(chr_id);
    std::filesystem::create_directory(chr_dir);
    // T start, end, cell, end_max; 
    return UnpackedFrags<FileUIntWriter> {
        FileUIntWriter((chr_dir / "start").c_str()),
        FileUIntWriter((chr_dir / "end").c_str()),
        FileUIntWriter((chr_dir / "cell").c_str()),
        FileUIntWriter((chr_dir / "end_max").c_str())
    };
}

PackedFrags<FileUIntWriter> FileFragmentsSaver::chrWriterPacked(uint32_t chr_id) {
    auto chr_dir = this->dir / "chromosomes" / std::to_string(chr_id);
    std::filesystem::create_directory(chr_dir);
    //    T start_data, start_idx, start_starts,
    //         end_data, end_idx, end_max,
    //         cell_data, cell_idx, count;
    return PackedFrags<FileUIntWriter> {
        FileUIntWriter((chr_dir / "start_data").c_str()),
        FileUIntWriter((chr_dir / "start_idx").c_str()),
        FileUIntWriter((chr_dir / "start_starts").c_str()),
        FileUIntWriter((chr_dir / "end_data").c_str()),
        FileUIntWriter((chr_dir / "end_idx").c_str()),
        FileUIntWriter((chr_dir / "end_max").c_str()),
        FileUIntWriter((chr_dir / "cell_data").c_str()),
        FileUIntWriter((chr_dir / "cell_idx").c_str()),
        FileUIntWriter((chr_dir / "count").c_str())
    };
}

void FileFragmentsSaver::writeCellNames(std::vector<std::string> cell_names) {
    std::ofstream f((this->dir / "cell_names.txt").c_str());
    for (auto n : cell_names) {
        if (n.find('\n') != std::string::npos) {
            throw std::runtime_error("Cell names must not contain newline characters");
        }
        f << n << '\n';
    }
}
void FileFragmentsSaver::writeChrNames(std::vector<std::string> chr_names) {
    std::ofstream f((this->dir / "chr_names.txt").c_str());
    for (auto n : chr_names) {
        if (n.find('\n') != std::string::npos) {
            throw std::runtime_error("Chromosme names must not contain newline characters");
        }
        f << n << '\n';
    }
}

FileFragmentsLoader::FileFragmentsLoader(std::string dir, bool is_packed) : dir(dir) {
    // Check that all the necessary files exist
    if (!std::filesystem::exists(dir)) {
        throw std::runtime_error(std::string("Missing directory: ") + dir);
    }

    std::vector<std::string> version = readLines(this->dir / "version");
    if (version.size() != 1 || version[0] != "v1")
        throw std::runtime_error("version file has incompatible version (not v1)");
    

    if (!std::filesystem::exists(this->dir / "cell_names.txt")) {
        throw std::runtime_error(std::string("Missing file: ") + (this->dir / "cell_names.txt").string());
    }

    // Reading the chr_names so we know how many to expect
    uint32_t chr_count = readChrNames().size();
    for (uint32_t i = 0; i < chr_count; i++) {
        std::filesystem::path d(this->dir / "chromosomes" / std::to_string(i));
        if (!std::filesystem::exists(d)) {
            throw std::runtime_error(std::string("Missing directory: ") + d.string());
        }
        
        if (is_packed) {
            const char *expected_files[] = {"start_data", "start_idx", "start_starts", "end_data", "end_idx", "end_max", "cell_data", "cell_idx", "count"};
            for (auto f : expected_files) {
                if (!std::filesystem::exists(d / f)) {
                    throw std::runtime_error(std::string("Missing file: ") + (d / f).c_str());
                }
            }
        } else {
            const char *expected_files[] = {"start", "end", "cell", "end_max"};
            for (auto f : expected_files) {
                if (!std::filesystem::exists(d / f)) {
                    throw std::runtime_error(std::string("Missing file: ") + (d / f).c_str());
                }
            }
        }
    }
}

UnpackedFrags<FileUIntReader> FileFragmentsLoader::chrReaderUnpacked(uint32_t chr_id) {
    std::filesystem::path d(this->dir / "chromosomes" / std::to_string(chr_id));
    return UnpackedFrags<FileUIntReader> {
        FileUIntReader((d / "start").c_str()),
        FileUIntReader((d / "end").c_str()),
        FileUIntReader((d / "cell").c_str()),
        FileUIntReader((d / "end_max").c_str())
    };
}

PackedFrags<FileUIntReader> FileFragmentsLoader::chrReaderPacked(uint32_t chr_id) {
    std::filesystem::path d(this->dir / "chromosomes" / std::to_string(chr_id));
    return PackedFrags<FileUIntReader> {
        FileUIntReader((d / "start_data").c_str()),
        FileUIntReader((d / "start_idx").c_str()),
        FileUIntReader((d / "start_starts").c_str()),
        FileUIntReader((d / "end_data").c_str()),
        FileUIntReader((d / "end_idx").c_str()),
        FileUIntReader((d / "end_max").c_str()),
        FileUIntReader((d / "cell_data").c_str()),
        FileUIntReader((d / "cell_idx").c_str()),
        FileUIntReader((d / "count").c_str())
    };
}

std::vector<std::string> readLines(std::filesystem::path path) {
    std::ifstream in;
    std::string line;
    std::vector<std::string> ret;

    in.open(path.c_str());
    if (!in) 
        throw std::runtime_error(std::string("Could not open file: ") + path.c_str());
    
    while (std::getline(in, line)) {
        ret.push_back(line);
    }
    return ret;
}


} // end namespace BPCells