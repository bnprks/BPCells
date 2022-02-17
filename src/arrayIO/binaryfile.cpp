#include "binaryfile.h"

namespace BPCells {

// This should equal UINT32v1 when printed on a little endian system
const static uint32_t FileHeader[2] = {0x544e4955, 0x31763233};

FileUIntWriter::FileUIntWriter(const char* path) {
    // Make sure we get exceptions when things fail
    file.exceptions(std::ofstream::failbit | std::ofstream::badbit);
    
    // Turn off I/O buffering (Removed because I can't get it to match reasonable performance when I do the buffering manually)
    // file.rdbuf()->pubsetbuf(NULL, 0); 
    
    file.open(path, std::ios_base::binary);
    file.write((char*) FileHeader, 8);
}

uint32_t FileUIntWriter::write(uint32_t *in, uint32_t count) {
    file.write((char *) in, count*sizeof(uint32_t));
    return count;
}

FileUIntReader::FileUIntReader(const char* path) {
    // Turn off I/O buffering (Removed because I can't get it to match reasonable performance when I do the buffering manually)
    // file.rdbuf()->pubsetbuf(NULL, 0); 

    file.open(path, std::ios_base::binary);
    if (!file) {
        throw std::runtime_error(std::string("Error opening file: ") + path);
    }
    
    uint32_t header[2];
    file.read((char *) header, 8);
    if (header[0] == FileHeader[0] && header[1] == FileHeader[1]) {
        byte_swap = false;
    } else if (__builtin_bswap32(header[0]) == FileHeader[0] &&
               __builtin_bswap32(header[1]) == FileHeader[1] ) {
        byte_swap = true;
    } else {
        throw std::invalid_argument(std::string("File header doesn't match magic number (UINT32v1 or byteswapped TNIU1v23): ") + path);
    }

    // Detect the file size & cache it
    uint32_t cur = file.tellg();
    file.seekg(0, file.end);
    total_size = (file.tellg() / sizeof(uint32_t)) - 2;
    file.seekg(cur);
}

uint32_t FileUIntReader::size() const {return total_size;}

void FileUIntReader::seek(uint32_t pos) {
    file.seekg(8 + pos * sizeof(uint32_t));
}

uint32_t FileUIntReader::load(uint32_t *out, uint32_t count) {
    file.read((char *) out, sizeof(uint32_t)*count);
    uint32_t read_count = file.gcount() / sizeof(uint32_t);
    if (byte_swap) {
        for (uint32_t i = 0; i < read_count; i++) {
            out[i] = __builtin_bswap32(out[i]);
        }
    }
    return read_count;
}

FileStringReader::FileStringReader(std::filesystem::path path) : data(readLines(path)) {}
const char* FileStringReader::get(uint32_t idx) const {
    if (idx < data.size()) return data[idx].c_str();
    return NULL;
}

uint32_t FileStringReader::size() const {return data.size();}

FileStringWriter::FileStringWriter(std::filesystem::path path) : path(path) {}
void FileStringWriter::write(const StringReader &reader) {
    std::ofstream f(path.c_str());
    uint32_t i = 0;
    while (true) {
        const char* s = reader.get(i);
        if (s == NULL) break;
        while(*s != '\0') {
            f.put(*s);
            s++;
        }
        f.put('\n');
        i += 1;
    }
}

FileWriterBuilder::FileWriterBuilder(std::string _dir, uint32_t buffer_size) :
    dir(_dir), buffer_size(buffer_size) {

    if (std::filesystem::exists(dir)) {
        throw std::runtime_error(std::string("Path already exists: ") + _dir);
    }

    std::filesystem::create_directories(dir);
}

UIntWriter FileWriterBuilder::createUIntWriter(std::string name) {
    return UIntWriter(
        std::make_unique<FileUIntWriter>((dir / name).c_str()), 
        buffer_size
    );
}

std::unique_ptr<StringWriter> FileWriterBuilder::createStringWriter(std::string name) {
    return std::make_unique<FileStringWriter>(dir / name);
}

void FileWriterBuilder::writeVersion(std::string version) {
    std::ofstream f((dir / "version").c_str());
    f << version << std::endl;
}

FileReaderBuilder::FileReaderBuilder(std::string _dir, uint32_t buffer_size, uint32_t read_size) :
    dir(_dir), buffer_size(buffer_size), read_size(read_size) {
    
    if (!std::filesystem::exists(dir)) {
        throw std::invalid_argument(std::string("Missing directory: ") + _dir);
    }
}

UIntReader FileReaderBuilder::openUIntReader(std::string name) {
    return UIntReader(
        std::make_unique<FileUIntReader>((dir / name).c_str()),
        buffer_size, read_size
    );
}

std::unique_ptr<StringReader> FileReaderBuilder::openStringReader(std::string name) {
    return std::make_unique<FileStringReader>(dir / name);
}

std::string FileReaderBuilder::readVersion() {
    std::vector<std::string> version = readLines(dir / "version");
    if (version.size() != 1)
        throw std::runtime_error("Version file does not have exactly one line");
    
    return version[0];
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