#include "binaryfile.h"

namespace BPCells {

FileStringReader::FileStringReader(std::filesystem::path path) : data(readLines(path)) {}
const char *FileStringReader::get(uint64_t idx) const {
    if (idx < data.size()) return data[idx].c_str();
    return NULL;
}

uint64_t FileStringReader::size() const { return data.size(); }

FileStringWriter::FileStringWriter(std::filesystem::path path) : path(path) {}
void FileStringWriter::write(const StringReader &reader) {
    std::ofstream f(path.c_str());
    uint64_t i = 0;
    while (true) {
        const char *s = reader.get(i);
        if (s == NULL) break;
        while (*s != '\0') {
            f.put(*s);
            s++;
        }
        f.put('\n');
        i += 1;
    }
}

FileWriterBuilder::FileWriterBuilder(std::string _dir, uint64_t buffer_size, bool allow_exists)
    : dir(_dir)
    , buffer_size(buffer_size) {

    if (!allow_exists && std::filesystem::exists(dir)) {
        throw std::runtime_error(std::string("Path already exists: ") + _dir);
    }

    std::filesystem::create_directories(dir);
}

UIntWriter FileWriterBuilder::createUIntWriter(std::string name) {
    return UIntWriter(
        std::make_unique<FileNumWriter<uint32_t>>((dir / name).string().c_str()), buffer_size
    );
}

FloatWriter FileWriterBuilder::createFloatWriter(std::string name) {
    return FloatWriter(
        std::make_unique<FileNumWriter<float>>((dir / name).string().c_str()), buffer_size
    );
}

ULongWriter FileWriterBuilder::createULongWriter(std::string name) {
    return ULongWriter(
        std::make_unique<FileNumWriter<uint64_t>>((dir / name).string().c_str()), buffer_size
    );
}
DoubleWriter FileWriterBuilder::createDoubleWriter(std::string name) {
    return DoubleWriter(
        std::make_unique<FileNumWriter<double>>((dir / name).string().c_str()), buffer_size
    );
}

std::unique_ptr<StringWriter> FileWriterBuilder::createStringWriter(std::string name) {
    return std::make_unique<FileStringWriter>(dir / name);
}

void FileWriterBuilder::writeVersion(std::string version) {
    std::ofstream f((dir / "version").c_str());
    f << version << std::endl;
}

void FileWriterBuilder::deleteWriter(std::string name) { std::filesystem::remove(dir / name); }

FileReaderBuilder::FileReaderBuilder(std::string _dir, uint64_t buffer_size, uint64_t read_size)
    : dir(_dir)
    , buffer_size(buffer_size)
    , read_size(read_size) {

    if (!std::filesystem::exists(dir)) {
        throw std::invalid_argument(std::string("Missing directory: ") + _dir);
    }
}

UIntReader FileReaderBuilder::openUIntReader(std::string name) {
    return UIntReader(
        std::make_unique<FileNumReader<uint32_t>>((dir / name).string().c_str()),
        buffer_size,
        read_size
    );
}

FloatReader FileReaderBuilder::openFloatReader(std::string name) {
    return FloatReader(
        std::make_unique<FileNumReader<float>>((dir / name).string().c_str()),
        buffer_size,
        read_size
    );
}

ULongReader FileReaderBuilder::openULongReader(std::string name) {
    return ULongReader(
        std::make_unique<FileNumReader<uint64_t>>((dir / name).string().c_str()),
        buffer_size,
        read_size
    );
}
DoubleReader FileReaderBuilder::openDoubleReader(std::string name) {
    return DoubleReader(
        std::make_unique<FileNumReader<double>>((dir / name).string().c_str()),
        buffer_size,
        read_size
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
    if (!in) throw std::runtime_error(std::string("Could not open file: ") + path.string());

    while (std::getline(in, line)) {
        ret.push_back(line);
    }
    return ret;
}

} // end namespace BPCells
