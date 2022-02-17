#include "vector.h"
#include <cstring>
namespace BPCells {


VecUIntWriter::VecUIntWriter(std::vector<uint32_t> &vec): vec(vec) {}

uint32_t VecUIntWriter::write(uint32_t *in, uint32_t count) {
    size_t initial_size = vec.size();
    vec.resize(vec.size() + count);
    std::memmove(vec.data() + initial_size, in, sizeof(uint32_t) * count);
    return count;
}

VecUIntReader::VecUIntReader(const uint32_t *vec, std::size_t capacity) :
    vec(vec), capacity(capacity) {}

uint32_t VecUIntReader::size() const {return capacity;}

void VecUIntReader::seek(uint32_t new_pos) {pos = new_pos;}

uint32_t VecUIntReader::load(uint32_t *out, uint32_t count) {
    std::memmove(out, vec+pos, sizeof(uint32_t)*count);
    pos += count;
    return count;
};

VecStringWriter::VecStringWriter(std::vector<std::string> &data) : data(data) {}
void VecStringWriter::write(const StringReader &reader) {
    uint32_t i = 0;
    data.resize(0);
    while (true) {
        const char* s = reader.get(i);
        if (s == NULL) break;
        data.push_back(s);
        i++;
    }
}

VecReaderWriterBuilder::VecReaderWriterBuilder(uint32_t chunk_size) : chunk_size(chunk_size) {}

UIntWriter VecReaderWriterBuilder::createUIntWriter(std::string name) {
    int_vecs[name] = std::vector<uint32_t>();
    return UIntWriter(
        std::make_unique<VecUIntWriter>(int_vecs.at(name)),
        chunk_size
    );
}
std::unique_ptr<StringWriter> VecReaderWriterBuilder::createStringWriter(std::string name) {
    string_vecs.emplace(name, std::vector<std::string>());
    return std::make_unique<VecStringWriter>(string_vecs.at(name));
}
void VecReaderWriterBuilder::writeVersion(std::string version) {this->version = version;} // Don't store version information

UIntReader VecReaderWriterBuilder::openUIntReader(std::string name) {
    std::vector<uint32_t> &v = int_vecs.at(name);
    return UIntReader(std::make_unique<VecUIntReader>(v.data(), v.size()), 1024, 1024);
}

std::unique_ptr<StringReader> VecReaderWriterBuilder::openStringReader(std::string name) {
    std::vector<std::string> &v = string_vecs.at(name);
    return std::make_unique<VecStringReader>(v);
}
std::string VecReaderWriterBuilder::readVersion() {return version;}

std::vector<uint32_t>& VecReaderWriterBuilder::getIntVec(std::string name) {return int_vecs.at(name);}
std::vector<std::string>& VecReaderWriterBuilder::getStringVec(std::string name) {return string_vecs.at(name);}
} // end namespace BPCells