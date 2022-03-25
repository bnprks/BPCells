#pragma once

#include <stdint.h>
#include <map>
#include <vector>
#include <string>

#include "array_interfaces.h"

namespace BPCells {
    
template<class T>
class VecNumWriter final : public BulkNumWriter<T> {
private:
    std::vector<T> &vec;
public:
    VecNumWriter(std::vector<T> &vec) : vec(vec) {}
    
    uint32_t write(T *in, uint32_t count) override {
        size_t initial_size = vec.size();
        vec.resize(vec.size() + count);
        std::memmove(vec.data() + initial_size, in, sizeof(uint32_t) * count);
        return count;
    }
};

using VecUIntWriter = VecNumWriter<uint32_t>;

template<class T>
class VecNumReader : public BulkNumReader<T> {
private:
    const T *vec;
    uint32_t capacity;
    uint32_t pos = 0;
public:
    VecNumReader(const T *vec, std::size_t capacity) : vec(vec), capacity(capacity) {}

    // Return total number of integers in the reader
    uint32_t size() const override {return capacity; }

    // Change the next load to start at index pos
    void seek(uint32_t new_pos) override {pos = new_pos; }

    // Copy up to `count` integers into `out`, returning the actual number copied.
    // Will always load >0 unless there is no more input
    uint32_t load(T *out, uint32_t count) override {
        std::memmove(out, vec+pos, sizeof(T)*count);
        pos += count;
        return count;
    }
};

using VecUIntReader = VecNumReader<uint32_t>;

class VecStringWriter : public StringWriter {
private:
    std::vector<std::string> &data;
public:
    VecStringWriter(std::vector<std::string> &data);
    void write(const StringReader &reader) override;
};

class VecReaderWriterBuilder : public WriterBuilder, public ReaderBuilder {
protected:
    std::map<std::string, std::vector<uint32_t>> int_vecs;
    std::map<std::string, std::vector<float>> float_vecs;
    std::map<std::string, std::vector<uint64_t>> long_vecs;
    std::map<std::string, std::vector<double>> double_vecs;
    std::map<std::string, std::vector<std::string>> string_vecs;
    std::string version;
    uint32_t chunk_size;
public:
    VecReaderWriterBuilder(uint32_t chunk_size = 1024);
    UIntWriter createUIntWriter(std::string name) override;
    ULongWriter createULongWriter(std::string name) override;
    FloatWriter createFloatWriter(std::string name) override;
    DoubleWriter createDoubleWriter(std::string name) override;

    std::unique_ptr<StringWriter> createStringWriter(std::string name) override;
    void writeVersion(std::string version) override;

    UIntReader openUIntReader(std::string name) override;
    ULongReader openULongReader(std::string name) override;
    FloatReader openFloatReader(std::string name) override;
    DoubleReader openDoubleReader(std::string name) override;
    
    std::unique_ptr<StringReader> openStringReader(std::string name) override;
    std::string readVersion() override;

    std::map<std::string, std::vector<uint32_t>>& getIntVecs(std::string name);
    std::map<std::string, std::vector<float>>& getFloatVecs(std::string name);
    std::map<std::string, std::vector<uint64_t>>& getLongVecs(std::string name);
    std::map<std::string, std::vector<double>>& getDoubleVecs(std::string name);
    std::map<std::string, std::vector<std::string>>& getStringVec(std::string name);
};

} // end namespace BPCells