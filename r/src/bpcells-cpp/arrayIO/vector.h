// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <map>
#include <stdint.h>
#include <string>
#include <vector>

#include "array_interfaces.h"

namespace BPCells {

template <class T> class VecNumWriter final : public BulkNumWriter<T> {
  private:
    std::vector<T> &vec;

  public:
    VecNumWriter(std::vector<T> &vec) : vec(vec) {}

    uint64_t write(T *in, uint64_t count) override {
        size_t initial_size = vec.size();
        vec.resize(vec.size() + count);
        std::memmove(vec.data() + initial_size, in, sizeof(T) * count);
        return count;
    }
};

using VecUIntWriter = VecNumWriter<uint32_t>;

template <class T> class PointerNumWriter final : public BulkNumWriter<T> {
  private:
    T *vec;
    uint64_t capacity;
    uint64_t pos = 0;
  public:
    // Capacity is given in units of sizeof(T), e.g. 4 bytes for uint32_t
    PointerNumWriter(T *vec, uint64_t capacity) : vec(vec), capacity(capacity) {}

    uint64_t write(T *in, uint64_t count) override {
        if (pos + count >= capacity) {
          throw std::runtime_error("PointerNumWriter: exceeded destination memory capacity");
        }
        std::memmove(vec + pos, in, sizeof(T)*count);
        pos += count;
        return count;
    }
};

template <class T> class VecNumReader : public BulkNumReader<T> {
  private:
    const T *vec;
    uint64_t capacity;
    uint64_t pos = 0;

  public:
    VecNumReader(const T *vec, std::size_t capacity) : vec(vec), capacity(capacity) {}

    // Return total number of integers in the reader
    uint64_t size() const override { return capacity; }

    // Change the next load to start at index pos
    void seek(uint64_t new_pos) override { pos = new_pos; }

    // Copy up to `count` integers into `out`, returning the actual number copied.
    // Will always load >0 unless there is no more input
    uint64_t load(T *out, uint64_t count) override {
        uint64_t load_size = std::min<uint64_t>(capacity - pos, count);
        std::memmove(out, vec + pos, sizeof(T) * load_size);
        pos += load_size;
        return load_size;
    }
};

using VecUIntReader = VecNumReader<uint32_t>;

class VecStringWriter : public StringWriter {
  private:
    std::vector<std::string> &data;

  public:
    VecStringWriter(std::vector<std::string> &data);
    void write(StringReader &reader) override;
};

class VecReaderWriterBuilder : public WriterBuilder, public ReaderBuilder {
  protected:
    std::map<std::string, std::vector<uint32_t>> int_vecs;
    std::map<std::string, std::vector<float>> float_vecs;
    std::map<std::string, std::vector<uint64_t>> long_vecs;
    std::map<std::string, std::vector<double>> double_vecs;
    std::map<std::string, std::vector<std::string>> string_vecs;
    std::string version;
    uint64_t chunk_size;

  public:
    VecReaderWriterBuilder(uint64_t chunk_size = 1024);
    UIntWriter createUIntWriter(std::string name) override;
    ULongWriter createULongWriter(std::string name) override;
    FloatWriter createFloatWriter(std::string name) override;
    DoubleWriter createDoubleWriter(std::string name) override;

    std::unique_ptr<StringWriter> createStringWriter(std::string name) override;
    void writeVersion(std::string version) override;
    // Delete all writers with a given name
    void deleteWriter(std::string name) override;

    UIntReader openUIntReader(std::string name) override;
    ULongReader openULongReader(std::string name) override;
    FloatReader openFloatReader(std::string name) override;
    DoubleReader openDoubleReader(std::string name) override;

    std::unique_ptr<StringReader> openStringReader(std::string name) override;
    std::string readVersion() override;

    std::map<std::string, std::vector<uint32_t>> &getIntVecs();
    std::map<std::string, std::vector<float>> &getFloatVecs();
    std::map<std::string, std::vector<uint64_t>> &getLongVecs();
    std::map<std::string, std::vector<double>> &getDoubleVecs();
    std::map<std::string, std::vector<std::string>> &getStringVecs();
};

} // end namespace BPCells