// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "vector.h"
#include <cstring>
namespace BPCells {

VecStringWriter::VecStringWriter(std::vector<std::string> &data) : data(data) {}
void VecStringWriter::write(StringReader &reader) {
    uint64_t i = 0;
    data.resize(0);
    while (true) {
        const char *s = reader.get(i);
        if (s == NULL) break;
        data.push_back(s);
        i++;
    }
}

VecReaderWriterBuilder::VecReaderWriterBuilder(uint64_t chunk_size) : chunk_size(chunk_size) {}

// UIntWriter VecReaderWriterBuilder::createUIntWriter(std::string name) {
//     int_vecs[name] = std::vector<uint32_t>();
//     return UIntWriter(
//         std::make_unique<VecUIntWriter>(int_vecs.at(name)),
//         chunk_size
//     );
// }

UIntWriter VecReaderWriterBuilder::createUIntWriter(std::string name) {
    int_vecs[name] = std::vector<uint32_t>();
    return UIntWriter(std::make_unique<VecUIntWriter>(int_vecs.at(name)), chunk_size);
}

FloatWriter VecReaderWriterBuilder::createFloatWriter(std::string name) {
    float_vecs[name] = std::vector<float>();
    return FloatWriter(std::make_unique<VecNumWriter<float>>(float_vecs.at(name)), chunk_size);
}

ULongWriter VecReaderWriterBuilder::createULongWriter(std::string name) {
    long_vecs[name] = std::vector<uint64_t>();
    return ULongWriter(std::make_unique<VecNumWriter<uint64_t>>(long_vecs.at(name)), chunk_size);
}

DoubleWriter VecReaderWriterBuilder::createDoubleWriter(std::string name) {
    double_vecs[name] = std::vector<double>();
    return DoubleWriter(std::make_unique<VecNumWriter<double>>(double_vecs.at(name)), chunk_size);
}

std::unique_ptr<StringWriter> VecReaderWriterBuilder::createStringWriter(std::string name) {
    string_vecs.emplace(name, std::vector<std::string>());
    return std::make_unique<VecStringWriter>(string_vecs.at(name));
}
void VecReaderWriterBuilder::writeVersion(std::string version) { this->version = version; }

void VecReaderWriterBuilder::deleteWriter(std::string name) {
    this->int_vecs.erase(name);
    this->float_vecs.erase(name);
    this->long_vecs.erase(name);
    this->double_vecs.erase(name);
    this->string_vecs.erase(name);
}

UIntReader VecReaderWriterBuilder::openUIntReader(std::string name) {
    std::vector<uint32_t> &v = int_vecs.at(name);
    return UIntReader(std::make_unique<VecUIntReader>(v.data(), v.size()), chunk_size, chunk_size);
}

FloatReader VecReaderWriterBuilder::openFloatReader(std::string name) {
    std::vector<float> &v = float_vecs.at(name);
    return FloatReader(
        std::make_unique<VecNumReader<float>>(v.data(), v.size()), chunk_size, chunk_size
    );
}

ULongReader VecReaderWriterBuilder::openULongReader(std::string name) {
    std::vector<uint64_t> &v = long_vecs.at(name);
    return ULongReader(
        std::make_unique<VecNumReader<uint64_t>>(v.data(), v.size()), chunk_size, chunk_size
    );
}

DoubleReader VecReaderWriterBuilder::openDoubleReader(std::string name) {
    std::vector<double> &v = double_vecs.at(name);
    return DoubleReader(
        std::make_unique<VecNumReader<double>>(v.data(), v.size()), chunk_size, chunk_size
    );
}

std::unique_ptr<StringReader> VecReaderWriterBuilder::openStringReader(std::string name) {
    std::vector<std::string> &v = string_vecs.at(name);
    return std::make_unique<VecStringReader>(v);
}
std::string VecReaderWriterBuilder::readVersion() { return version; }

std::map<std::string, std::vector<uint32_t>> &VecReaderWriterBuilder::getIntVecs() {
    return int_vecs;
}
std::map<std::string, std::vector<float>> &VecReaderWriterBuilder::getFloatVecs() {
    return float_vecs;
}
std::map<std::string, std::vector<uint64_t>> &VecReaderWriterBuilder::getLongVecs() {
    return long_vecs;
}
std::map<std::string, std::vector<double>> &VecReaderWriterBuilder::getDoubleVecs() {
    return double_vecs;
}
std::map<std::string, std::vector<std::string>> &VecReaderWriterBuilder::getStringVecs() {
    return string_vecs;
}

} // end namespace BPCells
