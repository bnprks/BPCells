// Copyright 2025 BPCells contributors
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "hdf5_threadsafe.h"

namespace BPCells {

std::recursive_mutex &hdf5_global_lock() {
    static std::recursive_mutex mut;
    return mut;
}

template <class T> void H5DataSet1D<T>::load(uint64_t pos, std::vector<T> &out, uint64_t count) {
    H5ScopeGuard guard;
    dataset.select({pos}, {count}).read(out);
}

// Specialized handling for fixed-sized and variable-sized string types
template <>
void H5DataSet1D<std::string>::load(uint64_t pos, std::vector<std::string> &out, uint64_t count) {
    H5ScopeGuard guard;

    HighFive::DataType type = dataset.getDataType();

    if (type.isVariableStr()) {
        // Workaround for HighFive bug: don't try reading an empty string vector
        if (count > 0) {
            dataset.select({pos}, {count}).read(out);
        } else {
            out.resize(0);
        }
    } else {
        uint64_t bytes = type.getSize();
        std::vector<char> char_data(bytes * count);
        dataset.select({pos}, {count}).read_raw(char_data.data(), type);
        out.resize(count);
        for (uint64_t i = 0; i < count; i++) {
            out[i] = std::string(char_data.data() + bytes * i, char_data.data() + bytes * (i + 1));
        }
    }
}

template <class T> void H5DataSet1D<T>::write(const std::vector<T> &data) {
    H5ScopeGuard guard;
    if (getDimension() != data.size()) {
        throw std::runtime_error("H5Dataset1D::write() data size does not match dimension size.");
    }
    if (data.size() > 0) dataset.write(data);
}

} // namespace BPCells