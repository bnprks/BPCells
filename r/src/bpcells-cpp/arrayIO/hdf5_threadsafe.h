// Copyright 2025 BPCells contributors
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <H5pubconf.h>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5PropertyList.hpp>

#ifndef H5_HAVE_THREADSAFE
#include <mutex>
#endif

// HDF5 Thread safety:
// According to the HDF5 docs, no HDF5 operations are thread safe by default:
//   https://support.hdfgroup.org/documentation/hdf5/latest/thread-safe-lib.html
// Even multi-threading on Linux seems to not have caused crashes, but we've gotten
// reports of crashes on Macs: https://github.com/bnprks/BPCells/issues/250

// This header provides a global lock and wrapper classes to help preserve thread
// safety in code that uses HDF5

namespace BPCells {

#ifndef H5_HAVE_THREADSAFE
std::recursive_mutex &hdf5_global_lock();
#endif

// Enforce thread safety and silence error messages for the duration of a scope
class H5ScopeGuard {
  private:
    HighFive::SilenceHDF5 s;
#ifndef H5_HAVE_THREADSAFE
    std::lock_guard<std::recursive_mutex> g;

  public:
    H5ScopeGuard() : g(hdf5_global_lock()) {}
#endif
};

template <typename T> class H5DataSet {
  protected:
    HighFive::DataSet dataset;
    HighFive::DataType datatype;

    static HighFive::DataSet copyH5DataSetSafe(const HighFive::DataSet &d) {
        H5ScopeGuard guard;
        return d;
    }

    template <class T> static HighFive::DataType createH5DataTypeSafe() {
        H5ScopeGuard guard;
        return HighFive::create_datatype<T>();
    }

  public:
    H5DataSet(HighFive::DataSet &&d) : dataset(std::move(d)), datatype(createH5DataTypeSafe<T>()) {}
    H5DataSet(const H5DataSet &other)
        : dataset(copyH5DataSetSafe(other.dataset))
        , datatype(createH5DataTypeSafe<T>()) {}
    H5DataSet &operator=(const H5DataSet &other) {
        H5ScopeGuard guard;
        dataset = other.dataset;
        datatype = other.datatype;
        return *this;
    }
    H5DataSet(H5DataSet &&other) = default;
    H5DataSet &operator=(H5DataSet &&other) = default;

    std::string getAttribute(const std::string &attribute_name) {
        H5ScopeGuard guard;
        std::string result;
        dataset.getAttribute(attribute_name).read(result);
        return result;
    }

    void setAttribute(const std::string &attribute_name, const std::string &value) {
        H5ScopeGuard guard;
        if (dataset.hasAttribute(attribute_name)) {
            dataset.getAttribute(attribute_name).write(value);
        } else {
            dataset.createAttribute<std::string>(attribute_name, value);
        }
    }
};

// Safe wrapper for a 1D HDF5 DataSet
template <typename T> class H5DataSet1D : H5DataSet<T> {
  public:
    using H5DataSet::H5DataSet;

    // Append to a 1D dataset
    void append(T *in, uint64_t count) {
        H5ScopeGuard guard;
        uint64_t cur_size = dataset.getDimensions()[0];
        dataset.resize({cur_size + count});
        dataset.select({cur_size}, {count}).write_raw(in, datatype);
    }

    
    // Read from a 1D dataset
    void load(uint64_t pos, T *out, uint64_t count) {
        H5ScopeGuard guard;
        dataset.select({pos}, {count}).read_raw(out, datatype);
    }
    
    // Read into a vector (with special support for reading variable/fixed-length strings)
    void load(uint64_t pos, std::vector<T> &out, uint64_t count);
    
    // Write to a 1D dataset. Should cover the full size
    void write(const std::vector<T> &data);

    size_t getDimension() const {
        H5ScopeGuard guard;
        return dataset.getDimensions()[0];
    }
};

template <typename T> class H5DataSet2D : H5DataSet<T> {
  public:
    using H5DataSet::H5DataSet;

    void write_row(size_t row, const std::vector<T> &data) {
        H5ScopeGuard guard;
        dataset.select({row, 0}, {1, data.size()}).write_raw(data.data(), datatype);
    }

    void write_col(size_t col, const std::vector<T> &data) {
        H5ScopeGuard guard;
        dataset.select({0, row}, {data.size(), 1}).write_raw(data.data(), datatype);
    }

    void read_partial_row(size_t row, size_t col, size_t count, T* out) {
        H5ScopeGuard guard;
        dataset.select({row, col}, {1, count}).read_raw(out, datatype);
    }
};

class H5Group {
  public:
    enum class OpenMode {
        // Read existing group
        Read,
        // Read/Write access for existing group
        Write,
        // Create a new group and open with Read/Write access
        Create,
        // Open/create a group with Read/Write access
        WriteOrCreate
    };

  private:
    HighFive::Group group;

    static HighFive::Group copyH5GroupSafe(const HighFive::Group &g) {
        H5ScopeGuard guard;
        return g;
    }

    // Try to open a file for read-write, then fall back to read only if needed.
    // If we first open a file ReadOnly, it prevents future opening with ReadWrite
    // (bad if we want to read + write the same file).
    // This retry makes it possible to still open a file if it's read-only though.
    static HighFive::File openH5File(const std::string &file_path, OpenMode mode) {
        H5ScopeGuard guard;

        switch (mode) {
        case OpenMode::Read:
            try {
                return HighFive::File(file_path, HighFive::File::ReadWrite);
            } catch (const HighFive::FileException &f) {
                return HighFive::File(file_path, HighFive::File::ReadOnly);
            }
        case OpenMode::Write:
            return HighFive::File(file_path, HighFive::File::ReadWrite);
        case OpenMode::Create: // Create refers to the group, so it's ok to have an existing file
        case OpenMode::WriteOrCreate:
            return HighFive::File(file_path, HighFive::File::OpenOrCreate);
        }
    }

    static HighFive::Group
    openH5Group(const std::string &file_path, const std::string &group_path, OpenMode mode) {
        H5ScopeGuard guard;

        HighFive::File f(openH5File(file_path, mode));

        switch (mode) {
        case OpenMode::Read:
        case OpenMode::Write:
            return f.getGroup(group_path);
        case OpenMode::Create:
            if (f.exist(group_path)) {
                HighFive::Group g(f.getGroup(group_path));
                if (g.getNumberObjects() != 0)
                    throw std::runtime_error("Requested hdf5 group not empty");
                return g;

            } else {
                return f.createGroup(group_path);
            }
        case OpenMode::WriteOrCreate:
            if (f.exist(group_path)) return f.getGroup(group_path);
            else return f.createGroup(group_path);
        }
    }

    // Create an H5 DataSet or overwrite existing
    // Assume that if any dim is 0, it is meant to be expandable
    // If chunk_size is 0, don't add chunking
    template <class T>
    HighFive::DataSet createH5DataSet(
        const std::string &dataset_path,
        const std::vector<size_t> &dims,
        uint64_t chunk_size,
        uint32_t gzip_level
    ) {
        H5ScopeGuard guard;

        std::vector<size_t> dim;
        std::vector<size_t> max_dim;
        std::vector<hsize_t> chunking;
        for (const size_t &x : dims) {
            dim.push_back(x);
            if (x == 0) {
                max_dim.push_back(HighFive::DataSpace::UNLIMITED);
            } else {
                max_dim.push_back(x);
            }
            chunking.push_back(std::min(chunk_size, max_dim.back()));
        }
        HighFive::DataSpace dataspace(dim, max_dim);

        HighFive::DataSetCreateProps props;
        if (chunk_size > 0) {
            props.add(HighFive::Chunking(chunking));
        }
        if (gzip_level > 0) {
            if constexpr (std::is_arithmetic_v<T>) props.add(HighFive::Shuffle());
            props.add(HighFive::Deflate(gzip_level));
        }

        if (group.exist(dataset_path)) group.unlink(dataset_path);

        return group.createDataSet<T>(dataset_path, dataspace, props);
    }

    H5Group(HighFive::Group &&g) : group(std::move(g)) {}

  public:
    H5Group(const H5Group &other) : group(copyH5GroupSafe(other.group)) {}
    H5Group &operator=(const H5Group &other) {
        group = copyH5GroupSafe(other.group);
        return *this;
    }
    H5Group(H5Group &&other) = default;
    H5Group &operator=(H5Group &&other) = default;

    H5Group(const std::string &file_path, const std::string &group_path, OpenMode mode)
        : group(openH5Group(file_path, group_path, mode)) {}

    bool exist(const std::string &group_path) const {
        H5ScopeGuard guard;
        return group.exist(group_path);
    }
    void unlink(const std::string &node_name) const {
        H5ScopeGuard guard;
        group.unlink(node_name);
    }

    /**
     * Create 1D dataset or overwrite existing
     * @param dim Size of dataset (or 0 for unlimited size)
     * @param chunk_size Chunk size to use (or 0 for no chunking on fixed-size dataset)
     * @param gzip_level Gzip compression level (or 0 for no compression)
     */
    template <class T>
    H5DataSet1D<T> createDataSet1D(
        const std::string &dataset_name, uint64_t dim, uint64_t chunk_size, uint32_t gzip_level
    ) {
        H5ScopeGuard guard;
        return H5DataSet1D<T>{createH5DataSet<T>(dataset_name, {dim}, chunk_size, gzip_level)};
    }

    /**
     * Create 2D dataset or overwrite existing
     * @param nrow Number of rows (or 0 for unlimited size)
     * @param ncol Number of columns (or 0 for unlimited size)
     * @param chunk_size Chunk size to use (or 0 for no chunking on fixed-size dataset)
     * @param gzip_level Gzip compression level (or 0 for no compression)
     */
    template <class T>
    H5DataSet2D<T> createDataSet2D(
        const std::string &dataset_name,
        uint64_t nrow,
        uint64_t ncol,
        uint64_t chunk_size,
        uint32_t gzip_level
    ) {
        H5ScopeGuard guard;
        return H5DataSet1D<T>{createH5DataSet<T>(dataset_name, {nrow, ncol}, chunk_size, gzip_level)};
    }

    template <class T> H5DataSet1D<T> openDataSet1D(const std::string &dataset_name) const {
        H5ScopeGuard guard;
        HighFive::DataSet ds(group.getDataSet(dataset_name));
        if (ds.getDimensions().size() != 1) {
            throw std::runtime_error(std::string("Expected 1D HDF5 dataset: ") + dataset_name);
        }
        return H5DataSet1D<T>{std::move(ds)};
    }

    template <class T> H5DataSet2D<T> openDataSet2D(const std::string &dataset_name) const {
        H5ScopeGuard guard;
        HighFive::DataSet ds(group.getDataSet(dataset_name));
        if (ds.getDimensions().size() != 2) {
            throw std::runtime_error(std::string("Expected 2D HDF5 dataset: ") + dataset_name);
        }
        return H5DataSet2D<T>{std::move(ds)};
    }

    std::string getAttribute(const std::string &attribute_name) {
        H5ScopeGuard guard;
        std::string result;
        group.getAttribute(attribute_name).read(result);
        return result;
    }

    void setAttribute(const std::string &attribute_name, const std::string &value) {
        H5ScopeGuard guard;
        if (group.hasAttribute(attribute_name)) {
            group.getAttribute(attribute_name).write(value);
        } else {
            group.createAttribute<std::string>(attribute_name, HighFive::DataSpace::From(value))
                .write(value);
        }
    }
};

} // namespace BPCells