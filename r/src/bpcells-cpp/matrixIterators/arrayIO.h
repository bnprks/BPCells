// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include "../utils/filesystem_compat.h"

#include "arrayIO/binaryfile.h"
#include "arrayIO/hdf5.h"

#include "PackedMatrix.h"

namespace BPCells {

PackedMatrix
openPackedMatrixH5(std::string file_path, std::string group_path, uint32_t buffer_size = 8192) {
    HighFive::SilenceHDF5 s;
    if (group_path == "") group_path = "/";

    HighFive::Group group(openH5ForReading(file_path).getGroup(group_path));

    std::string version;
    group.getAttribute("version").read(version);
    if (version != "v1-packed")
        throw std::runtime_error("HDF5 group does not have correct version attribute (v1-packed)");

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

PackedMatrixWriter createPackedMatrixH5(
    std::string file_path,
    std::string group_path,
    uint32_t buffer_size = 8192,
    uint32_t chunk_size = 1024
) {
    HighFive::SilenceHDF5 s;
    if (group_path == "") group_path = "/";

    std_fs::path path(file_path);
    if (path.has_parent_path() && !std_fs::exists(path.parent_path())) {
        std_fs::create_directories(path.parent_path());
    }

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

    std::string version("v1-packed");
    group.createAttribute<std::string>("version", HighFive::DataSpace::From(version))
        .write(version);

    return PackedMatrix(
        std::make_unique<ZH5UIntWriter>(group, "val_data", buffer_size),
        std::make_unique<ZH5UIntWriter>(group, "val_idx", buffer_size),
        std::make_unique<ZH5UIntWriter>(group, "row_data", buffer_size),
        std::make_unique<ZH5UIntWriter>(group, "row_starts", buffer_size),
        std::make_unique<ZH5UIntWriter>(group, "row_idx", buffer_size),
        std::make_unique<ZH5UIntWriter>(group, "col_ptr", buffer_size),
        std::make_unique<ZH5UIntWriter>(group, "row_count", buffer_size)
    );
}

} // end namespace BPCells