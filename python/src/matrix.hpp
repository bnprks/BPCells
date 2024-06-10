// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <optional>
#include <string>
#include <utility>

#include <Eigen/SparseCore>

#include "bpcells-cpp/arrayIO/vector.h"

namespace BPCells::py {


void write_matrix_dir_from_memory(const Eigen::SparseMatrix<uint32_t> in, std::string out_path, bool row_major = false);

void write_matrix_dir_from_concat(std::vector<std::string> in_paths, std::string out_path, bool concat_cols);

void write_matrix_dir_from_h5ad(std::string h5ad_path, std::string out_path, std::string group);

std::vector<Eigen::SparseMatrix<uint32_t>> load_matrix_dir_subset(
    std::string matrix_path,
    std::optional<std::vector<uint32_t>> rows,
    std::vector<uint32_t> columns,
    uint32_t threads
);

std::tuple<uint32_t, uint32_t> dims_matrix_dir(std::string matrix_path);

VecReaderWriterBuilder load_matrix_dir_to_memory(std::string matrix_path);

std::vector<Eigen::SparseMatrix<uint32_t>> load_matrix_memory_subset(
    VecReaderWriterBuilder &rb,
    std::optional<std::vector<uint32_t>> rows,
    std::vector<uint32_t> columns,
    uint32_t threads
);

} // namespace BPCells::py