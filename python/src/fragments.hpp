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
#include <vector>

#include <Eigen/Core>

namespace BPCells::py {

void import_10x_fragments(
    std::string input_10x,
    std::string output_bpcells,
    int shift_start,
    int shift_end,
    std::optional<std::vector<std::string>> keeper_cells
);

std::vector<uint32_t> echo_vec(std::vector<uint32_t> a);

std::vector<std::string> cell_names_fragments_dir(std::string input_bpcells);
std::vector<std::string> chr_names_fragments_dir(std::string input_bpcells);

Eigen::MatrixXi pseudobulk_coverage(
    std::string fragments_path,
    std::vector<std::string> chr,
    std::vector<uint32_t> start,
    std::vector<uint32_t> end,
    std::vector<int32_t> cell_groups,
    int bin_size
);

// Calculate an experimental CompressedSparseColumn matrix
void precalculate_pseudobulk_coverage(
    std::string fragments_path,
    std::string output_path,
    std::string tmp_path,
    std::vector<std::string> chr,
    std::vector<uint32_t> chr_size,
    std::vector<int32_t> cell_groups,
    int bin_size,
    int threads
);

Eigen::MatrixXi query_precalculated_pseudobulk_coverage(
    std::string mat_path,
    std::vector<uint32_t> range_starts,
    uint32_t range_len
);

} // namespace BPCells::py