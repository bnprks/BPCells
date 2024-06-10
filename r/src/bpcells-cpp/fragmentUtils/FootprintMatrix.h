// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <atomic>

#ifndef RCPP_EIGEN
#include <Eigen/Dense>
#else
#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <RcppEigen.h>
#endif

#include "../arrayIO/array_interfaces.h"
#include "../fragmentIterators/FragmentIterator.h"

namespace BPCells {

// Make a matrix for footprinting signal around motifs or TSS sites
// Arguments:
//    frags: Input fragment loader object
//    chr, center, strand: chromosome, center bp and strand information for footprint sites
//    flank_width: Number of bp to extend footprinting (symmetrically applied
//        both upstream and downstream)
//    cell_groups, cell_weights: Which group (output row) each cell should count
//        in, and what value should be added to the matrix from one fragment in
//        a cell.
//    chr_levels: List of chromosome levels
// Return:
//   Eigen matrix of dimensions max(cell_groups) x 2 * flank_width + 1.
//   Each entry[i,j] is the sum of the cell_weight[c] for each cell where
//   cell_groups[c] == i, summed across all fragments where the start/end coordinate
//   is: (j <= flank_width) => flank_width+1-j bases upstream; (j >= flank_width) j - flank_width
//   bases downstream;
Eigen::MatrixXd footprintMatrix(
    FragmentLoader &frags,
    const std::vector<uint32_t> &chr,
    const std::vector<uint32_t> &center,
    const std::vector<int32_t> &strand,
    uint32_t flank_width,
    std::unique_ptr<StringReader> &&chr_levels,
    const std::vector<uint32_t> &cell_groups,
    const std::vector<double> &cell_weights,
    std::atomic<bool> *user_interrupt
);

} // namespace BPCells