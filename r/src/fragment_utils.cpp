// Copyright 2021 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include <sstream>

#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <Rcpp.h>
#include <RcppEigen.h>

#include "R_array_io.h"
#include "R_xptr_wrapper.h"
#include "R_interrupts.h"

#include "bpcells-cpp/fragmentIterators/CellSelect.h"
#include "bpcells-cpp/fragmentIterators/ChrSelect.h"
#include "bpcells-cpp/fragmentIterators/FragmentIterator.h"
#include "bpcells-cpp/fragmentIterators/LengthSelect.h"
#include "bpcells-cpp/fragmentIterators/MergeFragments.h"
#include "bpcells-cpp/fragmentIterators/RegionSelect.h"
#include "bpcells-cpp/fragmentIterators/Rename.h"
#include "bpcells-cpp/fragmentIterators/ShiftCoords.h"
#include "bpcells-cpp/fragmentUtils/BedgraphWriter.h"
#include "bpcells-cpp/fragmentUtils/FootprintMatrix.h"

#include "bpcells-cpp/matrixIterators/PeakMatrix.h"
#include "bpcells-cpp/matrixIterators/TileMatrix.h"
// #include "matrixIterators/PeakMatrix2.h"
// #include "matrixIterators/PeakMatrix3.h"
// #include "matrixIterators/PeakMatrix4.h"
// #include "matrixIterators/PeakMatrix5.h"
// #include "matrixIterators/PeakMatrix6.h"
// #include "matrixIterators/TileMatrix.h"
// #include "matrixIterators/TileMatrix2.h"

using namespace Rcpp;
using namespace BPCells;

// [[Rcpp::export]]
NumericVector scan_fragments_cpp(SEXP fragments) {
    FragmentIterator iter(take_unique_xptr<FragmentLoader>(fragments));
    uint64_t start_sum = 0;
    uint64_t end_sum = 0;
    uint64_t cell_sum = 0;
    uint64_t len = 0;

    while (iter.nextChr()) {
        while (iter.nextFrag()) {
            len++;
            start_sum += iter.start();
            end_sum += iter.end();
            cell_sum += iter.cell();
        }
    }
    return {
        (double)len,
        (double)start_sum,
        (double)end_sum,
        (double)cell_sum,
        (double)(((start_sum % 104729) + (end_sum % 104729)) % 104729)};
}

// [[Rcpp::export]]
StringVector get_tile_names_cpp(
    IntegerVector chr_id,
    IntegerVector start,
    IntegerVector end,
    IntegerVector tile_width,
    StringVector chr_levels
) {
    int64_t total_tiles = 0;
    for (int64_t i = 0; i < chr_id.length(); i++) {
        total_tiles += (end[i] - start[i] + tile_width[i] - 1) / tile_width[i];
    }
    StringVector ret(total_tiles);
    size_t out_idx = 0;
    for (int64_t i = 0; i < chr_id.length(); i++) {
        for (int64_t base = start[i]; base < end[i]; base += tile_width[i]) {
            std::stringstream out;
            out << (const char *)chr_levels[chr_id[i]] << ":" << (base) << "-"
                << std::min((int64_t)end[i], base + tile_width[i]);
            ret[out_idx++] = out.str();
        }
    }
    return ret;
}

// [[Rcpp::export]]
List get_tile_ranges_cpp(
    IntegerVector chr_id,
    IntegerVector start,
    IntegerVector end,
    IntegerVector tile_width,
    StringVector chr_levels,
    NumericVector selection
) {
    // Calculate starting index for each tile range
    std::vector<int64_t> starting_idx;
    int64_t total_tiles = 0;
    for (int64_t i = 0; i < chr_id.length(); i++) {
        starting_idx.push_back(total_tiles);
        total_tiles += (end[i] - start[i] + tile_width[i] - 1) / tile_width[i];
    }

    // Calculate coordinates for each selection
    std::vector<int32_t> ret_chr_id, ret_start, ret_end;
    for (int64_t i = 0; i < selection.length(); i++) {
        int64_t idx = (int64_t)selection[i];
        int64_t range = std::upper_bound(starting_idx.begin(), starting_idx.end(), idx) - 1 -
                        starting_idx.begin();

        ret_chr_id.push_back(chr_id[range]);
        ret_start.push_back(start[range] + tile_width[range] * (idx - starting_idx[range]));
        ret_end.push_back(std::min(ret_start.back() + tile_width[range], end[range]));
    }
    return List::create(
        Named("chr") = ret_chr_id, Named("start") = ret_start, Named("end") = ret_end
    );
}

// [[Rcpp::export]]
List subset_tiles_cpp(
    IntegerVector chr_id,
    IntegerVector start,
    IntegerVector end,
    IntegerVector tile_width,
    StringVector chr_levels,
    NumericVector selection
) {
    // Check that selection is sorted
    double prev = -1;
    for (double &x : selection) {
        if (x < 0 || x <= prev) {
            stop("subset_tiles_cpp: Selection must be non-negative and sorted in increasing order");
        }
    }

    // Calculate starting index for each tile range
    std::vector<int64_t> starting_idx;
    int64_t total_tiles = 0;
    for (int64_t i = 0; i < chr_id.length(); i++) {
        starting_idx.push_back(total_tiles);
        total_tiles += (end[i] - start[i] + tile_width[i] - 1) / tile_width[i];
    }

    // Calculate coordinates for each selection
    std::vector<int32_t> ret_chr_id, ret_start, ret_end, ret_width;
    
    int64_t prev_idx = -100; 
    int64_t prev_range = -100;
    for (int64_t i = 0; i < selection.length(); i++) {
        int64_t idx = (int64_t)selection[i];
        int64_t range = std::upper_bound(starting_idx.begin(), starting_idx.end(), idx) - 1 -
                        starting_idx.begin();
        if (idx - 1 == prev_idx && range == prev_range) {
            ret_end.back() += tile_width[range];
        } else {
            ret_chr_id.push_back(chr_id[range]);
            ret_start.push_back(start[range] + tile_width[range] * (idx - starting_idx[range]));
            ret_end.push_back(std::min(ret_start.back() + tile_width[range], end[range]));
            ret_width.push_back(tile_width[range]);
        }
        prev_idx = idx;
        prev_range = range;
    }
     return List::create(
        Named("chr_id") = ret_chr_id, Named("start") = ret_start, Named("end") = ret_end, Named("tile_width") = ret_width
    );
}

// Compute number of fragments like ArchR does for cell stats
//  [[Rcpp::export]]
List nucleosome_counts_cpp(SEXP fragments, uint32_t nuc_width = 147) {
    FragmentIterator iter(take_unique_xptr<FragmentLoader>(fragments));

    std::vector<int> subNuc(0);
    std::vector<int> monoNuc(0);
    std::vector<int> multiNuc(0);

    if (iter.cellCount() != -1) {
        subNuc.resize(iter.cellCount(), 0);
        monoNuc.resize(iter.cellCount(), 0);
        multiNuc.resize(iter.cellCount(), 0);
    }

    uint32_t count = 0;
    while (iter.nextChr()) {
        while (iter.nextFrag()) {
            if (count++ % (1 << 14) == 0) Rcpp::checkUserInterrupt();
            uint32_t cell = iter.cell();
            if (cell >= subNuc.size()) {
                subNuc.resize(cell + 1, 0);
                monoNuc.resize(cell + 1, 0);
                multiNuc.resize(cell + 1, 0);
            }
            uint32_t size = (iter.end() - iter.start()) / nuc_width;
            if (size == 0) subNuc[cell]++;
            else if (size == 1) monoNuc[cell]++;
            else multiNuc[cell]++;
        }
    }

    return List::create(
        Named("subNucleosomal") = subNuc,
        Named("monoNucleosomal") = monoNuc,
        Named("multiNucleosomal") = multiNuc
    );
}

// Compute fragment length distribution
// [[Rcpp::export]]
std::vector<int> fragment_lengths_cpp(SEXP fragments) {
    FragmentIterator iter(take_unique_xptr<FragmentLoader>(fragments));

    std::vector<int> lengths(0);

    while (iter.nextChr()) {
        while (iter.nextFrag()) {
            uint32_t width = iter.end() - iter.start();
            if (width >= lengths.size()) lengths.resize(width + 1);
            lengths[width] += 1;
        }
    }

    return lengths;
}

// [[Rcpp::export]]
Eigen::MatrixXd footprint_matrix_cpp(
    SEXP fragments,
    std::vector<uint32_t> chr,
    std::vector<uint32_t> center,
    std::vector<int32_t> strand,
    uint32_t flank_width,
    StringVector chr_levels,
    std::vector<uint32_t> cell_groups,
    std::vector<double> cell_weights
) {
    auto frags = take_unique_xptr<FragmentLoader>(fragments);
    return run_with_R_interrupt_check(
        &footprintMatrix,
        std::ref(*frags),
        std::cref(chr),
        std::cref(center),
        std::cref(strand),
        flank_width,
        std::make_unique<RcppStringReader>(chr_levels),
        std::cref(cell_groups),
        std::cref(cell_weights)
    );
}

// [[Rcpp::export]]
void write_insertion_bedgraph_cpp(
    SEXP fragments,
    std::vector<uint32_t> cell_groups,
    std::vector<std::string> output_paths,
    std::string mode_string
) {
    BedgraphInsertionMode mode;
    if (mode_string == "both") {
        mode = BedgraphInsertionMode::Both;
    } else if (mode_string == "start_only") {
        mode = BedgraphInsertionMode::StartOnly;
    } else if (mode_string == "end_only") {
        mode = BedgraphInsertionMode::EndOnly;
    } else {
        throw std::runtime_error("write_bedgraph_cpp: invalid mode found: " + mode_string);
    }

    auto frags = take_unique_xptr<FragmentLoader>(fragments);
    
    run_with_R_interrupt_check(
        &writeInsertionBedgraph,
        std::ref(*frags),
        std::cref(cell_groups),
        std::cref(output_paths),
        mode
    );
}


// [[Rcpp::export]]
SEXP iterate_peak_matrix_cpp(
    SEXP fragments,
    std::vector<uint32_t> chr,
    std::vector<uint32_t> start,
    std::vector<uint32_t> end,
    StringVector chr_levels,
    std::string mode
) {
    auto loader = take_unique_xptr<FragmentLoader>(fragments);
    std::unique_ptr<StringReader> chr_levels_reader =
        std::make_unique<RcppStringReader>(chr_levels);
    uint32_t mode_int = 0;
    if (mode == "insertions") {
        mode_int = 0;
    } else if (mode == "fragments") {
        mode_int = 1;
    } else if (mode == "overlaps") {
        mode_int = 2;
    } else {
        throw std::invalid_argument("mode must be one of insertions, fragments, or overlaps");
    }
    return make_unique_xptr<PeakMatrix>(
        std::move(loader), chr, start, end, std::move(chr_levels_reader), mode_int
    );
}

// [[Rcpp::export]]
SEXP iterate_tile_matrix_cpp(
    SEXP fragments,
    std::vector<uint32_t> chr,
    std::vector<uint32_t> start,
    std::vector<uint32_t> end,
    std::vector<uint32_t> width,
    StringVector chr_levels,
    std::string mode
) {
    bool count_fragments;
    if (mode == "insertions") {
        count_fragments = false;
    } else if (mode == "fragments") {
        count_fragments = true;
    } else {
        throw std::invalid_argument("mode must be one of insertions or fragments");
    }
    return make_unique_xptr<TileMatrix>(
        take_unique_xptr<FragmentLoader>(fragments),
        chr,
        start,
        end,
        width,
        std::make_unique<RcppStringReader>(chr_levels),
        count_fragments
    );
}

// [[Rcpp::export]]
SEXP iterate_shift_cpp(SEXP fragments, int32_t shift_start, int32_t shift_end) {
    return make_unique_xptr<ShiftCoords>(
        take_unique_xptr<FragmentLoader>(fragments), shift_start, shift_end
    );
}

// [[Rcpp::export]]
SEXP iterate_length_select_cpp(SEXP fragments, uint32_t min_len, uint32_t max_len) {
    return make_unique_xptr<LengthSelect>(
        take_unique_xptr<FragmentLoader>(fragments), min_len, max_len
    );
}

// [[Rcpp::export]]
SEXP iterate_chr_index_select_cpp(SEXP fragments, std::vector<uint32_t> chr_selection) {
    return make_unique_xptr<ChrIndexSelect>(
        take_unique_xptr<FragmentLoader>(fragments), chr_selection
    );
}

// [[Rcpp::export]]
SEXP iterate_chr_name_select_cpp(SEXP fragments, std::vector<std::string> chr_selection) {
    return make_unique_xptr<ChrNameSelect>(
        take_unique_xptr<FragmentLoader>(fragments), chr_selection
    );
}

// [[Rcpp::export]]
SEXP iterate_cell_index_select_cpp(SEXP fragments, std::vector<uint32_t> cell_selection) {
    return make_unique_xptr<CellIndexSelect>(
        take_unique_xptr<FragmentLoader>(fragments), cell_selection
    );
}

// [[Rcpp::export]]
SEXP iterate_cell_name_select_cpp(SEXP fragments, std::vector<std::string> cell_selection) {
    return make_unique_xptr<CellNameSelect>(
        take_unique_xptr<FragmentLoader>(fragments), cell_selection
    );
}

// [[Rcpp::export]]
SEXP iterate_cell_merge_cpp(
    SEXP fragments, std::vector<uint32_t> group_ids, StringVector group_names
) {
    return make_unique_xptr<CellMerge>(
        take_unique_xptr<FragmentLoader>(fragments),
        group_ids,
        std::make_unique<RcppStringReader>(group_names)
    );
}

// [[Rcpp::export]]
SEXP iterate_chr_rename_cpp(SEXP fragments, const StringVector &chr_names) {
    return make_unique_xptr<RenameChrs>(
        take_unique_xptr<FragmentLoader>(fragments), std::make_unique<RcppStringReader>(chr_names)
    );
}

// [[Rcpp::export]]
SEXP iterate_cell_rename_cpp(SEXP fragments, const StringVector &cell_names) {
    return make_unique_xptr<RenameCells>(
        take_unique_xptr<FragmentLoader>(fragments), std::make_unique<RcppStringReader>(cell_names)
    );
}

// [[Rcpp::export]]
SEXP iterate_cell_prefix_cpp(SEXP fragments, std::string &prefix) {
    return make_unique_xptr<PrefixCells>(take_unique_xptr<FragmentLoader>(fragments), prefix);
}

// [[Rcpp::export]]
SEXP iterate_region_select_cpp(
    SEXP fragments,
    std::vector<uint32_t> chr,
    std::vector<uint32_t> start,
    std::vector<uint32_t> end,
    StringVector chr_levels,
    bool invert_selection
) {
    return make_unique_xptr<RegionSelect>(
        take_unique_xptr<FragmentLoader>(fragments),
        chr,
        start,
        end,
        std::make_unique<RcppStringReader>(chr_levels),
        invert_selection
    );
}

// [[Rcpp::export]]
SEXP iterate_merge_fragments_cpp(SEXP fragments_list, std::vector<std::string> chr_order) {
    std::vector<std::unique_ptr<FragmentLoader>> fragments_vec;
    List l = fragments_list;
    for (uint32_t i = 0; i < l.size(); i++) {
        SEXP elem = l[i];
        fragments_vec.push_back(take_unique_xptr<FragmentLoader>(elem));
    }

    return make_unique_xptr<MergeFragments>(std::move(fragments_vec), chr_order);
}
