// Copyright 2021 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <algorithm>

#include "../arrayIO/array_interfaces.h"
#include "../fragmentIterators/FragmentIterator.h"
#include "MatrixAccumulators.h"
#include "MatrixIterator.h"

namespace BPCells {

// Output cell x peak matrix (rows = cell_id, col = peak_id)
// Regions are given as half-open format
// Peaks can overlap, and columns are ordered by (end, start) coordinate

// Notes on MODE template:
// - The template is in order to easily swap out the ~30 lines that actually register
//   overlap counts without having to copy-paste a bunch of code
// - The cpp file still contains the method definitions, and explicit instatiations
//   of the valid MODE options
// - Valid MODE options:
//    - 0: Count insertions per peak. One fragment can count twice if it lands
//      in same peak twice
//    - 1: Count fragments per peak. Same as MODE == 0 but fragment cannot count
//      twice per peak
//    - 2: Count fragment overlaps. Same as MODE == 1 but fragments that fully
//      surround a peak will also be counted in that peak
class PeakMatrix : public MatrixLoader<uint32_t> {
  private:
    class Peak {
      public:
        uint32_t chr, start, end;
    };
    const uint32_t mode;

    std::unique_ptr<FragmentLoader> frags;
    std::unique_ptr<StringReader> chr_levels;
    MatrixAccumulator<uint32_t> accumulator;
    std::vector<uint32_t> end_sorted_lookup; // end_sorted_lookup[i] gives the index of end-sorted
                                             // peak i in the start-sorted list
    std::vector<Peak> sorted_peaks;
    std::vector<Peak> active_peaks;
    uint32_t next_completed_peak = 0;
    uint32_t current_output_peak = UINT32_MAX;
    uint32_t next_active_peak = 0;
    uint32_t n_peaks;

    std::string peak_name; // buffer to use to store the peak name

    void loadFragments();

  public:
    // Note: It's the caller's responsibility to make sure that
    // the FragmentLoader will not be deleted while this PeakMatrix is still alive
    // Arguments:
    // frags - fragments loader with the input fragments for the peak matrix
    // chr - list of chrIDs (following the chrIDs in frags)
    // start, end - list of start + end coordinates for the peaks (start inclusive, end exclusive)
    // chr_levels - list of expected chr levels, for safety checking that peaks are coming from the
    //    correct chromosomes
    PeakMatrix(
        std::unique_ptr<FragmentLoader> &&frags,
        const std::vector<uint32_t> &chr,
        const std::vector<uint32_t> &start,
        const std::vector<uint32_t> &end,
        std::unique_ptr<StringReader> &&chr_levels,
        uint32_t mode
    );

    uint32_t rows() const override;
    uint32_t cols() const override;

    const char *rowNames(uint32_t row) override;
    const char *colNames(uint32_t col) override;

    void restart() override;
    void seekCol(uint32_t col) override;

    bool nextCol() override;

    uint32_t currentCol() const override;

    bool load() override;

    uint32_t capacity() const override;

    uint32_t *rowData() override;
    uint32_t *valData() override;
};

} // end namespace BPCells