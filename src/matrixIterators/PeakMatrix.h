#pragma once

#include <algorithm>

#include "../arrayIO/array_interfaces.h"
#include "../bitpacking/simd_vec.h"
#include "../fragmentIterators/FragmentIterator.h"
#include "MatrixAccumulators.h"
#include "MatrixIterator.h"

namespace BPCells {

// Output cell x peak matrix (rows = cell_id, col = peak_id)
// Regions are given as half-open format
// Peaks can overlap, and columns are ordered by (end, start) coordinate
class PeakMatrix : public MatrixLoader<uint32_t> {
  private:
    class Peak {
      public:
        uint32_t chr, start, end;
    };

    FragmentLoader &frags;
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
        FragmentLoader &frags,
        const std::vector<uint32_t> &chr,
        const std::vector<uint32_t> &start,
        const std::vector<uint32_t> &end,
        std::unique_ptr<StringReader> &&chr_levels
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