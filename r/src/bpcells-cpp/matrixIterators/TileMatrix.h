// Copyright 2022 BPCells contributors
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
#include <libdivide/libdivide.h>

namespace BPCells {

// Output cell x tile matrix (rows = cell_id, col = tile_id)
// Regions are given as half-open format and cannot overlap
// Each region can select a different tile width
// The last tile in a region may be truncated if the region is not an even multiple of the tile
// width
// - Valid MODE options:
//    - 0: Count insertions per tile. One fragment can count twice if it lands
//      in same tile twice
//    - 1: Count fragments per tile. Same as MODE == 0 but fragment counts at most
//      once per tile
class TileMatrix : public MatrixLoader<uint32_t> {
  private:
    class Tile {
      public:
        uint32_t chr, start, end, output_idx;
        struct libdivide::libdivide_u32_t width;
    };

    std::unique_ptr<FragmentLoader> frags;
    std::unique_ptr<StringReader> chr_levels;
    MatrixAccumulator<uint32_t> accumulator;
    std::vector<Tile> sorted_tiles;
    std::vector<Tile> active_tiles;
    uint32_t next_completed_tile = 0; // All columns below this number are ready for output
    uint32_t current_output_tile =
        UINT32_MAX; // Index of the current column we are trying to output
    uint32_t next_active_tile =
        0; // Index of sorted_tiles for the next Tile where we'll find overlaps
    uint32_t n_tiles = 0;

    std::string tile_name; // buffer to use to store the tile name
    const bool count_fragments;

    void loadFragments();

  public:
    // Note: It's the caller's responsibility to make sure that
    // the FragmentLoader will not be deleted while this TileMatrix is still alive
    TileMatrix(
        std::unique_ptr<FragmentLoader> &&frags,
        const std::vector<uint32_t> &chr,
        const std::vector<uint32_t> &start,
        const std::vector<uint32_t> &end,
        const std::vector<uint32_t> &width,
        std::unique_ptr<StringReader> &&chr_levels,
        bool count_fragments = false
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