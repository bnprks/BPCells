// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include "../../vendor/libdivide/libdivide.h"
#include <cstdint>

namespace BPCells::simd {

/**
 * @brief Calculate tile overlaps
 *
 * Pre-conditions:
 *
 * - `cell_id_out` and `tile_idx_out` must have available space for at least `2*n` overlaps
 * - `cell_id_out` and `tile_idx_out` must not overlap in memory with any other parameters
 *
 * Input Data:
 *
 * - Fragments:`n` (cell_id, start, end) tuples passed as `cell_ids`, `starts`, `ends`. 
 *   Must be sorted by start coordinate
 * - Tile: Given with coordinates (tile_start, tile_end), pre-constructed tile_width, and base
 * output_idx
 *
 * Output:
 *
 * - For every overlap, write the cell_id to `cell_id_out`, and the tile index to `tile_idx_out`.
 *   Tile index of a coord is calculated as `tile_output_idx + (coord - tile_start) / tile_width`
 * - Overlap results are written contiguously in arbitrary order
 * - Return `m`, the total number of overlaps calculated
 *
 * Mode details:
 *
 *  - 0 == Insertion: Count insertions overlaps. One fragment can count twice if it lands
 *    in same region twice.
 *  - 1 == Fragment: Count fragments. Same as Insertion but fragment cannot count
 *    twice per region.
 *  - (Note: Mode 2 == Overlap is not supported for tile_overlaps, due to the possibility of an
 * unbounded number of overlaps created per-fragment for long fragments overlapping many tiles)
 */
uint32_t tile_overlaps(
    const uint32_t *cell_ids,
    const uint32_t *starts,
    const uint32_t *ends,
    uint32_t n,
    const uint32_t tile_start,
    const uint32_t tile_end,
    const uint32_t tile_output_idx,
    const struct libdivide::libdivide_u32_t *tile_width,
    uint32_t *cell_id_out,
    uint32_t *tile_idx_out,
    const uint32_t mode
);

/**
 * @brief Calculate peak overlaps
 *
 * Pre-conditions:
 *
 * - `cell_id_out` and `count_out` must have available space for at least `n` overlaps
 * - `cell_id_out` and `count_out` must not overlap in memory with any other parameters
 *
 * Input Data:
 *
 * - Fragments:`n` (cell_id, start, end) tuples passed as `cell_ids`, `starts`, `ends`
 *   Must be sorted by start coordinate
 * - Peak: Given with coordinates (peak_start, peak_end)
 *
 * Output:
 *
 * - For every overlap, write the cell_id to `cell_id_out`, and the count of the overlaps for the
 * fragment to `count_out`.
 * - Overlap results are written contiguously in arbitrary order
 * - Return `m`, the total number of overlaps calculated
 *
 * Mode details:
 *
 *  - 0 == Insertion: Count insertions overlaps. One fragment can count twice if it lands
 *    in same region twice.
 *  - 1 == Fragment: Count fragments. Same as Insertion but fragment cannot count
 *    twice per region.
 *  - 2 == Overlap: Count overlaps. Same as Fragment but fragments that fully
 *    surround a region will also be counted for that region.
 */
uint32_t peak_overlaps(
    const uint32_t *cell_ids,
    const uint32_t *starts,
    const uint32_t *ends,
    uint32_t n,
    const uint32_t peak_start,
    const uint32_t peak_end,
    uint32_t *cell_id_out,
    uint32_t *count_out,
    const uint32_t mode
);

} // namespace BPCells::simd