// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "simd/overlaps.cpp"
#include <hwy/foreach_target.h>

#include <algorithm>
#include <cassert>

#include "libdivide-inl.h"
#include <hwy/highway.h>

HWY_BEFORE_NAMESPACE();

namespace BPCells::simd::HWY_NAMESPACE {

template <int MODE>
uint32_t tile_overlaps(
    const uint32_t *cell_ids,
    const uint32_t *starts,
    const uint32_t *ends,
    uint32_t n,
    const uint32_t tile_start,
    const uint32_t tile_end,
    const uint32_t tile_output_idx,
    const struct libdivide::libdivide_u32_t *tile_width,
    uint32_t *HWY_RESTRICT cell_id_out,
    uint32_t *HWY_RESTRICT tile_idx_out
) {
    using namespace hwy::HWY_NAMESPACE;
    ScalableTag<uint32_t> d;

    static_assert(MODE == 0 || MODE == 1);

    uint32_t written = 0;

    size_t i = 0;
    const size_t N = Lanes(d);
    auto t_start = Set(d, tile_start);
    auto t_start_p1 = Set(d, tile_start + 1);
    auto t_end = Set(d, tile_end);
    auto t_output_idx = Set(d, tile_output_idx);

    // Vectorized loop
    for (; i + N <= n && starts[i] < tile_end; i += N) {
        auto cell_id = LoadU(d, cell_ids + i);
        auto start = LoadU(d, starts + i);
        auto end = LoadU(d, ends + i);

        auto start_overlap = And(Ge(start, t_start), Lt(start, t_end));
        auto start_idx = Add(t_output_idx, libdivide_u32(Sub(start, t_start), tile_width));

        CompressStore(cell_id, start_overlap, d, cell_id_out + written);
        written += CompressStore(start_idx, start_overlap, d, tile_idx_out + written);

        auto end_overlap = And(Gt(end, t_start), Le(end, t_end));

        if constexpr (MODE == 1) {
            end_overlap = AndNot(start_overlap, end_overlap);
        }

        auto end_idx = Add(t_output_idx, libdivide_u32(Sub(end, t_start_p1), tile_width));

        CompressStore(cell_id, end_overlap, d, cell_id_out + written);
        written += CompressStore(end_idx, end_overlap, d, tile_idx_out + written);
    }

    // Cleanup loop
    for (; i < n && starts[i] < tile_end; i++) {
        bool start_overlap = starts[i] >= tile_start && starts[i] < tile_end;
        uint32_t start_idx =
            tile_output_idx + libdivide::libdivide_u32_do(starts[i] - tile_start, tile_width);

        if (start_overlap) {
            cell_id_out[written] = cell_ids[i];
            tile_idx_out[written] = start_idx;
            written += 1;
        }

        bool end_overlap = ends[i] > tile_start && ends[i] <= tile_end;
        if constexpr (MODE == 1) {
            end_overlap = end_overlap && !start_overlap;
        }
        uint32_t end_idx =
            tile_output_idx + libdivide::libdivide_u32_do(ends[i] - 1 - tile_start, tile_width);

        if (end_overlap) {
            cell_id_out[written] = cell_ids[i];
            tile_idx_out[written] = end_idx;
            written += 1;
        }
    }

    return written;
}

uint32_t tile_overlaps_insertion(
    const uint32_t *cell_ids,
    const uint32_t *starts,
    const uint32_t *ends,
    uint32_t n,
    const uint32_t tile_start,
    const uint32_t tile_end,
    const uint32_t tile_output_idx,
    const struct libdivide::libdivide_u32_t *tile_width,
    uint32_t *HWY_RESTRICT cell_id_out,
    uint32_t *HWY_RESTRICT tile_idx_out
) {
    return tile_overlaps<0>(
        cell_ids,
        starts,
        ends,
        n,
        tile_start,
        tile_end,
        tile_output_idx,
        tile_width,
        cell_id_out,
        tile_idx_out
    );
}

uint32_t tile_overlaps_fragment(
    const uint32_t *cell_ids,
    const uint32_t *starts,
    const uint32_t *ends,
    uint32_t n,
    const uint32_t tile_start,
    const uint32_t tile_end,
    const uint32_t tile_output_idx,
    const struct libdivide::libdivide_u32_t *tile_width,
    uint32_t *HWY_RESTRICT cell_id_out,
    uint32_t *HWY_RESTRICT tile_idx_out
) {
    return tile_overlaps<1>(
        cell_ids,
        starts,
        ends,
        n,
        tile_start,
        tile_end,
        tile_output_idx,
        tile_width,
        cell_id_out,
        tile_idx_out
    );
}

template <int MODE>
uint32_t peak_overlaps(
    const uint32_t *cell_ids,
    const uint32_t *starts,
    const uint32_t *ends,
    uint32_t n,
    const uint32_t peak_start,
    const uint32_t peak_end,
    uint32_t *cell_id_out,
    uint32_t *count_out
) {
    using namespace hwy::HWY_NAMESPACE;
    ScalableTag<uint32_t> d;

    static_assert(MODE == 0 || MODE == 1 || MODE == 2);

    uint32_t written = 0;

    size_t i = 0;
    const size_t N = Lanes(d);
    auto p_start = Set(d, peak_start);
    auto p_end = Set(d, peak_end);
    auto one = Set(d, 1);
    auto two = Set(d, 2);
    // Vectorized loop
    for (; i + N <= n && starts[i] < peak_end; i += N) {
        auto cell_id = LoadU(d, cell_ids + i);
        auto start = LoadU(d, starts + i);
        auto end = LoadU(d, ends + i);

        if constexpr (MODE == 0 || MODE == 1) {
            auto start_overlap = And(Ge(start, p_start), Lt(start, p_end));
            auto end_overlap = And(Gt(end, p_start), Le(end, p_end));
            auto has_overlap = Or(start_overlap, end_overlap);
            if constexpr (MODE == 1) {
                written += CompressStore(cell_id, has_overlap, d, cell_id_out + written);
            }
            if constexpr (MODE == 0) {
                auto overlap_count = IfThenElse(And(start_overlap, end_overlap), two, one);
                CompressStore(cell_id, has_overlap, d, cell_id_out + written);
                written += CompressStore(overlap_count, has_overlap, d, count_out + written);
            }
        }
        if constexpr (MODE == 2) {
            auto has_overlap = And(Lt(start, p_end), Gt(end, p_start));
            written += CompressStore(cell_id, has_overlap, d, cell_id_out + written);
        }
    }

    // Cleanup loop
    for (; i < n && starts[i] < peak_end; i++) {
        if constexpr (MODE == 0 || MODE == 1) {
            bool start_overlap = starts[i] >= peak_start && starts[i] < peak_end;
            bool end_overlap = ends[i] > peak_start && ends[i] <= peak_end;
            if (start_overlap || end_overlap) {
                if constexpr (MODE == 1) {
                    cell_id_out[written] = cell_ids[i];
                    written += 1;
                }
                if constexpr (MODE == 0) {
                    cell_id_out[written] = cell_ids[i];
                    count_out[written] = start_overlap + end_overlap;
                    written += 1;
                }
            }
        }
        if constexpr (MODE == 2) {
            bool overlap = starts[i] < peak_end && ends[i] > peak_start;
            if (overlap) {
                cell_id_out[written] = cell_ids[i];
                written += 1;
            }
        }
    }

    if constexpr (MODE == 1 || MODE == 2) {
        std::fill_n(count_out, written, 1);
    }

    return written;
}

uint32_t peak_overlaps_insertion(
    const uint32_t *cell_ids,
    const uint32_t *starts,
    const uint32_t *ends,
    uint32_t n,
    const uint32_t peak_start,
    const uint32_t peak_end,
    uint32_t *cell_id_out,
    uint32_t *count_out
) {
    return peak_overlaps<0>(
        cell_ids, starts, ends, n, peak_start, peak_end, cell_id_out, count_out
    );
}

uint32_t peak_overlaps_fragment(
    const uint32_t *cell_ids,
    const uint32_t *starts,
    const uint32_t *ends,
    uint32_t n,
    const uint32_t peak_start,
    const uint32_t peak_end,
    uint32_t *cell_id_out,
    uint32_t *count_out
) {
    return peak_overlaps<1>(
        cell_ids, starts, ends, n, peak_start, peak_end, cell_id_out, count_out
    );
}

uint32_t peak_overlaps_overlap(
    const uint32_t *cell_ids,
    const uint32_t *starts,
    const uint32_t *ends,
    uint32_t n,
    const uint32_t peak_start,
    const uint32_t peak_end,
    uint32_t *cell_id_out,
    uint32_t *count_out
) {
    return peak_overlaps<2>(
        cell_ids, starts, ends, n, peak_start, peak_end, cell_id_out, count_out
    );
}

} // namespace BPCells::simd::HWY_NAMESPACE
HWY_AFTER_NAMESPACE();

#if HWY_ONCE

namespace BPCells::simd {

HWY_EXPORT(tile_overlaps_insertion);
HWY_EXPORT(tile_overlaps_fragment);

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
) {
    assert(mode == 0 || mode == 1);

    if (mode == 0) {
        return HWY_DYNAMIC_DISPATCH(tile_overlaps_insertion)(
            cell_ids,
            starts,
            ends,
            n,
            tile_start,
            tile_end,
            tile_output_idx,
            tile_width,
            cell_id_out,
            tile_idx_out
        );
    } else { // (mode == 1)
        return HWY_DYNAMIC_DISPATCH(tile_overlaps_fragment)(
            cell_ids,
            starts,
            ends,
            n,
            tile_start,
            tile_end,
            tile_output_idx,
            tile_width,
            cell_id_out,
            tile_idx_out
        );
    }
}

HWY_EXPORT(peak_overlaps_insertion);
HWY_EXPORT(peak_overlaps_fragment);
HWY_EXPORT(peak_overlaps_overlap);

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
) {
    assert(mode == 0 || mode == 1 || mode == 2);

    if (mode == 0) {
        return HWY_DYNAMIC_DISPATCH(peak_overlaps_insertion)(
            cell_ids, starts, ends, n, peak_start, peak_end, cell_id_out, count_out
        );
    } else if (mode == 1) {
        return HWY_DYNAMIC_DISPATCH(peak_overlaps_fragment)(
            cell_ids, starts, ends, n, peak_start, peak_end, cell_id_out, count_out
        );
    } else {
        return HWY_DYNAMIC_DISPATCH(peak_overlaps_overlap)(
            cell_ids, starts, ends, n, peak_start, peak_end, cell_id_out, count_out
        );
    }
}

} // namespace BPCells::simd
#endif