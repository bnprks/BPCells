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
#include "FragmentIterator.h"

namespace BPCells {

// Filter a FragmentsLoader to include only fragments that
// overlap (or don't overlap) a given set of genomic regions.
// Note that if the read spans a given region it will also count as overlapping even
// if neither the start or end are in the region
class RegionSelect : public FragmentLoaderWrapper {
  private:
    class Region {
      public:
        uint32_t chr, start, end;
    };

    std::vector<Region> sorted_regions; // All regions sorted by (chr, start)
    uint32_t active_region = 0;         // Index in sorted_regions of current active region
    uint32_t loaded, current_chr_id;
    bool did_seek_active_region = false; // Set to true once we have seeked to the active region
    bool invert_selection; // If true, exclude fragments overlapping regions rather than including

    std::unique_ptr<StringReader> chr_levels;

    // Return the index of the first region in sorted_regions where r.chr >= chr and r.end > base
    uint32_t computeNextActiveRegion(uint32_t chr, uint32_t base) const;

    // Return the index in chr_levels of the given chr_name. Return UINT32_MAX if chromosome not
    // found
    uint32_t findChrIDTranslation(const char *chr_name) const;

  public:
    RegionSelect(
        std::unique_ptr<FragmentLoader> &&loader,
        const std::vector<uint32_t> &chr,
        const std::vector<uint32_t> &start,
        const std::vector<uint32_t> &end,
        std::unique_ptr<StringReader> &&chr_levels,
        bool invert_selection
    );

    ~RegionSelect() = default;

    void seek(uint32_t chr_id, uint32_t base) override;
    void restart() override;

    bool nextChr() override;

    bool load() override;
    uint32_t capacity() const override;
};

} // end namespace BPCells