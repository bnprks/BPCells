// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "RegionSelect.h"

namespace BPCells {

RegionSelect::RegionSelect(
    std::unique_ptr<FragmentLoader> &&loader,
    const std::vector<uint32_t> &chr,
    const std::vector<uint32_t> &start,
    const std::vector<uint32_t> &end,
    std::unique_ptr<StringReader> &&chr_levels,
    bool invert_selection
)
    : FragmentLoaderWrapper(std::move(loader))
    , invert_selection(invert_selection)
    , chr_levels(std::move(chr_levels)) {

    if (chr.size() != start.size() || chr.size() != end.size())
        throw std::invalid_argument("chr, start, and end must all be same length");

    for (size_t i = 0; i < chr.size(); i++) {
        if (chr[i] >= this->chr_levels->size())
            throw std::invalid_argument("chr has values higher than length of chr_levels");
        Region r;
        r.start = start[i];
        r.end = end[i];
        r.chr = chr[i];
        sorted_regions.push_back(r);
    }

    // Sort regions by start coord since that's the order we'll see them
    std::sort(sorted_regions.begin(), sorted_regions.end(), [](Region a, Region b) {
        if (a.chr != b.chr) return a.chr < b.chr;
        else if (a.start != b.start) return a.start < b.start;
        else return a.end < b.end;
    });

    // Combine overlapping regions
    uint32_t out_idx = 0;
    for (uint32_t i = 1; i < sorted_regions.size(); i++) {
        Region &a = sorted_regions[out_idx];
        Region &b = sorted_regions[i];
        if (a.chr == b.chr && a.end >= b.start) {
            a.end = std::max(a.end, b.end);
        } else {
            out_idx += 1;
            sorted_regions[out_idx] = b;
        }
    }
    sorted_regions.resize(out_idx + 1);

    // Sentinel value at end of sorted_regions. Use chr UINT32_MAX-1 to avoid conflicts with
    // findChrIDTranslation
    sorted_regions.push_back({UINT32_MAX - 1, UINT32_MAX, UINT32_MAX});
}

void RegionSelect::seek(uint32_t chr_id, uint32_t base) {
    loader->seek(chr_id, base);
    if ((int64_t) chr_id < loader->chrCount()) {
        current_chr_id = findChrIDTranslation(loader->chrNames(loader->currentChr()));
        active_region = computeNextActiveRegion(current_chr_id, base);
    }
    did_seek_active_region = false;
}

void RegionSelect::restart() {
    active_region = 0;
    did_seek_active_region = false;
    loader->restart();
}

bool RegionSelect::nextChr() {
    bool ret = loader->nextChr();
    current_chr_id = findChrIDTranslation(loader->chrNames(loader->currentChr()));
    if (ret) {
        active_region = computeNextActiveRegion(current_chr_id, 0);
        did_seek_active_region = false;
    } else {
        active_region = sorted_regions.size() - 1;
    }
    return ret;
}

bool RegionSelect::load() {
    // Overlap procedure:
    //  1. Scan until frag.start >= region.start, marking any overlaps that happen from frag.end
    //  2. Binary search for frag.start >= region.end to find all the remaining overlaps
    //  3. Increment active_region and continue scan from step 1
    loaded = 0;
    while (loaded == 0) {
        if (!loader->load()) return false;
        uint32_t capacity = loader->capacity();
        uint32_t *start = loader->startData();
        uint32_t *end = loader->endData();
        uint32_t *cell = loader->cellData();
        uint32_t i = 0;
        while (i < capacity) {
            Region r = sorted_regions[active_region];
            // Check if we've gotten to a region beyond the current chromosome
            if (r.chr != current_chr_id) {
                if (invert_selection) {
                    std::memmove(&cell[loaded], &cell[i], sizeof(uint32_t) * (capacity - i));
                    std::memmove(&end[loaded], &end[i], sizeof(uint32_t) * (capacity - i));
                    std::memmove(&start[loaded], &start[i], sizeof(uint32_t) * (capacity - i));
                    loaded += capacity - i;
                    return true;
                } else {
                    return loaded > 0;
                }
            }
            // 1. Scan until frag.start >= region.start, marking any overlaps that happen from
            // frag.end
            while (i < capacity && start[i] < r.start) {
                cell[loaded] = cell[i];
                start[loaded] = start[i];
                end[loaded] = end[i];
                loaded += invert_selection !=
                          (end[i] > r.start); // Don't compare end[i] < r.end, since spanning the
                                              // region counts as overlapping
                i++;
            }
            if (i >= capacity) break;
            // 2. Binary search for frag.start >= region.end to find all the remaining overlaps
            auto easy_overlaps = std::lower_bound(&start[i], &start[capacity], r.end) - &start[i];
            if (!invert_selection) {
                std::memmove(&cell[loaded], &cell[i], sizeof(uint32_t) * easy_overlaps);
                std::memmove(&end[loaded], &end[i], sizeof(uint32_t) * easy_overlaps);
                std::memmove(&start[loaded], &start[i], sizeof(uint32_t) * easy_overlaps);
                loaded += easy_overlaps;
            }
            i += easy_overlaps;
            //  3. Increment active_region and continue scan from step 1
            if (i < capacity) {
                active_region++;
                did_seek_active_region = false;
            }
        }
        // If loaded == 0, try seeking to get closer to viable fragments, assuming we haven't
        // already done a seek on the current region (to account for slop in the seeking process)
        if (loaded == 0 && !did_seek_active_region) {
            if (invert_selection)
                loader->seek(loader->currentChr(), sorted_regions[active_region].end);
            else loader->seek(loader->currentChr(), sorted_regions[active_region].start);
            did_seek_active_region = true;
        }
    }
    return loaded > 0;
}

uint32_t RegionSelect::capacity() const { return loaded; }

uint32_t RegionSelect::computeNextActiveRegion(uint32_t chr, uint32_t base) const {
    // Find first region where end > start
    auto it = std::upper_bound(
        sorted_regions.begin(),
        sorted_regions.end(),
        std::pair{chr, base},
        [](std::pair<uint32_t, uint32_t> value, Region r) {
            if (value.first != r.chr) return value.first < r.chr;
            return value.second < r.end;
        }
    );
    uint32_t pos = it - sorted_regions.begin();
    // If we're on the wrong chromosme, our search will return sorted_regions.size(), but
    // we want to be looking at the sentinel region
    return pos - (pos == sorted_regions.size());
}

// Return the index in chr_levels of the given chr_name. Return UINT32_MAX if chromosome not found
uint32_t RegionSelect::findChrIDTranslation(const char *chr_name) const {
    if (chr_name == NULL)
        throw std::runtime_error("RegionSelect saw NULL chrName from fragment loader");
    for (uint32_t i = 0; i < chr_levels->size(); i++) {
        const char *chr_level = chr_levels->get(i);
        if (chr_level == NULL)
            throw std::runtime_error("RegionSelect saw NULL chrName from fragment loader");
        if (strcmp(chr_name, chr_level) == 0) {
            return i;
        }
    }
    return UINT32_MAX;
}

} // end namespace BPCells
