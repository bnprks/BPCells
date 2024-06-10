// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "PeakMatrix.h"
#include "../simd/math.h"
#include "../simd/overlaps.h"

namespace BPCells {


PeakMatrix::PeakMatrix(
    std::unique_ptr<FragmentLoader> &&frags,
    const std::vector<uint32_t> &chr,
    const std::vector<uint32_t> &start,
    const std::vector<uint32_t> &end,
    std::unique_ptr<StringReader> &&chr_levels,
    uint32_t mode
)
    : mode(mode)
    , frags(std::move(frags))
    , chr_levels(std::move(chr_levels))
    , end_sorted_lookup(start.size())
    , n_peaks(start.size()) {
    if (this->frags->cellCount() < 0)
        throw std::invalid_argument(
            "frags must have a known cell count. Consider using a cell selection to define the "
            "number of cells."
        );

    if (chr.size() != start.size() || chr.size() != end.size())
        throw std::invalid_argument("chr, start, and end must all be same length");

    // Check that chr name matches for all the available chrNames in frags
    for (uint32_t i = 0; i < this->chr_levels->size(); i++) {
        const char *chr_name_frag = this->frags->chrNames(i);
        const char *chr_name_args = this->chr_levels->get(i);
        if (chr_name_frag != NULL &&
            (chr_name_args == NULL || strcmp(chr_name_frag, chr_name_args) != 0)) {
            throw std::runtime_error(
                std::string("PeakMatrix encountered fragment with incorrect chrLevel: ") +
                std::string(chr_name_frag) + std::string(" expected: ") + std::string(chr_name_args)
            );
        }
    }

    Peak prev;
    for (size_t i = 0; i < chr.size(); i++) {
        if (chr[i] >= this->chr_levels->size())
            throw std::invalid_argument("chr has values higher than length of chr_levels");
        Peak p;
        p.start = start[i];
        p.end = end[i];
        p.chr = chr[i];
        sorted_peaks.push_back(p);
        if (i > 0) {
            bool ordered = true;
            if (prev.chr != p.chr) ordered = prev.chr < p.chr;
            else if (prev.end != p.end) ordered = prev.end < p.end;
            else ordered = prev.start <= p.start;
            if (!ordered) {
                throw std::invalid_argument("Peaks are not in sorted order by (chr, end, start)");
            }
        }
        prev = p;
    }

    // Make a lookup so that end_sorted_lookup[i] gives the end-sorted index of start-sorted index i
    // This lets us translate from the order we see new peaks into the order we'll output them
    for (size_t i = 0; i < chr.size(); i++) {
        end_sorted_lookup[i] = i;
    }

    const std::vector<Peak> &peaks_ref = sorted_peaks;
    std::sort(
        end_sorted_lookup.begin(),
        end_sorted_lookup.end(),
        [&peaks_ref](uint32_t a, uint32_t b) {
            Peak pa = peaks_ref[a];
            Peak pb = peaks_ref[b];
            if (pa.chr != pb.chr) return pa.chr < pb.chr;
            else if (pa.start != pb.start) return pa.start < pb.start;
            else return pa.end < pb.end;
        }
    );

    // Sort peaks by start coord since that's the order we'll see them
    std::sort(sorted_peaks.begin(), sorted_peaks.end(), [](Peak a, Peak b) {
        if (a.chr != b.chr) return a.chr < b.chr;
        else if (a.start != b.start) return a.start < b.start;
        else return a.end < b.end;
    });

    // Sentinel value at end of sorted_peaks
    sorted_peaks.push_back({UINT32_MAX, UINT32_MAX, UINT32_MAX});

    // Call nextChr to start fragments (analagous to what's done in loadFragments)
    if (!this->frags->nextChr()) {
        next_completed_peak = sorted_peaks.size() - 1;
        active_peaks.clear();
        return;
    }
    // Check that chr name matches
    const char *chr_name_frag = this->frags->chrNames(this->frags->currentChr());
    const char *chr_name_args = this->chr_levels->get(this->frags->currentChr());
    if (chr_name_frag == NULL || chr_name_args == NULL ||
        strcmp(chr_name_frag, chr_name_args) != 0) {
        throw std::runtime_error(
            std::string("PeakMatrix encountered fragment with incorrect chrLevel: ") +
            std::string(chr_name_frag) + std::string(" expected: ") + std::string(chr_name_args)
        );
    }
    while (sorted_peaks[next_completed_peak].chr < this->frags->currentChr()) {
        next_completed_peak++;
    }
    next_active_peak = next_completed_peak;
}

uint32_t PeakMatrix::rows() const { return frags->cellCount(); }
uint32_t PeakMatrix::cols() const { return n_peaks; }

const char *PeakMatrix::rowNames(uint32_t row) {
    return frags->cellNames(row);
}
const char *PeakMatrix::colNames(uint32_t col) {
    if (col >= n_peaks) return NULL;
    auto peak = sorted_peaks[end_sorted_lookup[col]];
    peak_name.clear();
    peak_name += frags->chrNames(peak.chr);
    peak_name += ":";
    peak_name += std::to_string(peak.start);
    peak_name += "-";
    peak_name += std::to_string(peak.end);
    return peak_name.c_str();
}

void PeakMatrix::restart() {
    accumulator.clear();
    active_peaks.clear();
    next_completed_peak = 0;
    next_active_peak = 0;
    current_output_peak = UINT32_MAX;
}
void PeakMatrix::seekCol(uint32_t col) {
    if (!frags->isSeekable())
        throw std::runtime_error("Can't seek a PeakMatrix if the fragments aren't seekable");
    next_active_peak = end_sorted_lookup[col];
    next_completed_peak = 0;
    current_output_peak = col - 1;
    active_peaks.clear();
    accumulator.clear();
    // Don't need to do any frags->seek call here since loadFragments will do it
    nextCol();
}

bool PeakMatrix::nextCol() {
    current_output_peak += 1;

    if (current_output_peak >= n_peaks) {
        current_output_peak -= 1;
        return false;
    }
    if (current_output_peak >= next_completed_peak) loadFragments();
    accumulator.discard_until(current_output_peak);
    return true;
}

uint32_t PeakMatrix::currentCol() const {
    return current_output_peak;
}

bool PeakMatrix::load() {
    return accumulator.load(current_output_peak, 1024);
}

uint32_t PeakMatrix::capacity() const {
    return accumulator.capacity();
}

uint32_t *PeakMatrix::rowData() { return accumulator.rowData(); }

uint32_t *PeakMatrix::valData() { return accumulator.valData(); }

void PeakMatrix::loadFragments() {
    // Load fragments data until we hit enough loaded for output
    // Post-conditions:
    //  - accumulator is ready for loading, and has complete data on peaks at least until
    //  current_output_peak

    // - If no peaks active, seek to the next peak
    // - Infinite loop
    //   1. load fragments
    //     - if false, load next chromosome, confirm its name, and
    //       adjust the next peak to load
    //   2. activate new peaks if relevant
    //   3. iterate through the available fragments, tallying overlaps
    //      - Iterate peak outside & fragments inside
    //   4. break if we're ready to accumulate and have next_completed_peak > current_output_peak
    if (accumulator.ready_for_loading() && next_completed_peak > current_output_peak) return;

    if (next_active_peak == sorted_peaks.size()) return;

    if (active_peaks.size() == 0 && frags->isSeekable()) {
        frags->seek(sorted_peaks[next_active_peak].chr, sorted_peaks[next_active_peak].start);
    }

    while (true) {
        // Load fragments, and check for end of chromosomes
        while (!frags->load()) {
            uint32_t prev_chr_id = frags->currentChr();
            if (!frags->nextChr()) {
                next_completed_peak = sorted_peaks.size() - 1;
                active_peaks.clear();
                return;
            }
            if (frags->currentChr() <= prev_chr_id) {
                throw std::runtime_error(
                    "PeakMatrix encountered fragments with out of order chromosome IDs. Please "
                    "save + load fragments before passing to PeakMatrix to fix this issue."
                );
            }
            // Check that chr name matches
            const char *chr_name_frag = frags->chrNames(frags->currentChr());
            const char *chr_name_args = chr_levels->get(frags->currentChr());
            if (chr_name_frag == NULL || chr_name_args == NULL ||
                strcmp(chr_name_frag, chr_name_args) != 0) {
                throw std::runtime_error(
                    std::string("PeakMatrix encountered fragment with incorrect chrLevel: ") +
                    std::string(chr_name_frag) + std::string(" expected: ") +
                    std::string(chr_name_args)
                );
            }
            while (sorted_peaks[next_completed_peak].chr < frags->currentChr()) {
                next_completed_peak++;
            }
            next_active_peak = next_completed_peak;
            active_peaks.clear();
        }
        uint32_t capacity = frags->capacity();
        const uint32_t *start_data = frags->startData();
        const uint32_t *end_data = frags->endData();
        const uint32_t *cell_data = frags->cellData();

        uint32_t i = 0;
        uint32_t end_max = 0;

        const uint32_t max_items = 256;
        uint32_t overlap_cell[max_items];
        uint32_t overlap_count[max_items];
        // Loop through reads in blocks of 128 at a time
        while (i < capacity) {
            uint32_t items = std::min(max_items, capacity - i);
            end_max = std::max(end_max, simd::max(end_data + i, items));

            // Check for new peaks to activate
            while (sorted_peaks[next_active_peak].chr == frags->currentChr() &&
                   sorted_peaks[next_active_peak].start < end_max) {

                active_peaks.push_back(sorted_peaks[next_active_peak]);
                active_peaks.back().chr =
                    end_sorted_lookup[next_active_peak]; // Sloppy, but overload active_peaks chr
                                                         // field to be output index
                next_active_peak += 1;
            }

            // For each active peak, iterate through the fragments & tally overlaps
            for (uint32_t j = 0; j < active_peaks.size(); j++) {
                Peak p = active_peaks[j];
                
                uint32_t n_overlaps = simd::peak_overlaps(
                    cell_data + i,
                    start_data + i,
                    end_data + i,
                    items,
                    p.start,
                    p.end,
                    overlap_cell,
                    overlap_count,
                    mode
                );

                for (uint32_t k = 0; k < n_overlaps; k++) {
                    accumulator.add_one(
                        p.chr,
                        overlap_cell[k],
                        overlap_count[k]
                    );
                }

                // Remove the peak from active_peaks if we're done, and mark the
                // next completed peak on our list
                if (start_data[i + items - 1] >= p.end) {
                    while (sorted_peaks[next_completed_peak].chr == frags->currentChr() &&
                           sorted_peaks[next_completed_peak].end <= p.end &&
                           sorted_peaks[next_completed_peak].start <= p.start) {
                        next_completed_peak += 1;
                    }
                    std::swap(active_peaks.back(), active_peaks[j]);
                    active_peaks.pop_back();
                    j -= 1;
                }
            }
            i += items;
        }

        // break if we're ready to accumulate and have next_completed_peak >current_output_peak
        if (accumulator.ready_for_loading() && next_completed_peak > current_output_peak) {
            break;
        }
    }
}

} // end namespace BPCells
