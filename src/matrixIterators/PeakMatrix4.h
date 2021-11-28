#pragma once
#include <algorithm>
#include <cassert>

#include "../fragmentIterators/FragmentsIterator.h"
#include "../fragmentIterators/OverlapIterator.h"
#include "MatrixIterator.h"

namespace BPCells {

// Output cell x peak matrix (rows = cell_id, col = peak_id)
// Regions are given as half-open format
class PeakMatrix4 : public MatrixLoader<uint32_t> {
public:
    PeakMatrix4(FragmentsLoader &frags, 
            std::vector<uint32_t> &chr, std::vector<uint32_t> &start, std::vector<uint32_t> &end, 
            std::vector<std::string>  &chr_levels);
    
    // Reset the loader to start from the beginning
    void restart() override;

    // Return the count of rows and columns
    uint32_t rows() const override { return frags.cellCount(); }
    uint32_t cols() const override { return output_peak_index.size(); }

    bool nextCol() override;
    uint32_t currentCol() const override;

    int32_t load(uint32_t count, SparseVector<uint32_t> buffer) override;
private:
    class Count {
    public:
        uint32_t cell;
        uint32_t count;
    };
    std::vector<std::vector<Count>> inactive_buffers; // Object pool of vectors with non-zero capacity
    std::vector<std::vector<Count>> peak_buffers; // Count buffers for each peak, in order by sorted peaks
    
    FragmentsLoader &frags;
    OverlapIterator overlaps;

    std::vector<uint32_t> output_peak_index; // Output column index for each peak, in same order as peak_buffers

    uint32_t output_peak;
    std::vector<uint32_t> output_counts; // Length cells, with counts to output
    std::vector<uint32_t> output_cells; // List of non-zero indices in output_counts
    
    // Load the next peak and populate the output_counts + output_cells buffers
    bool loadPeak();

    void beginPeakBuffer(std::vector<Count> &buf) {
        if (inactive_buffers.empty()) return;
        std::swap(buf, inactive_buffers.back());
        inactive_buffers.pop_back();
    }
    void tallyPeakBuffer(std::vector<Count> &buf, uint32_t peak_idx) {
        for(Count c : buf) {
            if (output_counts[c.cell] == 0)
                output_cells.push_back(c.cell);
            output_counts[c.cell] += c.count;
        }

        output_peak = output_peak_index[peak_idx];
    };
    void clearPeakBuffer(std::vector<Count> &buf) {
        buf.resize(0);
        if (buf.capacity() > 0) {
            inactive_buffers.push_back(std::move(buf));
            buf = std::vector<Count>();
        }
    };

    // Get sorted regions, but adjust indexing to be the sorted order
    // This is a bit of a workaround so I can initialize overlaps in the
    // constructor's initializer list
    static OverlapIterator::RegionList reindexSortRegions(uint32_t n_chrs, std::vector<uint32_t> &chr, std::vector<uint32_t> &start, std::vector<uint32_t> &end) {
        OverlapIterator::RegionList sorted_regions = OverlapIterator::sortRegions(n_chrs, chr, start, end);
        size_t peak_index = 0;
        for (size_t i = 0; i < sorted_regions.size(); i++) {
            for (size_t j = 0; j < sorted_regions[i].size(); j++) {
                //output_peak_index[peak_index] = sorted_regions[i][j].index;
                sorted_regions[i][j].index = peak_index;
                peak_index++;
            }
        }
        return sorted_regions;
    }
};


PeakMatrix4::PeakMatrix4(FragmentsLoader &frags, 
        std::vector<uint32_t> &chr, std::vector<uint32_t> &start, std::vector<uint32_t> &end, 
        std::vector<std::string>  &chr_levels) :
        peak_buffers(chr.size()), 
        frags(frags),
        overlaps(frags, reindexSortRegions(chr_levels.size(), chr, start, end), chr_levels),
        output_peak_index(chr.size()) {

    if (frags.cellCount() < 0)
        throw std::invalid_argument("frags must have a known cell count. Consider using a cell selection to define the number of cells.");
    
    output_counts.resize(frags.cellCount(), 0);
    OverlapIterator::RegionList sorted_regions = 
        OverlapIterator::sortRegions(chr_levels.size(), chr, start, end);
    
    size_t peak_index = 0;
    for (size_t i = 0; i < sorted_regions.size(); i++) {
        for (size_t j = 0; j < sorted_regions[i].size(); j++) {
            output_peak_index[peak_index] = sorted_regions[i][j].index;
            //sorted_regions[i][j].index = peak_index;
            peak_index++;
        }
    }
}

void PeakMatrix4::restart() {
    // Clear out the active peak buffers
    for (auto b : peak_buffers) {
        clearPeakBuffer(b);
    }

    output_cells.resize(0);
    for (size_t i = 0; i < output_counts.size(); i++)
        output_counts[i] = 0;
    
    overlaps.restart();
}

bool PeakMatrix4::nextCol() {
    while (!output_cells.empty()) {
        output_counts[output_cells.back()] = 0;
        output_cells.pop_back();
    }
    if (!loadPeak()) return false;
    return true;
}

uint32_t PeakMatrix4::currentCol() const {
    return output_peak;
}

int32_t PeakMatrix4::load(uint32_t count, SparseVector<uint32_t> buffer) {
    for (int i = 0; i < count; i++) {
        if (output_cells.empty()) return i;
        buffer.idx[i] = output_cells.back();
        buffer.val[i] = output_counts[output_cells.back()];
        output_counts[output_cells.back()] = 0;
        output_cells.pop_back();
    }
    return count;
}

// After loadPeak
// - output_cells and output_counts are populated
bool PeakMatrix4::loadPeak() {
    while (!overlaps.hasCompletedRegions() && overlaps.nextOverlap()) {
        if(peak_buffers[overlaps.region()].empty()) {
            beginPeakBuffer(peak_buffers[overlaps.region()]);
        }
        uint32_t overlap_count = 
            (overlaps.fragStart() >= overlaps.regionStart() && overlaps.fragStart() < overlaps.regionEnd()) +
            (overlaps.fragEnd() > overlaps.regionStart() && overlaps.fragEnd() <= overlaps.regionEnd());
        peak_buffers[overlaps.region()].push_back(Count{overlaps.cell(), overlap_count});
    }
    if (overlaps.hasCompletedRegions()) {
        OverlapIterator::Region completed_region = overlaps.popCompletedRegion();
        
        tallyPeakBuffer(peak_buffers[completed_region.index], completed_region.index);
        clearPeakBuffer(peak_buffers[completed_region.index]);
        return true;
    }
    return false;
}

} // end namespace BPCells