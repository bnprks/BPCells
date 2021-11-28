#pragma once
#include <algorithm>
#include <cassert>

#include "../fragmentIterators/FragmentsIterator.h"
#include "../fragmentIterators/InsertionsIterator.h"
#include "MatrixIterator.h"

namespace BPCells {

// Output cell x peak matrix (rows = cell_id, col = peak_id)
// Regions are given as half-open format
class PeakMatrix5 : public MatrixLoader<uint32_t> {
public:
    PeakMatrix5(FragmentsLoader &frags, 
            std::vector<uint32_t> &chr, std::vector<uint32_t> &start, std::vector<uint32_t> &end, 
            std::vector<std::string>  &chr_levels);
    
    // Reset the loader to start from the beginning
    void restart() override;

    // Return the count of rows and columns
    uint32_t rows() const override { return insertions.cellCount(); }
    uint32_t cols() const override { 
        uint32_t count = 0;
        for (auto peakVec : sorted_peaks)
            count += peakVec.size();
        return count;
    }

    bool nextCol() override;
    uint32_t currentCol() const override;

    int32_t load(uint32_t count, SparseVector<uint32_t> buffer) override;
private:
    bool loadPeak();
    
    class Peak {
    public:
        uint32_t start, end, index;
        friend bool operator<(const Peak &p1, const Peak &p2) {
            return p1.start < p2.start;
        }
    };

    class PeakBuffer {
    public:
        std::vector<uint32_t> cell_buffer;
        Peak peak;

        inline void add_count(uint32_t cell_index) {
            cell_buffer.push_back(cell_index);
        }

        void reset(uint32_t max_capacity) {
            if (cell_buffer.capacity() > max_capacity) {
                cell_buffer.resize(max_capacity);
                cell_buffer.shrink_to_fit();
            }
            cell_buffer.resize(0);
        }

    };


    InsertionsIterator insertions;
    const std::vector<std::string> chr_levels;

    std::vector<std::vector<Peak> > sorted_peaks;
    uint32_t current_chr, next_peak; 

    // active_peaks is sorted by end_coord in descending order,
    // such that the next peak to tally is always at the end
    std::vector<PeakBuffer> active_peaks;
    std::vector<PeakBuffer> inactive_buffers;
    uint32_t n_completed_peaks;

    uint32_t output_peak;
    std::vector<uint32_t> output_counts; // Length cells, with counts to output
    std::vector<uint32_t> output_cells; // List of non-zero indices in output_counts

    void beginPeakBuffer(PeakBuffer &buf) {
        if (inactive_buffers.empty()) return;
        std::swap(buf, inactive_buffers.back());
        //inactive_buffers.pop_back();
    }
    void tallyPeakBuffer(PeakBuffer &buf, uint32_t output_idx) {
        for(auto c : buf.cell_buffer) {
            if (output_counts[c] == 0)
                output_cells.push_back(c);
            output_counts[c] += 1;
        }

        output_peak = output_idx;
    };
    void clearPeakBuffer(PeakBuffer &buf) {
        buf.reset(1024);
        if (buf.cell_buffer.capacity() > 0) {
            inactive_buffers.push_back(std::move(buf));
        }
    };


};


PeakMatrix5::PeakMatrix5(FragmentsLoader &frags, 
        std::vector<uint32_t> &chr, std::vector<uint32_t> &start, std::vector<uint32_t> &end, 
        std::vector<std::string>  &chr_levels) :
        insertions(frags), chr_levels(chr_levels),
        current_chr(0), next_peak(UINT32_MAX), 
        n_completed_peaks(0) {

    if (frags.cellCount() < 0)
        throw std::invalid_argument("frags must have a known cell count. Consider using a cell selection to define the number of cells.");
    
    output_counts.resize(frags.cellCount(), 0);

    if (chr.size() != start.size() || chr.size() != end.size())
        throw std::invalid_argument("chr, start, and end must all be same length");
    
    sorted_peaks.resize(chr_levels.size());
    for (size_t i = 0; i < chr.size(); i++) {
        if (chr[i] >= chr_levels.size())
            throw std::invalid_argument("chr has values higher than length of chr_levels");
        Peak p;
        p.start = start[i];
        p.end = end[i];
        p.index = i;
        sorted_peaks[chr[i]].push_back(p);
    }
    for (size_t i = 0; i < sorted_peaks.size(); i++) {
        std::sort(sorted_peaks[i].begin(), sorted_peaks[i].end());
    }
}

void PeakMatrix5::restart() {
    current_chr = 0;
    next_peak = UINT32_MAX;

    // Clear out the active peak buffers
    for (auto &buf : active_peaks) {
        clearPeakBuffer(buf);
    }
    active_peaks.resize(0);

    n_completed_peaks = 0;
    output_cells.resize(0);
    for (size_t i = 0; i < output_counts.size(); i++)
        output_counts[i] = 0;
}

bool PeakMatrix5::nextCol() {
    while (!output_cells.empty()) {
        output_counts[output_cells.back()] = 0;
        output_cells.pop_back();
    }
    if (!loadPeak()) return false;

    //printf("Outputting peak %d\n", active_peaks.back().peak.index);
    tallyPeakBuffer(active_peaks.back(), active_peaks.back().peak.index);
    clearPeakBuffer(active_peaks.back());
    active_peaks.pop_back();
    n_completed_peaks -= 1;

    return true;
}

uint32_t PeakMatrix5::currentCol() const {
    return output_peak;
}

int32_t PeakMatrix5::load(uint32_t count, SparseVector<uint32_t> buffer) {
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
// - The last n_completed_peaks entries of active_peaks are reeady to be tallied
bool PeakMatrix5::loadPeak() {
    if (n_completed_peaks > 0) {
        return true;
    }

    if (active_peaks.empty() && next_peak >= sorted_peaks[current_chr].size()) {
        //printf("Advancing to next chr\n");
        // Continue to the next chromosome for which we have peaks.
        // Match based on chromosome name rather than chromosome ID.
        while(true) {
            if (!insertions.nextChr()) return false;
            std::string chr_name = insertions.chrNames(insertions.currentChr());
            auto idx = std::find(chr_levels.begin(), chr_levels.end(), chr_name);
            if (idx != chr_levels.end()) {
                current_chr = idx - chr_levels.begin();
                if (sorted_peaks[current_chr].size() > 0)
                    break;
            }
        }
        next_peak = 0;
    }

    // If there are no active peaks, try to seek to just before the next peek begins
    if (active_peaks.empty() && insertions.isSeekable()) {
        insertions.seek(insertions.currentChr(), sorted_peaks[current_chr][next_peak].start);
    }

    while (insertions.nextInsertion()) {
        // Check for any new peaks we've hit
        //printf("Insertion at %d, compared to next_peak at %s:%d-%d\n",
        //     insertions.coord(),
        //     insertions.chrNames(current_chr),
        //     sorted_peaks[current_chr][next_peak].start,
        //     sorted_peaks[current_chr][next_peak].end
        // );
        while (next_peak < sorted_peaks[current_chr].size() &&
                (insertions.coord() >= sorted_peaks[current_chr][next_peak].start)) {
            PeakBuffer b;
            b.peak = sorted_peaks[current_chr][next_peak];
            active_peaks.push_back(std::move(b));
            
            // Insert the peak into position in the sorted order
            for(int i = active_peaks.size() - 1; i > 0; i--) {
                if (active_peaks[i].peak.end <= active_peaks[i-1].peak.end) break;
                std::swap(active_peaks[i], active_peaks[i-1]);
            }
            next_peak += 1;
        }

        for (auto &buf : active_peaks) {
            if (insertions.coord() >= buf.peak.end) {
                n_completed_peaks += 1;
                continue;
            }
            if (insertions.coord() >= buf.peak.start) {
                buf.add_count(insertions.cell());
            }
        }

        // We've completed an active peak, so we're ready to return
        if (n_completed_peaks > 0) return true;
    }

    //We've hit the end of a chromosome, so finalize all active peaks, and mark ready for next chromosome
    n_completed_peaks += active_peaks.size();

    next_peak = UINT32_MAX;

    return loadPeak();
}

} // end namespace BPCells