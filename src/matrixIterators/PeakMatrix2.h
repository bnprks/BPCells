#pragma once
#include <algorithm>
#include <cassert>

#include "../fragmentIterators/FragmentsIterator.h"
#include "../fragmentIterators/InsertionsIterator.h"
#include "MatrixIterator.h"

namespace BPCells {

// Output cell x peak matrix (rows = cell_id, col = peak_id)
// Regions are given as half-open format
class PeakMatrix2 : public MatrixLoader<uint32_t> {
public:
    PeakMatrix2(FragmentsLoader &frags, 
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
        
        // counts_buffer [i] = number counts seen so far for cell i
        // cell_buffer = list of cells with non-zero counts so far
        std::vector<uint32_t> counts_buffer, cell_buffer;
        bool active = false;
        bool completed = false;
    public:
        Peak peak;
        // Anticipated usage:
        // 0. Check that isActive returns false, otherwise the buffer is being used for another peak
        // 1. Repeatedly call add_count giving a cell index and count
        // 2. Call finalize once there are no more reads that will overlap the peak
        // 3. Repeatedly call get_count until it returns false
        friend bool operator<(const PeakBuffer &b1, const PeakBuffer &b2) {
            return b1.peak.end < b2.peak.end;
        }

        void reset() {
            assert(!active && completed);
            active = false;
            completed = false;
        }

        void begin() {
            assert(!active && !completed);
            active = true;
            completed = false;
        }

        inline void add_count(uint32_t cell_index, uint32_t count) {
            assert(active);
            if (count == 0) return;
            if (cell_index >= counts_buffer.size()) {
                counts_buffer.resize(cell_index+1, 0);
            }
            if (counts_buffer[cell_index] == 0)
                cell_buffer.push_back(cell_index);
            counts_buffer[cell_index] += count;            
        }
        
        void finalize() {
            assert(active);
            //Profiling says this sort takes 2/3 of computational time on my all-cells benchmark
            //std::sort(cell_buffer.begin(), cell_buffer.end(), std::greater<uint32_t>());
            active = false;
            completed = true;
        }

        inline bool has_count() const {
            assert(completed);
            return cell_buffer.size() > 0; 
        }

        inline void get_count(uint32_t &cell_out, uint32_t &count_out) {
            assert(has_count());
            cell_out = cell_buffer.back();
            cell_buffer.pop_back();
            count_out = counts_buffer[cell_out];
            counts_buffer[cell_out] = 0;
        }

        bool isActive() const {return active; } 
        bool isCompleted() const {return completed; }
    };

    InsertionsIterator insertions;
    const std::vector<std::string> chr_levels;

    std::vector<std::vector<Peak> > sorted_peaks;
    uint32_t current_chr, next_peak; 

    // The range up to n_active_peaks is sorted by end_coord
    std::vector<PeakBuffer> active_peaks;
    uint32_t n_active_peaks;
};


PeakMatrix2::PeakMatrix2(FragmentsLoader &frags, 
        std::vector<uint32_t> &chr, std::vector<uint32_t> &start, std::vector<uint32_t> &end, 
        std::vector<std::string>  &chr_levels) :
        insertions(frags), chr_levels(chr_levels),
        current_chr(0), next_peak(UINT32_MAX), 
        n_active_peaks(0) {

    if (frags.cellCount() < 0)
        throw std::invalid_argument("frags must have a known cell count. Consider using a cell selection to define the number of cells.");
    
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

void PeakMatrix2::restart() {
    current_chr = 0;
    next_peak = UINT32_MAX;

    // Clear out the active peak buffers
    uint32_t cell_out, count_out;
    for (int i = 0; i < n_active_peaks; i++) {
        if (active_peaks[i].isActive()) active_peaks[i].finalize();
        while(active_peaks[i].has_count()) {
            active_peaks[i].get_count(cell_out, count_out);
        }
        active_peaks[i].reset();
    }
    n_active_peaks = 0;
}

bool PeakMatrix2::nextCol() {
    // Clean up the last active peak
    if (n_active_peaks > 0) {
        active_peaks[0].reset();
        
        n_active_peaks -= 1;
        for (int i = 0; i < n_active_peaks; i++) {
            std::swap(active_peaks[i], active_peaks[i+1]);
        }
    }

    if (!loadPeak()) return false;
    return true;
}

uint32_t PeakMatrix2::currentCol() const {
    return active_peaks[0].peak.index;
}

int32_t PeakMatrix2::load(uint32_t count, SparseVector<uint32_t> buffer) {
    if (n_active_peaks == 0) {return 0;} 
    for (int i = 0; i < count; i++) {
        if (!active_peaks[0].has_count()) {
            return i;
        }
        active_peaks[0].get_count(buffer.idx[i], buffer.val[i]);
    }
    return count;
}


// After loadPeak
// - active_peaks[0:n_completed_peaks] are finalized and ready to iterate
bool PeakMatrix2::loadPeak() {
    if(n_active_peaks > 0 && active_peaks[0].isCompleted()) return true;

    if (n_active_peaks == 0 && next_peak >= sorted_peaks[current_chr].size()) {
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
    if (n_active_peaks == 0 && insertions.isSeekable()) {
        insertions.seek(insertions.currentChr(), sorted_peaks[current_chr][next_peak].start);
    }

    while (insertions.nextInsertion()) {
        // Check for any new peaks we've hit
        while (next_peak < sorted_peaks[current_chr].size() &&
                (insertions.coord() >= sorted_peaks[current_chr][next_peak].start)) {
            n_active_peaks += 1;
            if (n_active_peaks > active_peaks.size()) active_peaks.resize(n_active_peaks);
            active_peaks[n_active_peaks - 1].peak = sorted_peaks[current_chr][next_peak];
            active_peaks[n_active_peaks - 1].begin();
            // Insert the peak into position in the sorted order
            for(int i = n_active_peaks - 1; i > 0; i--) {
                if (!(active_peaks[i] < active_peaks[i-1])) break;
                std::swap(active_peaks[i], active_peaks[i-1]);
            }
            next_peak += 1;
        }

        for (uint32_t i = 0; i < n_active_peaks; i++) {
            if (insertions.coord() >= active_peaks[i].peak.end) {
                active_peaks[i].finalize();
                continue;
            }
            uint32_t overlap_count = insertions.coord() >= active_peaks[i].peak.start;
            active_peaks[i].add_count(insertions.cell(), overlap_count);
        }

        // We've completed an active peak, so we're ready to return
        if (n_active_peaks > 0 && active_peaks[0].isCompleted()) return true;
    }

    //We've hit the end of a chromosome, so finalize all active peaks, and mark ready for next chromosome
    for (uint32_t i = 0; i < n_active_peaks; i++) {
        active_peaks[i].finalize();
    }

    next_peak = UINT32_MAX;
    return true;
}

} // end namespace BPCells