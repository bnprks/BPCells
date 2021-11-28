#pragma once
#include <algorithm>
#include <cassert>

#include "../fragmentIterators/FragmentsIterator.h"
#include "MatrixIterator.h"

namespace BPCells {

// Output cell x peak matrix (rows = cell_id, col = peak_id)
// Regions are given as half-open format
class PeakMatrix3 : public MatrixLoader<uint32_t> {
public:
    PeakMatrix3(FragmentsLoader &frags, 
            std::vector<uint32_t> &chr, std::vector<uint32_t> &start, std::vector<uint32_t> &end, 
            std::vector<std::string>  &chr_levels);
    
    // Reset the loader to start from the beginning
    void restart() override;

    // Return the count of rows and columns
    uint32_t rows() const override { return frags.cellCount(); }
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
    class Peak {
    public:
        uint32_t start, end, index;
        friend bool operator<(const Peak &p1, const Peak &p2) {
            return p1.start < p2.start;
        }
    };
    
    class PeakBuffer {
    public:
        //enum class State {Inactive, Active, Completed};

        std::vector<uint32_t> counts_buffer, cell_buffer;
        //State s = State::Inactive;
        Peak peak;

        inline void add_count(uint32_t cell_index, uint32_t count) {
            counts_buffer.push_back(count);
            cell_buffer.push_back(cell_index);
        }
        //inline State getState() {return s;}
        //void setState(State new_state) {s = new_state;}
        void reset(uint32_t max_capacity) {
            if (counts_buffer.capacity() > max_capacity) {
                counts_buffer.resize(max_capacity);
                counts_buffer.shrink_to_fit();
            }
            if (cell_buffer.capacity() > max_capacity) {
                cell_buffer.resize(max_capacity);
                cell_buffer.shrink_to_fit();
            }
            cell_buffer.resize(0);
            counts_buffer.resize(0);
            //s = State::Inactive;
        }

    };
    FragmentsIterator frags;
    const std::vector<std::string> chr_levels;

    std::vector<std::vector<Peak> > sorted_peaks;
    uint32_t current_chr, next_peak; 

    std::vector<PeakBuffer> active_peaks;
    uint32_t n_active_peaks, n_completed_peaks;

    uint32_t output_peak;
    std::vector<uint32_t> output_counts; // Length cells, with counts to output
    std::vector<uint32_t> output_cells; // List of non-zero indices in output_counts
    
    // Load the next peak and populate the output_counts + output_cells buffers
    bool loadPeak();

    // Transfer the counts from the first buffer in active_peaks to
    // the output_counts and output_cells vectors,
    // then reset the peak buffer, and reorder active_peaks & update
    // n_active_peaks and n_completed_peaks accordingly
    void tallyFirstBuffer() {
        PeakBuffer &buf = active_peaks[0];
        assert(buf.counts_buffer.size() == buf.cell_buffer.size());

        for(size_t i = 0; i < buf.counts_buffer.size(); i++) {
            uint32_t cell = buf.cell_buffer[i];
            uint32_t count = buf.counts_buffer[i];
            
            if (output_counts[cell] == 0)
                output_cells.push_back(cell);
            output_counts[cell] += count;
        }

        output_peak = buf.peak.index;

        buf.reset(1024); // Try resetting to a fixed size here, rather than adaptive in number of cells, not sure if good or bad idea
        assert(n_active_peaks > 0 && n_completed_peaks > 0);
        // Swap so that the just finished peak goes to the end of the active_peaks vector,
        // and if there are more completed peaks we will have another completed peak at the
        // front of the vector as well
        std::swap(active_peaks[0], active_peaks[n_active_peaks - 1]);
        if (n_active_peaks > n_completed_peaks)
            std::swap(active_peaks[0], active_peaks[n_completed_peaks - 1]);
        if (n_active_peaks > 0)
        n_completed_peaks -= 1;
        n_active_peaks -= 1;
    }
    
};


PeakMatrix3::PeakMatrix3(FragmentsLoader &frags, 
        std::vector<uint32_t> &chr, std::vector<uint32_t> &start, std::vector<uint32_t> &end, 
        std::vector<std::string>  &chr_levels) :
        frags(frags), chr_levels(chr_levels),
        current_chr(0), next_peak(UINT32_MAX), 
        n_active_peaks(0), n_completed_peaks(0) {

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

void PeakMatrix3::restart() {
    current_chr = 0;
    next_peak = UINT32_MAX;

    // Clear out the active peak buffers
    for (int i = 0; i < n_active_peaks; i++) {
        active_peaks[i].reset(1024); // Reset down to a constant, manageable size
    }
    n_active_peaks = 0;
    n_completed_peaks = 0;
    output_cells.resize(0);
    for (size_t i = 0; i < output_counts.size(); i++)
        output_counts[i] = 0;
}

bool PeakMatrix3::nextCol() {
    while (!output_cells.empty()) {
        output_counts[output_cells.back()] = 0;
        output_cells.pop_back();
    }
    if (!loadPeak()) return false;
    return true;
}

uint32_t PeakMatrix3::currentCol() const {
    return output_peak;
}

int32_t PeakMatrix3::load(uint32_t count, SparseVector<uint32_t> buffer) {
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
// - active_peaks[0:n_completed_peaks] are finalized and ready to iterate
bool PeakMatrix3::loadPeak() {
    if(n_completed_peaks > 0) {
        tallyFirstBuffer();
        return true;
    }

    if (n_active_peaks == 0 && next_peak >= sorted_peaks[current_chr].size()) {
        // Continue to the next chromosome for which we have peaks.
        // Match based on chromosome name rather than chromosome ID.
        while(true) {
            if (!frags.nextChr()) return false;
            std::string chr_name = frags.chrNames(frags.currentChr());
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
    if (n_active_peaks == 0 && frags.isSeekable()) {
        frags.seek(frags.currentChr(), sorted_peaks[current_chr][next_peak].start);
    }

    while (frags.nextFrag()) {
        // Check for any new peaks we've hit
        while (next_peak < sorted_peaks[current_chr].size() &&
                (frags.end() >= sorted_peaks[current_chr][next_peak].start ||
                frags.start() >= sorted_peaks[current_chr][next_peak].start)) {
            n_active_peaks += 1;
            if (n_active_peaks > active_peaks.size()) active_peaks.resize(n_active_peaks);
            active_peaks[n_active_peaks - 1].peak = sorted_peaks[current_chr][next_peak];
            next_peak += 1;
        }

        for (uint32_t i = 0; i < n_active_peaks; i++) {
            if (frags.start() >= active_peaks[i].peak.end) {
                std::swap(active_peaks[n_completed_peaks], active_peaks[i]);
                n_completed_peaks += 1;
                continue;
            }
            uint32_t overlap_count = 
                (frags.start() >= active_peaks[i].peak.start && frags.start() < active_peaks[i].peak.end) +
                (frags.end() > active_peaks[i].peak.start && frags.end() <= active_peaks[i].peak.end);
            if (overlap_count)
                active_peaks[i].add_count(frags.cell(), overlap_count);
        }
        // We've completed an active peak, so we're ready to return
        if (n_completed_peaks) {
            tallyFirstBuffer();
            return true;
        }
    }
    //We've hit the end of a chromosome, so finalize all active peaks, and mark ready for next chromosome
    for (uint32_t i = 0; i < n_active_peaks; i++) {
        n_completed_peaks += 1;
    }

    next_peak = UINT32_MAX;
    if (n_completed_peaks > 0) {
        tallyFirstBuffer();
        return true;
    } else {
        // We get here if we tried to seek to a peak which was past the last
        // fragment in the data
        return loadPeak();
    }
}

} // end namespace BPCells