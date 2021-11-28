#pragma once
#include <algorithm>
#include <cassert>

#include "../fragmentIterators/FragmentsIterator.h"
#include "../fragmentIterators/InsertionsIterator2.h"
#include "MatrixIterator.h"

namespace BPCells {

// Output cell x peak matrix (rows = cell_id, col = peak_id)
// Regions are given as half-open format
class PeakMatrix6 : public MatrixLoader<uint32_t> {
public:
    PeakMatrix6(FragmentsLoader &frags, 
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
    };


    InsertionsIterator2 insertions;
    const std::vector<std::string> chr_levels;

    std::vector<std::vector<Peak> > sorted_peaks;
    uint32_t current_chr, next_peak; 

    // active_peaks is sorted by end_coord in descending order,
    // such that the next peak to tally is always at the end
    std::vector<PeakBuffer> active_bufs, inactive_bufs, finished_bufs;

    uint32_t output_peak;
    std::vector<uint32_t> output_counts; // Length cells, with counts to output
    std::vector<uint32_t> output_cells; // List of non-zero indices in output_counts


    size_t max_buffer_capacity = 0;
    void printBuffers() {
        printf("# inactive=%zu, active=", inactive_bufs.size());
        for (auto b : active_bufs) {
            printf("%d ", b.peak.index);
        }
        printf(" finished=");
        for (auto b : finished_bufs) {
            printf("%d ", b.peak.index);
        }
        printf("\n");
    }
    // Add a PeakBuffer for the peak p to active_bufs
    void beginPeakBuffer(Peak p) {
        //printf("Creating new peak buffer on peak %d ", p.index);
        //printBuffers();
        if (inactive_bufs.empty()) {
            printf("ALLOCATING NEW\n");
            active_bufs.emplace_back();
        } else {
            active_bufs.push_back(std::move(inactive_bufs.back()));
            inactive_bufs.pop_back();
        }
        active_bufs.back().peak = p;
        //printBuffers();
    }
    // Move a buffer from active_bufs to finished_bufs
    void finalizePeakBuffer(PeakBuffer &buf) {
        //printf("Finalizing peak buffer on peak %d ", buf.peak.index);
        //printBuffers();
        //printf("Pre-swap, active_bufs.back = %p, %d,%d-%d ", &active_bufs.back().cell_buffer[0], active_bufs.back().peak.index, active_bufs.back().peak.start, active_bufs.back().peak.end);
        std::swap(buf, active_bufs.back());
        //printf("post-swap, active_bufs.back = %p, %d,%d-%d\n", &active_bufs.back().cell_buffer[0], active_bufs.back().peak.index, active_bufs.back().peak.start, active_bufs.back().peak.end);
        finished_bufs.push_back(std::move(active_bufs.back()));
        active_bufs.pop_back();
        //printBuffers();
    }

    // Tally the results from the back of finished_buf into output buffers
    void tallyPeakBuffer() {
        //printf("Tallying peak buffer on peak %d ", finished_bufs.back().peak.index);
        //printBuffers();

        for(auto c : finished_bufs.back().cell_buffer) {
            if (output_counts[c] == 0)
                output_cells.push_back(c);
            output_counts[c] += 1;
        }
        finished_bufs.back().cell_buffer.resize(0);
        output_peak = finished_bufs.back().peak.index;

        inactive_bufs.push_back(std::move(finished_bufs.back()));
        finished_bufs.pop_back();
        //printBuffers();
    };
};


PeakMatrix6::PeakMatrix6(FragmentsLoader &frags, 
        std::vector<uint32_t> &chr, std::vector<uint32_t> &start, std::vector<uint32_t> &end, 
        std::vector<std::string>  &chr_levels) :
        insertions(frags), chr_levels(chr_levels),
        current_chr(0), next_peak(UINT32_MAX) {

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
        sorted_peaks[i].push_back({UINT32_MAX, UINT32_MAX, UINT32_MAX});
    }
}

void PeakMatrix6::restart() {
    current_chr = 0;
    next_peak = UINT32_MAX;

    // Clear out the active and finished peak buffers
    // This could potentially be improved by just moving
    // the peakbuffers to inactive rather than deleting them entirely
    active_bufs.resize(0);
    finished_bufs.resize(0);

    output_cells.resize(0);
    for (size_t i = 0; i < output_counts.size(); i++)
        output_counts[i] = 0;
}

bool PeakMatrix6::nextCol() {
    while (!output_cells.empty()) {
        output_counts[output_cells.back()] = 0;
        output_cells.pop_back();
    }
    if (!loadPeak()) return false;

    tallyPeakBuffer();

    //if (output_peak % 10000 == 0) printf("nextCol output_peak = %d", output_peak);

    return true;
}

uint32_t PeakMatrix6::currentCol() const {
    return output_peak;
}

int32_t PeakMatrix6::load(uint32_t count, SparseVector<uint32_t> buffer) {
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
bool PeakMatrix6::loadPeak() {
    //printf("Load Peak 1\n");
    size_t buffer_capacity = 0;
    for (auto b : inactive_bufs) {
        buffer_capacity += b.cell_buffer.capacity();
    }
    for (auto b : active_bufs) {
        buffer_capacity += b.cell_buffer.capacity();
    }
    for (auto b : finished_bufs) {
        buffer_capacity += b.cell_buffer.capacity();
    }
    if (buffer_capacity > max_buffer_capacity) {
        max_buffer_capacity = buffer_capacity;
        printf("New max capacity: %zu\n", buffer_capacity);
    }

    if (!finished_bufs.empty()) {
        return true;
    }
    
    // FOR DEBUGGING
    //if (next_peak > 5000 && next_peak != UINT32_MAX) return false;

    if (active_bufs.empty()) {
        if(next_peak >= sorted_peaks[current_chr].size() || 
           sorted_peaks[current_chr][next_peak].start == UINT32_MAX) {
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
        if (insertions.isSeekable()) {
            insertions.seek(insertions.currentChr(), sorted_peaks[current_chr][next_peak].start);
        }   
    }

    while (insertions.nextInsertion()) {
        while (insertions.coord() >= sorted_peaks[current_chr][next_peak].start) {
            if (insertions.coord() < sorted_peaks[current_chr][next_peak].end) {
                beginPeakBuffer(sorted_peaks[current_chr][next_peak]);
            }
            next_peak += 1;
        }

        size_t i = 0;
        while (i < active_bufs.size()) {
            if (insertions.coord() >= active_bufs[i].peak.end) {
                finalizePeakBuffer(active_bufs[i]);
                continue;
            }
            if (insertions.coord() >= active_bufs[i].peak.start) {
                active_bufs[i].add_count(insertions.cell());
            }
            i++;
        }

        // We've completed an active peak, so we're ready to return
        if (!finished_bufs.empty()) return true;
    }

    //We've hit the end of a chromosome, so finalize all active peaks, and mark ready for next chromosome
    while (!active_bufs.empty()) finalizePeakBuffer(active_bufs.back());

    next_peak = UINT32_MAX;

    return loadPeak();
}

} // end namespace BPCells