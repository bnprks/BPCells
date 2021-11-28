#pragma once
#include <algorithm>
#include <cassert>

#include "FragmentsIterator.h"
#include "../matrixIterators/MatrixIterator.h"

namespace BPCells {

// Helper class for iterating over sorted fragments overlapping a given set of genomic regions.
// Usage:
// while (it.nextOverlap()) {
//      uint32_t frag_start = it.fragStart();
//      uint32_t frag_end = it.fragEnd();
//      uint32_t cell_id = it.cell();
//      uint32_t chr = it.chr();
//      uint32_t region_id = it.region();
//      uint32_t region_start = it.regionStart();
//      uint32_t region_end = it.regionEnd();
// }
// 
class OverlapIterator {
public:
    class Region {
    public:
        uint32_t start, end, index;
        friend bool operator<(const Region &r1, const Region &r2) {
            return r1.start < r2.start;
        }
    };

    using RegionList = std::vector<std::vector<OverlapIterator::Region>>;

    // High-level constructor: pass unsorted chromosome ID, start, end triplets
    OverlapIterator(FragmentsLoader &frags, 
            std::vector<uint32_t> &chr, std::vector<uint32_t> &start, std::vector<uint32_t> &end, 
            std::vector<std::string>  &chr_levels);

    // Low-level constructor: pass pre-sorted regions
    OverlapIterator(FragmentsLoader &frags, 
            std::vector<std::vector<Region>> sorted_regions, 
            std::vector<std::string>  &chr_levels);
    
    // Reset the iterator to start from the beginning
    void restart();

    // Advance through overlaps
    inline bool nextOverlap();

    // Access information on the current overlapping fragment
    inline uint32_t fragStart() const {return frags.start();}
    inline uint32_t fragEnd() const {return frags.end(); }
    inline uint32_t cell() const {return frags.cell(); }
    inline uint32_t chr() const {return current_chr; }
    inline uint32_t region() const {return active_regions[current_active_region].index; }
    inline uint32_t regionStart() const {return active_regions[current_active_region].start; }
    inline uint32_t regionEnd() const {return active_regions[current_active_region].end; }

    // Check for regions that were completed
    inline bool hasCompletedRegions() {return !completed_regions.empty(); }
    inline Region popCompletedRegion() {
        Region r = completed_regions.back();
        completed_regions.pop_back();
        return r;
    }

    static inline std::vector<std::vector<Region>> sortRegions(uint32_t n_chrs, 
        std::vector<uint32_t> &chr, std::vector<uint32_t> &start, std::vector<uint32_t> &end);
private:
    
    
    FragmentsIterator frags;
    const std::vector<std::string> chr_levels;

    std::vector<std::vector<Region> > sorted_regions;
    uint32_t current_chr, next_region; 

    std::vector<Region> active_regions;
    std::vector<Region> completed_regions;


    uint32_t current_active_region;
    bool frags_finished = false;

    // Load the next chromosome in frags until we get to a 
    // chromosome where we have regions listed. 
    // Postconditions:
    // - current_chr is updated to match the index of current chromosome in
    //   sorted_regions
    // - next_region is 0
    // - frags.nextFrag() will get the first fragment of the next chromosome
    inline bool loadNextChr();
    
};

// High-level constructor: pass unsorted chromosome ID, start, end triplets
OverlapIterator::OverlapIterator(FragmentsLoader &frags, 
        std::vector<uint32_t> &chr, std::vector<uint32_t> &start, std::vector<uint32_t> &end, 
        std::vector<std::string>  &chr_levels) : 
        OverlapIterator(frags, sortRegions(chr_levels.size(), chr, start, end), chr_levels) {}

OverlapIterator::RegionList OverlapIterator::sortRegions(
    uint32_t n_chrs, std::vector<uint32_t> &chr, std::vector<uint32_t> &start, std::vector<uint32_t> &end) {

    if (chr.size() != start.size() || chr.size() != end.size())
    throw std::invalid_argument("chr, start, and end must all be same length");

    RegionList ret(n_chrs);
    for (size_t i = 0; i < chr.size(); i++) {
        if (chr[i] >= n_chrs)
            throw std::invalid_argument("chr has values higher than the number of chromosomes");
        Region r;
        r.start = start[i];
        r.end = end[i];
        r.index = i;
        ret[chr[i]].push_back(r);
    }
    for (size_t i = 0; i < ret.size(); i++) {
        std::sort(ret[i].begin(), ret[i].end());
    }
    return ret;
}

// Low-level constructor: pass pre-sorted regions
OverlapIterator::OverlapIterator(FragmentsLoader &frags, 
        RegionList sorted_regions, 
        std::vector<std::string>  &chr_levels) :
        frags(frags), chr_levels(chr_levels), sorted_regions(sorted_regions),
        current_chr(0), next_region(UINT32_MAX), current_active_region(0) {
    
    if (sorted_regions.size() != chr_levels.size())
        throw std::invalid_argument("sorted_regions must have same length as chr_levels");
    
    for (size_t i = 0; i < sorted_regions.size(); i++) {
        if (!std::is_sorted(sorted_regions[i].begin(), sorted_regions[i].end()))
            throw std::invalid_argument("sorted_regions must be sorted by start coordinate");
    }

}

inline void OverlapIterator::restart() {
    current_chr = 0;
    next_region = UINT32_MAX;
    active_regions.resize(0);
}

// Load the next chromosome in frags until we get to a 
// chromosome where we have regions listed. 
// Postconditions:
// - current_chr is updated to match the index of current chromosome in
//   sorted_regions
// - next_region is 0
// - frags.nextFrag() will get the first fragment of the next chromosome
// - active_regions is empty, and any remaining regions are appended to completed_regions
inline bool OverlapIterator::loadNextChr() {
    while (!active_regions.empty()) {
        completed_regions.push_back(std::move(active_regions.back()));
        active_regions.pop_back();
    }
    while(true) {
        if (!frags.nextChr()) {
            frags_finished = true;
            return false;
        }
        std::string chr_name = frags.chrNames(frags.currentChr());
        auto idx = std::find(chr_levels.begin(), chr_levels.end(), chr_name);
        if (idx != chr_levels.end()) {
            current_chr = idx - chr_levels.begin();
            if (sorted_regions[current_chr].size() > 0)
                break;
        }
    }
    next_region = 0;
    return true;
}

// Find the next fragment overlapping a region
// Returns false if there are no more overlaps with the input fragments
// Postconditions:
// - active_regions[current_active_region] overlaps with region [frags.start(), frags.end())
inline bool OverlapIterator::nextOverlap() {
    if (frags_finished) return false;
    current_active_region++;
    
    while (true) {
        // Process current fragment, removing competed regions and breaking
        // when we find an overlap
        while (current_active_region < active_regions.size()) {
            Region &current_region = active_regions[current_active_region];
            if (frags.start() >= current_region.end) {
                // Remove the active region and continue without incrementing
                // current_active_region
                std::swap(current_region, active_regions.back());
                completed_regions.push_back(std::move(active_regions.back()));
                active_regions.pop_back();
                continue;
            } 
            // Check for overlap and break if there is one
            bool overlap = 
                (frags.start() >= current_region.start && 
                 frags.start() < current_region.end) ||
                (frags.end() > current_region.start && 
                 frags.end() <= current_region.end);
            
            if (overlap) { 
                return true;
            }

            current_active_region++;
        }
        current_active_region = 0;
        
        if (active_regions.empty()) {
            // If we're done with regions on this chromosome, advance to next chromosome
            if (next_region >= sorted_regions[current_chr].size()) {
                if(!loadNextChr()) return false;
            }
            // Seek to next region and iterate until we hit the next region
            if (frags.isSeekable()) {     
                frags.seek(frags.currentChr(), sorted_regions[current_chr][next_region].start);
            }
            while (frags.nextFrag()) {
                if (frags.end() >= sorted_regions[current_chr][next_region].start) break;
            }
        } else {
            // Load next fragment, and if we hit end of chromosome,
            // go to next chromosome with regions
            while (!frags.nextFrag()) {
                if (!loadNextChr()) return false;
            }
        }
        
        // Check for any new peaks that need to be added for new fragment
        while (next_region < sorted_regions[current_chr].size() &&
                (frags.end() >= sorted_regions[current_chr][next_region].start ||
                frags.start() >= sorted_regions[current_chr][next_region].start)) {
            active_regions.push_back(sorted_regions[current_chr][next_region]);
            next_region += 1;
        }
    }
}

} // end namespace BPCells