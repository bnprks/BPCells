#pragma once
#include <algorithm>
#include <cassert>

#include "../fragmentIterators/FragmentIterator.h"
#include "../fragmentIterators/OverlapIterator.h"
#include "MatrixIterator.h"

namespace BPCells {

// Output cell x tile matrix (rows = cell_id, col = peak_id)
// Regions are given as half-open format, with a tiling width.
// Each region is divided into ceiling((end-start)/tile_width) tiles, which
// are output in order. (So if we have 200 tiles per region, region 1 will have
// columns 0-199, region 2 will have columns 200-399, etc.). Note that
// tile widths do not need to be all identical.
class TileMatrix : public MatrixLoader<uint32_t> {
public:
    TileMatrix(FragmentsLoader &frags, 
            std::vector<uint32_t> &chr, std::vector<uint32_t> &start, std::vector<uint32_t> &end, 
            std::vector<uint32_t> &tile_widths,
            std::vector<std::string>  &chr_levels);
    
    // Reset the loader to start from the beginning
    void restart() override;

    // Return the count of rows and columns
    uint32_t rows() const override { return frags.cellCount(); }
    uint32_t cols() const override { return total_tiles; }

    bool nextCol() override;
    uint32_t currentCol() const override;

    int32_t load(uint32_t count, SparseVector<uint32_t> buffer) override;
private:
    class Count {
    public:
        uint32_t cell; // It's assumed the counts are all 1, so no need to store explicitly
    };
    using CountVector = std::vector<Count>;

    class RegionInfo {
    public:
        uint32_t output_index; // Column index of 0th tile output
        uint32_t tile_width; // width of each tile
        uint32_t next_tile = 0; // Next tile index to output
        uint32_t completed_tiles = 0; // Number of tiles which have been completed
    };

    std::vector<CountVector> inactive_tile_buffers; // Object pool of vectors with non-zero capacity
    std::vector<std::vector<CountVector>> inactive_region_buffers; // Object pool of vectors with non-zero capacity

    std::vector<std::vector<CountVector>> region_buffers; // list of tile buffers for each region, sorted by region coordinate
    
    FragmentsLoader &frags;
    OverlapIterator overlaps;

    std::vector<RegionInfo> region_info;

    uint32_t output_col;
    std::vector<uint32_t> output_counts; // Length cells, with counts to output
    std::vector<uint32_t> output_cells; // List of non-zero indices in output_counts
    
    uint32_t total_tiles;
    // Index of a completed region for which multiple tiles need to be output, 
    // or UINT32_MAX when no region is currently finished
    uint32_t completed_region = UINT32_MAX; 
    
    
    // Load the next tile and populate the output_counts + output_cells buffers
    bool loadTile();

    void beginTileBuffer(CountVector &buf) {
        if (inactive_tile_buffers.empty()) return;
        std::swap(buf, inactive_tile_buffers.back());
        inactive_tile_buffers.pop_back();
    }
    void tallyTileBuffer(CountVector &buf, uint32_t output_idx) {
        for(Count c : buf) {
            if (output_counts[c.cell] == 0)
                output_cells.push_back(c.cell);
            output_counts[c.cell] += 1;
        }

        output_col = output_idx;
    };
    void clearTileBuffer(CountVector &buf) {
        buf.resize(0);
        if (buf.capacity() > 0) {
                // if (buf.capacity() > 128) {
                //     // Try freeing memory if we have a large buffer to stores
                //     buf.resize(128);
                //     buf.shrink_to_fit();
                //     buf.resize(0);
                // }
            inactive_tile_buffers.push_back(std::move(buf));
            buf = CountVector();
        }
    };

    void beginRegionBuffer(std::vector<CountVector> &buf) {
        if (inactive_region_buffers.empty()) return;
        std::swap(buf, inactive_region_buffers.back());
        inactive_region_buffers.pop_back();
    };
    void clearRegionBuffer(std::vector<CountVector> &buf) {
        for (auto vec : buf)
            clearTileBuffer(vec);
        buf.resize(0);
        if (buf.capacity() > 0) {
            inactive_region_buffers.push_back(std::move(buf));
            buf = std::vector<CountVector>();
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


TileMatrix::TileMatrix(FragmentsLoader &frags, 
        std::vector<uint32_t> &chr, std::vector<uint32_t> &start, std::vector<uint32_t> &end, 
        std::vector<uint32_t> &tile_widths,
        std::vector<std::string>  &chr_levels) :
        region_buffers(chr.size()), 
        frags(frags),
        overlaps(frags, reindexSortRegions(chr_levels.size(), chr, start, end), chr_levels),
        region_info(chr.size()) {

    if (frags.cellCount() < 0)
        throw std::invalid_argument("frags must have a known cell count. Consider using a cell selection to define the number of cells.");
    
    if (chr.size() != tile_widths.size())
        throw std::invalid_argument("tile_widths must be same length as chr, start, and end");


    output_counts.resize(frags.cellCount(), 0);
    OverlapIterator::RegionList sorted_regions = 
        OverlapIterator::sortRegions(chr_levels.size(), chr, start, end);

    std::vector<uint32_t> output_index(chr.size());
    total_tiles = 0;
    for (size_t i = 0; i < tile_widths.size(); i++) {
        output_index[i] = total_tiles;

        uint32_t region_width = end[i] - start[i];
        uint32_t tile_width = tile_widths[i];
        uint32_t tile_count = (region_width + tile_width - 1) / tile_width;

        total_tiles += tile_count;
    }
    
    size_t region_index = 0;
    for (size_t i = 0; i < sorted_regions.size(); i++) {
        for (size_t j = 0; j < sorted_regions[i].size(); j++) {
            region_info[region_index].output_index = output_index[sorted_regions[i][j].index];
            region_info[region_index].tile_width = tile_widths[sorted_regions[i][j].index];
            //sorted_regions[i][j].index = peak_index;
            region_index++;
        }
    }
    //printf("Completed TileMatrix constructor. total_tiles=%d\n", total_tiles);
}

void TileMatrix::restart() {
    // Clear out the active peak buffers
    for (auto b : region_buffers) {
        clearRegionBuffer(b);
    }

    output_cells.resize(0);
    for (size_t i = 0; i < output_counts.size(); i++)
        output_counts[i] = 0;
    
    overlaps.restart();
}

bool TileMatrix::nextCol() {
    while (!output_cells.empty()) {
        output_counts[output_cells.back()] = 0;
        output_cells.pop_back();
    }
    if (!loadTile()) return false;
    //printf("Loaded tile, column = %d\n", output_col);
    return true;
}

uint32_t TileMatrix::currentCol() const {
    return output_col;
}

int32_t TileMatrix::load(uint32_t count, SparseVector<uint32_t> buffer) {
    for (int i = 0; i < count; i++) {
        if (output_cells.empty()) return i;
        buffer.idx[i] = output_cells.back();
        buffer.val[i] = output_counts[output_cells.back()];
        output_counts[output_cells.back()] = 0;
        output_cells.pop_back();
    }
    return count;
}

// After loadTile
// - output_cells and output_counts are populated
bool TileMatrix::loadTile() {
    //printf("LoadTile 1\n");
    // Check if we already have a region with tiles ready to output
    if (completed_region != UINT32_MAX) {
        RegionInfo &info = region_info[completed_region];
        // No more tiles to output from this region right now
        if (info.next_tile >= info.completed_tiles) {
            if (info.next_tile >= region_buffers[completed_region].size())
                clearRegionBuffer(region_buffers[completed_region]);
        } else {
            while (info.next_tile < info.completed_tiles) {
                // Only ouptut for tiles with counts
                if (region_buffers[completed_region][info.next_tile].empty()) {
                    info.next_tile++;
                    continue;
                }
                //printf("Outputting index %d, tile %d\n", completed_region, info.next_tile);
                // Output tile from this region
                tallyTileBuffer(region_buffers[completed_region][info.next_tile], info.output_index + info.next_tile);
                info.next_tile++;
                return true;
            }
            
        }
    }
    if (overlaps.hasCompletedRegions()) {
        OverlapIterator::Region finished_region = overlaps.popCompletedRegion();
        completed_region = finished_region.index;
        region_info[completed_region].completed_tiles = region_buffers[completed_region].size();
        return loadTile();
    }
    // Look for new overlaps
    while (!overlaps.hasCompletedRegions() && overlaps.nextOverlap()) {
        RegionInfo &info = region_info[overlaps.region()];

        // Start of a new region
        if(region_buffers[overlaps.region()].empty()) {
            beginRegionBuffer(region_buffers[overlaps.region()]);
            uint32_t region_width = overlaps.regionEnd() - overlaps.regionStart();
            uint32_t tile_count = (region_width + info.tile_width - 1) / info.tile_width;
            region_buffers[overlaps.region()].resize(tile_count);
            //printf("Starting new region %s:%d-%d (%dbp)\n", frags.chrNames(frags.currentChr()), overlaps.regionStart(), overlaps.regionEnd(), info.tile_width);
        }

        // Calculate overlap tiles

        // printf("Fragment %s:%d-%d, Region %s:%d-%d (%dbp) ",
        //     frags.cellNames(overlaps.cell()), overlaps.fragStart(), overlaps.fragEnd(),
        //     frags.chrNames(frags.currentChr()), overlaps.regionStart(), overlaps.regionEnd(), info.tile_width
        // );

        if (overlaps.fragEnd() <= overlaps.regionEnd()) {
            int tile_end = (int) (overlaps.fragEnd() - 1 - overlaps.regionStart()) / (int) info.tile_width;
            CountVector & buf = region_buffers[overlaps.region()][tile_end];
            if (buf.empty()) beginTileBuffer(buf);
            buf.push_back(Count{overlaps.cell()});
            // printf("tile_end %d", tile_end);
        }
        if (overlaps.fragStart() >= overlaps.regionStart()) {
            int tile_start = (int) (overlaps.fragStart() - overlaps.regionStart()) / (int) info.tile_width;
            CountVector & buf = region_buffers[overlaps.region()][tile_start];
            if (buf.empty()) beginTileBuffer(buf);
            buf.push_back(Count{overlaps.cell()});
            // printf("tile_start %d\n", tile_start);

            if (tile_start > info.completed_tiles) {
                info.completed_tiles = tile_start; // If tile_start = 1, then we have completed 1 tile (tile 0), etc.
                completed_region = overlaps.region();
                return loadTile();
            }
        }
        // } else {
        //     printf("\n");
        // }
    }
    if (overlaps.hasCompletedRegions()) {
        OverlapIterator::Region finished_region = overlaps.popCompletedRegion();
        completed_region = finished_region.index;
        region_info[completed_region].completed_tiles = region_buffers[completed_region].size();
        return loadTile();
    }
    // printf("LoadTile 3\n");
    return false;
}

} // end namespace BPCells