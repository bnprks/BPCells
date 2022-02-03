#pragma once

namespace BPCells {

template<class T>
class UnpackedFrags {
public:
    // end_max[i] = max(end[0:128*(i+1)]). This enables seeking at the cost of loading
    // 1% of the data size in to memory up-front
    T start, end, cell, end_max; 
    UnpackedFrags() = default;
    UnpackedFrags(UnpackedFrags&&) = default;
    UnpackedFrags(const UnpackedFrags&) = default;
    UnpackedFrags& operator=(UnpackedFrags&&) = default;
    UnpackedFrags& operator=(const UnpackedFrags&) = delete;
};

template<class T>
class PackedFrags {
public:
    // Original design considered using FOR packing for cell_ids.
    // Since we are currently just storing pos-sorted without
    // much cell chunking I'll disable the FOR encoding
    T start_data, start_idx, start_starts,
        end_data, end_idx, end_max,
        cell_data, cell_idx, count;
    PackedFrags() = default;
    PackedFrags(PackedFrags&&) = default;
    PackedFrags(const PackedFrags&) = default;
    PackedFrags& operator=(PackedFrags&&) = default;
    PackedFrags& operator=(const PackedFrags&) = delete;
};

} // end namespace BPCells