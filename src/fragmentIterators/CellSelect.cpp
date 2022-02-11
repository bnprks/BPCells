#include "CellSelect.h"

namespace BPCells {

// cell_indices -- vector with length <= the number of chromosomes in the input
//     FragmentsLoader. The output cell `i` will come from input cell
//     `chr_assignments[i]`. The entries of cell_indices must be unique
CellIndexSelect::CellIndexSelect(FragmentsLoader &loader, const std::vector<uint32_t> cell_indices) : 
    FragmentsLoaderWrapper(loader),
    cell_indices(cell_indices) {
    for (int i = 0; i < cell_indices.size(); i++) {
        if (reverse_indices.size() <= cell_indices[i]) 
            reverse_indices.resize(cell_indices[i] + 1, UINT32_MAX);
        if (reverse_indices[cell_indices[i]] < UINT32_MAX)
            throw std::invalid_argument("CellSelect maps same input cell to two output IDs");
        reverse_indices[cell_indices[i]] = i;
    }
};

int CellIndexSelect::cellCount() const { return cell_indices.size(); };

const char* CellIndexSelect::cellNames(uint32_t cell_id) const {
    if (cell_id >= cell_indices.size()) return NULL;
    return loader.cellNames(cell_indices[cell_id]); 
};

// Return number of items loaded. Should repeatedly return 0 at the end of a chromosome.
// Return -1 for error
int32_t CellIndexSelect::load(uint32_t count, FragmentArray buffer) {
    // load and filter until we load without filtering out everything
    while (true) {
        int32_t ret = loader.load(count, buffer);
        if (ret <= 0) return ret;
        uint32_t in, out;
        for (in = 0, out=0; in < ret; in++) {
            if (reverse_indices[buffer.cell[in]] != UINT32_MAX) {
                buffer.cell[out] = reverse_indices[buffer.cell[in]];
                buffer.start[out] = buffer.start[in];
                buffer.end[out] = buffer.end[in];
                out++;
            }
        }
        if (out > 0) return out;
    }
};


uint32_t CellNameSelect::getOutputCellID(uint32_t input_cell_id) {
    // Update the reverse_indices map up to input_cell_id
    if (input_cell_id >= reverse_indices.size()) {
        auto old_size = reverse_indices.size();
        reverse_indices.resize(input_cell_id+1, UINT32_MAX);
        for (auto i = old_size; i < reverse_indices.size(); i++) {
            auto res = output_index.find(loader.cellNames(i));
            if (res != output_index.end()) {
                reverse_indices[i] = res->second;
            }
        }
    }
    return reverse_indices[input_cell_id];
}

// cell_names -- vector with length <= the number of chromosomes in the input
//     FragmentsLoader. The output cell `i` will come from input cell
//     `chr_assignments[i]`. The entries of cell_names must be unique
CellNameSelect::CellNameSelect(FragmentsLoader &loader, const std::vector<std::string> cell_names) :
    FragmentsLoaderWrapper(loader),
    cell_names(cell_names) {
    for (int i = 0; i < cell_names.size(); i++) {    
        if(output_index.find(cell_names[i]) != output_index.end())
            throw std::invalid_argument("CellSelect maps same input cell to two output IDs");
        output_index[cell_names[i]] = i; 
    }
}

int CellNameSelect::cellCount() const { return cell_names.size(); };

const char* CellNameSelect::cellNames(uint32_t cell_id) const {
    if (cell_id >= cell_names.size()) return NULL;
    return cell_names[cell_id].c_str(); 
};

int32_t CellNameSelect::load(uint32_t count, FragmentArray buffer) {
    // load and filter until we load without filtering out everything
    while (true) {
        int32_t ret = loader.load(count, buffer);
        if (ret <= 0) return ret;
        uint32_t in, out;
        for (in = 0, out=0; in < ret; in++) {
            uint32_t new_cell_id = getOutputCellID(buffer.cell[in]);
            if (new_cell_id != UINT32_MAX) {
                buffer.cell[out] = new_cell_id;
                buffer.start[out] = buffer.start[in];
                buffer.end[out] = buffer.end[in];
                out++;
            }
        }
        if (out > 0) {
            return out;
        }
    }
};

} // end namespace BPCells