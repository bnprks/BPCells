#pragma once
#include <unordered_map>
#include "FragmentsIterator.h"

namespace BPCells {

// Transform a fragments iterator by renaming or filtering cell IDs 
class CellIndexSelect : public FragmentsLoaderWrapper {
private:
    const std::vector<uint32_t> cell_indices;
    // Reverse lookup for cell indices -- reverse_indices[i] gives the output cell_id
    // for input cell_id i
    std::vector<uint32_t> reverse_indices;
public:
    // cell_indices -- vector with length <= the number of chromosomes in the input
    //     FragmentsIterator. The output cell `i` will come from input cell
    //     `chr_assignments[i]`. The entries of cell_indices must be unique
    CellIndexSelect(FragmentsLoader &loader, const std::vector<uint32_t> cell_indices);

    ~CellIndexSelect() = default;

    // Return the number of cells/chromosomes, or return -1 if this number is 
    // not known ahead of time
    int cellCount() const override;

    const char* cellNames(uint32_t cell_id) const override;

    // Return number of items loaded. Should repeatedly return 0 at the end of a chromosome.
    // Return -1 for error
    int32_t load(uint32_t count, FragmentArray &buffer) override;
};


// Transform a fragments iterator by renaming or filtering cell IDs 
class CellNameSelect : public FragmentsLoaderWrapper {
private:
    const std::vector<std::string> cell_names;
    std::unordered_map<std::string, uint32_t> output_index;
    // Reverse lookup for cell indices -- reverse_indices[i] gives the output cell_id
    // for input cell_id i
    std::vector<uint32_t> reverse_indices;

    // Return output cell ID given an input cell id, or UINT32_MAX if the cell
    // shouldn't be output
    uint32_t getOutputCellID(uint32_t input_cell_id); 
public:
    // cell_names -- vector with length <= the number of chromosomes in the input
    //     FragmentsLoader. The output cell `i` will come from input cell
    //     `chr_assignments[i]`. The entries of cell_names must be unique
    CellNameSelect(FragmentsLoader &loader, const std::vector<std::string> cell_names);

    ~CellNameSelect() = default;
    
    int cellCount() const override;

    const char* cellNames(uint32_t cell_id) const override;
    
    // Return number of items loaded. Should repeatedly return 0 at the end of a chromosome.
    // Return -1 for error
    int32_t load(uint32_t count, FragmentArray &buffer) override;
};


} // end namespace BPCells