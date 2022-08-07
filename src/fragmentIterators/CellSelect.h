#pragma once
#include <unordered_map>
#include "FragmentIterator.h"

namespace BPCells {

// Transform a fragments iterator by renaming or filtering cell IDs 
class CellIndexSelect : public FragmentLoaderWrapper {
private:
    uint32_t loaded = 0;
    const std::vector<uint32_t> cell_indices;
    // Reverse lookup for cell indices -- reverse_indices[i] gives the output cell_id
    // for input cell_id i
    std::vector<uint32_t> reverse_indices;
public:
    // cell_indices -- vector with length <= the number of chromosomes in the input
    //     FragmentIterator. The output cell `i` will come from input cell
    //     `chr_assignments[i]`. The entries of cell_indices must be unique
    CellIndexSelect(FragmentLoader &loader, const std::vector<uint32_t> cell_indices);

    ~CellIndexSelect() = default;

    // Return the number of cells/chromosomes, or return -1 if this number is 
    // not known ahead of time
    int cellCount() const override;

    const char* cellNames(uint32_t cell_id) override;

    bool load() override;
    uint32_t capacity() const override;


};


// Transform a fragments iterator by renaming or filtering cell IDs 
class CellNameSelect : public FragmentLoaderWrapper {
private:
    uint32_t loaded;
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
    //     FragmentLoader. The output cell `i` will come from input cell
    //     `chr_assignments[i]`. The entries of cell_names must be unique
    CellNameSelect(FragmentLoader &loader, const std::vector<std::string> cell_names);

    ~CellNameSelect() = default;
    
    int cellCount() const override;

    const char* cellNames(uint32_t cell_id) override;

    bool load() override;
    uint32_t capacity() const override;
};


} // end namespace BPCells