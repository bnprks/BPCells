#pragma once

#include <algorithm>
#include "../arrayIO/array_interfaces.h"
#include "FragmentIterator.h"

namespace BPCells {

// Transform a fragments loader by renaming chromosomes
// Errors may throw at load time if the number of chromosome names provided is not
// equal to the number of chromosome names in the wrapped loader
class RenameChrs : public FragmentLoaderWrapper {
private:
    const std::unique_ptr<StringReader> chr_names;
public:
    RenameChrs(FragmentLoader &loader, std::unique_ptr<StringReader> &&chr_names);

    const char* chrNames(uint32_t chr_id) const override;
};



// Transform a fragments loader by renaming cells
// Errors may throw at load time if the number of cell names provided is not
// equal to the number of cell names in the wrapped loader
class RenameCells : public FragmentLoaderWrapper {
private:
    const std::unique_ptr<StringReader> cell_names;
public:
    RenameCells(FragmentLoader &loader, std::unique_ptr<StringReader> &&cell_names);

    const char* cellNames(uint32_t cell_id) const override;
};


} // end namespace BPCells