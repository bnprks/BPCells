#include "Rename.h"

namespace BPCells {

RenameChrs::RenameChrs(FragmentLoader &loader, std::unique_ptr<StringReader> &&chr_names) :
    FragmentLoaderWrapper(loader), chr_names(std::move(chr_names)) {}

const char* RenameChrs::chrNames(uint32_t chr_id) const {
    const char* ret = chr_names->get(chr_id);
    if ((ret == NULL) != (loader.chrNames(chr_id) == NULL)) {
        throw std::runtime_error(std::string("RenameChrs: mismatched number of named chromosomes, chr_id=") + std::to_string(chr_id));
    }
    return ret;
}

RenameCells::RenameCells(FragmentLoader &loader, std::unique_ptr<StringReader> &&cell_names) :
    FragmentLoaderWrapper(loader), cell_names(std::move(cell_names)) {}

const char* RenameCells::cellNames(uint32_t cell_id) const {
    const char* ret = cell_names->get(cell_id);
    if ((ret == NULL) != (loader.cellNames(cell_id) == NULL)) {
        throw std::runtime_error(std::string("RenameCells: mismatched number of named cells, cell_id=") + std::to_string(cell_id));
    }
    return ret;
}


} // end namespace BPCells