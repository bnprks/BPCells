#include "Rename.h"

namespace BPCells {

RenameChrs::RenameChrs(FragmentLoader &loader, std::unique_ptr<StringReader> &&chr_names)
    : FragmentLoaderWrapper(loader)
    , chr_names(std::move(chr_names)) {}

const char *RenameChrs::chrNames(uint32_t chr_id) {
    const char *ret = chr_names->get(chr_id);
    if ((ret == NULL) != (loader.chrNames(chr_id) == NULL)) {
        throw std::runtime_error(
            std::string("RenameChrs: mismatched number of named chromosomes, chr_id=") +
            std::to_string(chr_id)
        );
    }
    return ret;
}

RenameCells::RenameCells(FragmentLoader &loader, std::unique_ptr<StringReader> &&cell_names)
    : FragmentLoaderWrapper(loader)
    , cell_names(std::move(cell_names)) {}

const char *RenameCells::cellNames(uint32_t cell_id) {
    const char *ret = cell_names->get(cell_id);
    if ((ret == NULL) != (loader.cellNames(cell_id) == NULL)) {
        throw std::runtime_error(
            std::string("RenameCells: mismatched number of named cells, cell_id=") +
            std::to_string(cell_id)
        );
    }
    return ret;
}

// Tranform a fragments loader by adding a prefix to all the cell names
// This is a limited alternative to RenameCells for situations where
// the cell names are not known ahead-of-time

PrefixCells::PrefixCells(FragmentLoader &loader, std::string prefix)
    : FragmentLoaderWrapper(loader)
    , prefix(prefix) {

    name_buffer.resize(prefix.size() + 1);
    strncpy(name_buffer.data(), prefix.c_str(), prefix.size());
}

const char *PrefixCells::cellNames(uint32_t cell_id) {
    const char *raw = this->loader.cellNames(cell_id);
    if (raw == NULL) return NULL;
    if (strlen(raw) + 1 + prefix.size() >= name_buffer.size())
        name_buffer.resize(strlen(raw) + 1 + prefix.size());
    strcpy(name_buffer.data() + prefix.size(), raw);
    return name_buffer.data();
}

} // end namespace BPCells
