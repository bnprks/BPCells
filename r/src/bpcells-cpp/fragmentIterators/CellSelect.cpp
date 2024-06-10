// Copyright 2021 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "CellSelect.h"

namespace BPCells {

// cell_indices -- vector with length <= the number of chromosomes in the input
//     FragmentLoader. The output cell `i` will come from input cell
//     `cell_indices[i]`. The entries of cell_indices must be unique
CellIndexSelect::CellIndexSelect(std::unique_ptr<FragmentLoader> &&loader, const std::vector<uint32_t> cell_indices)
    : FragmentLoaderWrapper(std::move(loader))
    , cell_indices(cell_indices) {
    for (uint32_t i = 0; i < cell_indices.size(); i++) {
        if (reverse_indices.size() <= cell_indices[i])
            reverse_indices.resize(cell_indices[i] + 1, UINT32_MAX);
        if (reverse_indices[cell_indices[i]] < UINT32_MAX)
            throw std::invalid_argument("CellSelect maps same input cell to two output IDs");
        reverse_indices[cell_indices[i]] = i;
    }
}

int CellIndexSelect::cellCount() const { return cell_indices.size(); }

const char *CellIndexSelect::cellNames(uint32_t cell_id) {
    if (cell_id >= cell_indices.size()) return NULL;
    return loader->cellNames(cell_indices[cell_id]);
}

bool CellIndexSelect::load() {
    loaded = 0;
    // load and filter until we load without filtering out everything
    while (loaded == 0) {
        if (!loader->load()) return false;

        uint32_t *cell = loader->cellData();
        uint32_t *start = loader->startData();
        uint32_t *end = loader->endData();
        uint32_t capacity = loader->capacity();
        for (uint32_t i = 0; i < capacity; i++) {
            cell[loaded] = cell[i] < reverse_indices.size() ? reverse_indices[cell[i]] : UINT32_MAX;
            start[loaded] = start[i];
            end[loaded] = end[i];
            loaded += cell[loaded] != UINT32_MAX;
        }
    }
    return true;
}

uint32_t CellIndexSelect::capacity() const { return loaded; }

uint32_t CellNameSelect::getOutputCellID(uint32_t input_cell_id) {
    // Update the reverse_indices map up to input_cell_id
    if (input_cell_id >= reverse_indices.size()) {
        auto old_size = reverse_indices.size();
        reverse_indices.resize(input_cell_id + 1, UINT32_MAX);
        for (auto i = old_size; i < reverse_indices.size(); i++) {
            auto res = output_index.find(loader->cellNames(i));
            if (res != output_index.end()) {
                reverse_indices[i] = res->second;
            }
        }
    }
    return reverse_indices[input_cell_id];
}

// cell_names -- vector with length <= the number of chromosomes in the input
//     FragmentLoader. The output cell `i` will come from input cell
//     `cell_names[i]`. The entries of cell_names must be unique
CellNameSelect::CellNameSelect(std::unique_ptr<FragmentLoader> &&loader, const std::vector<std::string> cell_names)
    : FragmentLoaderWrapper(std::move(loader))
    , cell_names(cell_names) {
    for (uint32_t i = 0; i < cell_names.size(); i++) {
        if (output_index.find(cell_names[i]) != output_index.end())
            throw std::invalid_argument("CellSelect maps same input cell to two output IDs");
        output_index[cell_names[i]] = i;
    }
}

int CellNameSelect::cellCount() const { return cell_names.size(); }

const char *CellNameSelect::cellNames(uint32_t cell_id) {
    if (cell_id >= cell_names.size()) return NULL;
    return cell_names[cell_id].c_str();
}

bool CellNameSelect::load() {
    loaded = 0;
    // load and filter until we load without filtering out everything
    while (loaded == 0) {
        if (!loader->load()) return false;

        uint32_t *cell = loader->cellData();
        uint32_t *start = loader->startData();
        uint32_t *end = loader->endData();
        uint32_t capacity = loader->capacity();
        for (uint32_t i = 0; i < capacity; i++) {
            uint32_t new_cell_id = getOutputCellID(cell[i]);
            cell[loaded] = new_cell_id;
            start[loaded] = start[i];
            end[loaded] = end[i];
            loaded += new_cell_id != UINT32_MAX;
        }
    }
    return true;
}

uint32_t CellNameSelect::capacity() const { return loaded; }

CellMerge::CellMerge(
    std::unique_ptr<FragmentLoader> &&loader,
    const std::vector<uint32_t> group_ids,
    std::unique_ptr<StringReader> &&group_names
)
    : FragmentLoaderWrapper(std::move(loader))
    , group_ids(group_ids)
    , group_names(std::move(group_names)) {

    group_count = 0;
    for (auto id : group_ids) {
        group_count = std::max(group_count, id + 1);
    }

    if (group_count > this->group_names->size()) {
        throw std::invalid_argument("CellMerge has more groups given than group_names");
    }
    if ((int64_t) group_ids.size() != this->loader->cellCount()) {
        throw std::invalid_argument("CellMerge number of input cells != length of group_ids");
    }
    group_count = this->group_names->size();
}

int CellMerge::cellCount() const { return group_count; }

const char *CellMerge::cellNames(uint32_t cell_id) { return group_names->get(cell_id); }

bool CellMerge::load() {
    if (!loader->load()) return false;
    uint32_t capacity = loader->capacity();
    uint32_t *cell = loader->cellData();
    for (uint32_t i = 0; i < capacity; i++) {
        cell[i] = group_ids[cell[i]];
    }
    return true;
}

} // end namespace BPCells
