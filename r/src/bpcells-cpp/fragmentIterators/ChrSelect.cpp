// Copyright 2021 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "ChrSelect.h"

namespace BPCells {

// chr_assignments -- vector with length <= the number of chromosomes in the input
//     FragmentsIterator. The output chromosome `i` will come from input chromosome
//     `chr_assignments[i]`. The entries of chr_assignments must be unique
ChrIndexSelect::ChrIndexSelect(
    std::unique_ptr<FragmentLoader> &&loader, const std::vector<uint32_t> chr_assignments
)
    : FragmentLoaderWrapper(std::move(loader))
    , chr_assignments(chr_assignments) {

    std::vector<uint32_t> seen_id;
    for (auto new_id : chr_assignments) {
        if (seen_id.size() <= new_id) seen_id.resize(new_id + 1, 0);
        if (seen_id[new_id])
            throw std::invalid_argument("ChrSelect maps same input chromosome to two output IDs");
        seen_id[new_id]++;
    }
}

void ChrIndexSelect::restart() {
    current_chr = UINT32_MAX;
    loader->restart();
}

int ChrIndexSelect::chrCount() const { return chr_assignments.size(); }

const char *ChrIndexSelect::chrNames(uint32_t chr_id) {
    if (chr_id >= chr_assignments.size()) return NULL;
    return loader->chrNames(chr_assignments[chr_id]);
}

bool ChrIndexSelect::nextChr() {
    // If input is seekable, just load chromosomes in output order
    if (loader->isSeekable()) {
        current_chr += 1;
        if ((int64_t)current_chr >= chrCount()) {
            current_chr -= 1;
            return false;
        }
        loader->seek(chr_assignments[current_chr], 0);
        return true;
    }
    // Otherwise, just load from input order and say the right chromosome IDs going out
    bool res = loader->nextChr();
    while (res) {
        uint32_t chr_id = loader->currentChr();
        auto result = std::find(chr_assignments.begin(), chr_assignments.end(), chr_id);
        if (result != chr_assignments.end()) break;
        res = loader->nextChr();
    }
    return res;
}

uint32_t ChrIndexSelect::currentChr() const {
    auto result = std::find(chr_assignments.begin(), chr_assignments.end(), loader->currentChr());
    if (result == chr_assignments.end())
        throw std::invalid_argument("ChrSelect does not have a chromosome assigned to requested ID"
        );

    return result - chr_assignments.begin();
}

void ChrIndexSelect::seek(uint32_t chr_id, uint32_t base) {
    current_chr = chr_id;
    if (chr_id < chr_assignments.size()) {
        loader->seek(chr_assignments[chr_id], base);
    } else {
        loader->seek(UINT32_MAX, base);
    }
}

ChrNameSelect::ChrNameSelect(
    std::unique_ptr<FragmentLoader> &&loader, const std::vector<std::string> chr_names
)
    : FragmentLoaderWrapper(std::move(loader))
    , chr_names(chr_names) {
    for (uint32_t i = 0; i < chr_names.size(); i++) {
        if (output_index.find(chr_names[i]) != output_index.end())
            throw std::invalid_argument("ChrSelect maps same input chromosome to two output IDs");
        output_index[chr_names[i]] = i;
    }

    if (this->loader->isSeekable()) {
        input_index.resize(chr_names.size());
        // An important assumption here is that if we seek to UINT32_MAX then the loader will just
        // return false on next load() call
        for (auto &i : input_index) {
            i = UINT32_MAX;
        }
        int32_t chr_count = this->loader->chrCount();
        for (int i = 0; i < chr_count; i++) {
            if (output_index.find(this->loader->chrNames(i)) != output_index.end()) {
                input_index[output_index[this->loader->chrNames(i)]] = i;
            }
        }
    }
}

void ChrNameSelect::restart() {
    current_chr = UINT32_MAX;
    loader->restart();
}

int ChrNameSelect::chrCount() const { return chr_names.size(); }

const char *ChrNameSelect::chrNames(uint32_t chr_id) {
    if (chr_id >= chr_names.size()) return NULL;
    return chr_names[chr_id].c_str();
}

bool ChrNameSelect::nextChr() {
    // If input is seekable, just load chromosomes in output order
    if (loader->isSeekable()) {
        current_chr += 1;
        if ((int64_t)current_chr >= chrCount()) {
            current_chr -= 1;
            return false;
        }
        loader->seek(input_index[current_chr], 0);
        return true;
    }
    // Otherwise, just load from input order and say the right chromosome IDs going out
    bool res = loader->nextChr();
    while (res) {
        uint32_t chr_id = loader->currentChr();
        auto result = output_index.find(loader->chrNames(chr_id));
        if (result != output_index.end()) break;
        res = loader->nextChr();
    }
    return res;
}

uint32_t ChrNameSelect::currentChr() const {
    if (loader->isSeekable()) {
        return current_chr;
    }
    auto res = output_index.at(loader->chrNames(loader->currentChr()));
    return res;
}

void ChrNameSelect::seek(uint32_t chr_id, uint32_t base) {
    current_chr = chr_id;
    if (chr_id < input_index.size()) {
        loader->seek(input_index[chr_id], base);
    } else {
        loader->seek(UINT32_MAX, base);
    }
}

} // end namespace BPCells
