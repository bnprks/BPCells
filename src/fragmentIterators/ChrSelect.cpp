#include "ChrSelect.h"

namespace BPCells {

// chr_assignments -- vector with length <= the number of chromosomes in the input
//     FragmentsIterator. The output chromosome `i` will come from input chromosome
//     `chr_assignments[i]`. The entries of chr_assignments must be unique
ChrIndexSelect::ChrIndexSelect(FragmentLoader &loader, const std::vector<uint32_t> chr_assignments)
    : FragmentLoaderWrapper(loader)
    , chr_assignments(chr_assignments) {

    std::vector<uint32_t> seen_id;
    for (auto new_id : chr_assignments) {
        if (seen_id.size() <= new_id) seen_id.resize(new_id + 1, 0);
        if (seen_id[new_id])
            throw std::invalid_argument("ChrSelect maps same input chromosome to two output IDs");
        seen_id[new_id]++;
    }
}

int ChrIndexSelect::chrCount() const { return chr_assignments.size(); }

const char *ChrIndexSelect::chrNames(uint32_t chr_id) {
    if (chr_id >= chr_assignments.size()) return NULL;
    return loader.chrNames(chr_assignments[chr_id]);
}

bool ChrIndexSelect::nextChr() {
    bool res = loader.nextChr();
    while (res) {
        uint32_t chr_id = loader.currentChr();
        auto result = std::find(chr_assignments.begin(), chr_assignments.end(), chr_id);
        if (result != chr_assignments.end()) break;
        res = loader.nextChr();
    }
    return res;
}

uint32_t ChrIndexSelect::currentChr() const {
    auto result = std::find(chr_assignments.begin(), chr_assignments.end(), loader.currentChr());
    if (result == chr_assignments.end())
        throw std::invalid_argument("ChrSelect does not have a chromosome assigned to requested ID"
        );

    return result - chr_assignments.begin();
}

void ChrIndexSelect::seek(uint32_t chr_id, uint32_t base) {
    loader.seek(chr_assignments[chr_id], base);
}

ChrNameSelect::ChrNameSelect(FragmentLoader &loader, const std::vector<std::string> chr_names)
    : FragmentLoaderWrapper(loader)
    , chr_names(chr_names) {
    for (uint32_t i = 0; i < chr_names.size(); i++) {
        if (output_index.find(chr_names[i]) != output_index.end())
            throw std::invalid_argument("ChrSelect maps same input chromosome to two output IDs");
        output_index[chr_names[i]] = i;
    }

    if (loader.isSeekable()) {
        input_index.resize(chr_names.size());
        int32_t chr_count = loader.chrCount();
        for (int i = 0; i < chr_count; i++) {
            if (output_index.find(loader.chrNames(i)) != output_index.end()) {
                input_index[output_index[loader.chrNames(i)]] = i;
            }
        }
    }
}

int ChrNameSelect::chrCount() const { return chr_names.size(); }

const char *ChrNameSelect::chrNames(uint32_t chr_id) {
    if (chr_id >= chr_names.size()) return NULL;
    return chr_names[chr_id].c_str();
}

bool ChrNameSelect::nextChr() {
    bool res = loader.nextChr();
    ;
    while (res) {
        uint32_t chr_id = loader.currentChr();
        auto result = output_index.find(loader.chrNames(chr_id));
        if (result != output_index.end()) break;
        res = loader.nextChr();
    }
    return res;
}

uint32_t ChrNameSelect::currentChr() const {
    auto res = output_index.at(loader.chrNames(loader.currentChr()));
    return res;
}

void ChrNameSelect::seek(uint32_t chr_id, uint32_t base) { loader.seek(input_index[chr_id], base); }

} // end namespace BPCells