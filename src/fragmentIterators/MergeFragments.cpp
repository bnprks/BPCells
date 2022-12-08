#include "MergeFragments.h"

namespace BPCells {

// Merge several fragment sources into a single sorted loader.
// Cell IDs will be adjusted to be sequential between fragments
// Chromosome names & IDs must match between fragments and come in the same ordering.
//   (Use ChrSelect before merging if necessary to match ordering)

MergeFragments::MergeFragments(const std::vector<FragmentLoader *> &fragments, uint32_t load_size)
    : frags_completed(fragments.size())
    , start(load_size)
    , end(load_size)
    , cell(load_size) {

    if (fragments.size() < 2) throw std::runtime_error("Must have >= 2 fragments to merge");

    cell_id_offset.push_back(0);
    for (uint32_t i = 0; i < fragments.size(); i++) {
        frags.push_back(FragmentIterator(*fragments[i]));

        if (i < fragments.size() - 1 && fragments[i]->cellCount() == -1)
            throw std::runtime_error(
                "Cannot merge fragments where cell count is not known ahead-of-time. Select cells "
                "first."
            );
        cell_id_offset.push_back(cell_id_offset.back() + fragments[i]->cellCount());

        frags[i].nextChr(); // Make sure this->nextChr doesn't get confused on the first chromosome
    }
}

bool MergeFragments::isSeekable() const {
    for (auto f : frags) {
        if (!f.isSeekable()) return false;
    }
    return true;
}

void MergeFragments::seek(uint32_t chr_id, uint32_t base) {
    for (auto f : frags)
        f.seek(chr_id, base);
    heap.clear();
}

void MergeFragments::restart() {
    for (auto f : frags) {
        f.restart();
        f.nextChr(); // Make sure this->nextChr doesn't get confused on the first chromosome
    }
    heap.clear();
    current_chr = UINT32_MAX;
    for (auto &c : frags_completed)
        c = 0;
}

int MergeFragments::chrCount() const {
    int count = frags.front().chrCount();
    for (auto f : frags) {
        if (f.chrCount() == -1) return -1;
        if (f.chrCount() != count)
            throw std::runtime_error(
                "Not all merged fragments have the same chrCount. Select identical chromosomes "
                "before merging."
            );
    }
    return count;
}

int MergeFragments::cellCount() const {
    if (frags.back().cellCount() == -1) return -1;
    return cell_id_offset.back();
}

const char *MergeFragments::chrNames(uint32_t chr_id) {
    const char *name = NULL;
    for (auto f : frags) {
        const char *f_name = f.chrNames(chr_id);
        if (f_name == NULL) continue;
        if (name == NULL) name = f_name;
        if (strcmp(name, f_name) != 0) {
            throw std::runtime_error(
                std::string("MergeFragments: Names for chromosome ID ") + std::to_string(chr_id) +
                std::string(" mismatch: ") + std::string(f_name) + std::string(" vs. ") +
                std::string(name)
            );
        }
    }
    return name;
}
const char *MergeFragments::cellNames(uint32_t cell_id) {
    auto it = std::upper_bound(cell_id_offset.begin(), cell_id_offset.end(), cell_id);
    uint32_t idx = it - cell_id_offset.begin() - 1;

    if (idx == frags.size()) idx--;

    return frags[idx].cellNames(cell_id - cell_id_offset[idx]);
}

bool MergeFragments::nextChr() {
    heap.clear();

    // Needed checks:
    //  - Each of frags[i] increases chrID in order
    //  - Each of frags[i] agrees about the current chrName
    // Curveballs:
    //   - It's okay if some of frags[i] *skip* chromosomes, and we need to
    //     handle that gracefully

    bool any_ready =
        false; // Check that we have some fragments ready to read from the current chromosome
    while (!any_ready) {
        current_chr += 1;

        bool any_remaining = false; // Check that we have some fragments that aren't completed
        for (uint32_t i = 0; i < frags.size(); i++) {
            if (frags_completed[i]) continue;
            any_remaining = true;
            if (frags[i].currentChr() > current_chr) continue;
            if (frags[i].currentChr() == current_chr) {
                any_ready = true;
                continue;
            }
            // frags[i].currentChr() < current_chr
            uint32_t prev_chr = frags[i].currentChr();
            if (!frags[i].nextChr()) {
                frags_completed[i] = 1;
                continue;
            }
            if (prev_chr >= frags[i].currentChr())
                throw std::runtime_error(
                    "Fragments to merge have out-of-order chromosomes: " + std::to_string(i)
                );
            if (frags[i].currentChr() == current_chr) any_ready = true;
        }
        if (!any_remaining) {
            current_chr -= 1;
            return false;
        }
    }

    // Check that all the chrNames agree
    const char *chr_name = frags.front().chrNames(current_chr);
    for (uint32_t i = 0; i < frags.size(); i++) {
        const char *f_name = frags[i].chrNames(current_chr);
        if (f_name == NULL || chr_name == NULL || strcmp(chr_name, f_name) != 0) {
            throw std::runtime_error(
                std::string("MergeFragments: Names for chromosome ID ") +
                std::to_string(current_chr) + std::string(" mismatch: ") + std::string(f_name) +
                std::string(" vs. ") + std::string(chr_name)
            );
        }
    }

    return true;
}

uint32_t MergeFragments::currentChr() const { return current_chr; }

bool MergeFragments::load() {
    if (heap.empty()) {
        // Either initialize heap, or we are done loading
        for (uint32_t i = 0; i < frags.size(); i++) {
            if (frags[i].currentChr() == current_chr && frags[i].nextFrag()) {
                heap.push_back(i);
            }
        }
        if (heap.empty()) return false;
        dary_heap::make_heap<4>(
            heap.begin(),
            heap.end(),
            std::bind(&MergeFragments::compare, this, std::placeholders::_1, std::placeholders::_2)
        );
    }
    uint32_t i;
    for (i = 0; i < start.size() && !heap.empty(); i++) {
        uint32_t idx = heap.front();
        start[i] = frags[idx].start();
        end[i] = frags[idx].end();
        cell[i] = frags[idx].cell() + cell_id_offset[idx];
        dary_heap::pop_heap<4>(
            heap.begin(),
            heap.end(),
            std::bind(&MergeFragments::compare, this, std::placeholders::_1, std::placeholders::_2)
        );
        if (frags[idx].nextFrag()) {
            // Just add fragment source back into the heap
            dary_heap::push_heap<4>(
                heap.begin(),
                heap.end(),
                std::bind(
                    &MergeFragments::compare, this, std::placeholders::_1, std::placeholders::_2
                )
            );
        } else {
            heap.pop_back();
        }
    }
    loaded = i;
    return loaded > 0;
}

uint32_t MergeFragments::capacity() const { return loaded; }

uint32_t *MergeFragments::cellData() { return cell.data(); }
uint32_t *MergeFragments::startData() { return start.data(); }
uint32_t *MergeFragments::endData() { return end.data(); }

} // end namespace BPCells
