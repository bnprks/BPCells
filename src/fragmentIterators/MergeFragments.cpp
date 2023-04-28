#include "MergeFragments.h"

namespace BPCells {

// Merge several fragment sources into a single sorted loader.
// Cell IDs will be adjusted to be sequential between fragments
// Chromosome names & IDs must match between fragments and come in the same ordering.
//   (Use ChrSelect before merging if necessary to match ordering)

MergeFragments::MergeFragments(
    std::vector<std::unique_ptr<FragmentLoader>> &&fragments,
    const std::vector<std::string> &chr_order,
    uint32_t load_size
)
    : chr_order(chr_order)
    , start(load_size)
    , end(load_size)
    , cell(load_size) {

    if (fragments.size() < 2) throw std::runtime_error("Must have >= 2 fragments to merge");

    cell_id_offset.push_back(0);
    for (uint32_t i = 0; i < fragments.size(); i++) {
        auto f = FragmentIterator(std::move(fragments[i]));
        frags.push_back(std::move(f));

        if (frags[i].cellCount() == -1 || frags[i].chrCount() == -1 || !frags[i].isSeekable())
            throw std::runtime_error(
                "MergeFragments Error: all input fragments to merge must have known cell + chr "
                "counts and be seekable. Convert inputs to BPCells format first if needed."
            );
        cell_id_offset.push_back(cell_id_offset.back() + frags[i].cellCount());
    }

    std::unordered_map<std::string, uint32_t> chr_order_lookup;
    for (uint32_t j = 0; j < chr_order.size(); j++) {
        chr_order_lookup[chr_order[j]] = j;
    }
    source_chr.resize(fragments.size());
    for (uint32_t i = 0; i < fragments.size(); i++) {
        auto &f = frags[i];
        source_chr[i].resize(chr_order.size());
        for (auto &x : source_chr[i]) {
            x = UINT32_MAX;
        }
        int32_t chr_count = f.chrCount();
        for (int32_t j = 0; j < chr_count; j++) {
            if (chr_order_lookup.find(f.chrNames(j)) != chr_order_lookup.end()) {
                source_chr[i][chr_order_lookup[f.chrNames(j)]] = j;
            } else {
                throw std::runtime_error(
                    "MergeFragments Error: Input index " + std::to_string(i) + " has chromosome " +
                    std::string(f.chrNames(j)) + " which is not included in the output ordering."
                );
            }
        }
    }
}

bool MergeFragments::isSeekable() const { return true; }

void MergeFragments::seek(uint32_t chr_id, uint32_t base) {
    for (uint32_t i = 0; i < frags.size(); i++) {
        if (chr_id < chr_order.size()) {
            frags[i].seek(source_chr[i][chr_id], base);
        } else {
            frags[i].seek(UINT32_MAX, base);
        }
    }

    heap.clear();
    current_chr = chr_id;
}

void MergeFragments::restart() {
    for (auto &&f : frags) {
        f.restart();
    }
    heap.clear();
    current_chr = UINT32_MAX;
}

int MergeFragments::chrCount() const {
    return chr_order.size();
}

int MergeFragments::cellCount() const {
    if (frags.back().cellCount() == -1) return -1;
    return cell_id_offset.back();
}

const char *MergeFragments::chrNames(uint32_t chr_id) {
    if (chr_id >= chr_order.size()) return NULL;
    return chr_order[chr_id].c_str();
}

const char *MergeFragments::cellNames(uint32_t cell_id) {
    auto it = std::upper_bound(cell_id_offset.begin(), cell_id_offset.end(), cell_id);
    uint32_t idx = it - cell_id_offset.begin() - 1;

    if (idx == frags.size()) idx--;

    return frags[idx].cellNames(cell_id - cell_id_offset[idx]);
}

bool MergeFragments::nextChr() {
    heap.clear();

    current_chr += 1;
    if ((int64_t) current_chr >= chrCount()) {
        current_chr -= 1;
        return false;
    }

    for (uint32_t i = 0; i < frags.size(); i++) {
        if (source_chr[i][current_chr] != UINT32_MAX) {
            frags[i].seek(source_chr[i][current_chr], 0);
        } 
    }
    return true;
}

uint32_t MergeFragments::currentChr() const { return current_chr; }

bool MergeFragments::load() {
    if (heap.empty()) {
        if (current_chr >= chr_order.size()) { return false; }
        // Either initialize heap, or we are done loading
        for (uint32_t i = 0; i < frags.size(); i++) {
            if (source_chr[i][current_chr] != UINT32_MAX && frags[i].nextFrag()) {
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
