#include "UnpackedFragments.h"

namespace BPCells {

UnpackedFragments::UnpackedFragments(
    const std::vector<vec_uint32_t> start, const std::vector<vec_uint32_t> end,
    const std::vector<vec_uint32_t> cell_id,  
    const std::vector<std::string> cell_names, const std::vector<std::string> chr_names) : 
    start(start), end(end), cell_id(cell_id),
    chr_names(chr_names), cell_names(cell_names),
    current_chr(-1), current_fragment(0)  {

    if (start.size() != end.size() || start.size() != cell_id.size() || 
        start.size() != chr_names.size()) {
        throw std::invalid_argument("Start, end, cell_id, and chr_names must all be same length");
    }
    for (size_t i = 0; i < start.size(); i++) {
        if (start[i].size != end[i].size || start[i].size != cell_id[i].size) {
            throw std::invalid_argument("Start, end, cell_id must have matching lengths in each element");
        }
    }
}



bool UnpackedFragments::isSeekable() const {return false;};

void UnpackedFragments::seek(uint32_t chr_id, uint32_t base) {
    throw std::invalid_argument("Cannot seek UnpackedFragments");
}

void UnpackedFragments::restart() { current_chr = 0; current_fragment=0; };


// Return the number of cells/chromosomes, or return -1 if this number is 
// not known ahead of time
int UnpackedFragments::chrCount() const { return chr_names.size(); };
int UnpackedFragments::cellCount() const { return cell_names.size(); };

// Return name for a given chr_id or cell_id. Only valid to call
// for chromosme or cell_ids that have been actually returned by the iterator
const char* UnpackedFragments::chrNames(uint32_t chr_id) const {
    if (chr_id >= chr_names.size()) return NULL;
    return chr_names[chr_id].c_str(); 
};
const char* UnpackedFragments::cellNames(uint32_t cell_id) const {
    if (cell_id >= cell_names.size()) return NULL;
    return cell_names[cell_id].c_str(); 
};

// Advance the iterator to the next chromosome. Return false if there are no more chromosomes
bool UnpackedFragments::nextChr() {
    current_chr++;
    if (current_chr >= start.size()) {
        current_chr -= 1;
        return false;
    }
    current_fragment = 0;
    return true;
};
// Return chromosome ID of current fragment
uint32_t UnpackedFragments::currentChr() const {return current_chr;}


int32_t UnpackedFragments::load(uint32_t count, FragmentArray &buf) {
    size_t output = 0;
    while(output < count && current_fragment < start[current_chr].size) {
        buf.start[output] = start[current_chr][current_fragment];
        buf.end[output] = end[current_chr][current_fragment];
        buf.cell[output] = cell_id[current_chr][current_fragment];
        
        current_fragment += 1;
        output += 1;
    }
    
    return output;
}


const std::vector<UnpackedFragmentsWriter::unpacked_frags> *UnpackedFragmentsWriter::getFrags() {return &frags; }
std::vector<std::string> UnpackedFragmentsWriter::getChrNames() {return chr_names; }
std::vector<std::string> UnpackedFragmentsWriter::getCellNames() {return cell_names; }

bool UnpackedFragmentsWriter::write(FragmentsIterator &fragments, void (*checkInterrupt)(void)) {
    uint32_t chr_id;
    uint32_t max_cell_id = 0;

    while (fragments.nextChr()) {
        chr_id = fragments.currentChr();
        struct unpacked_frags out;

        while (fragments.nextFrag()) {
            out.len++;
            out.start_data.push_back(fragments.start());
            out.end_data.push_back(fragments.end());
            out.cell_data.push_back(fragments.cell());
            max_cell_id = std::max(max_cell_id, fragments.cell());

            if (checkInterrupt != NULL && (out.len % 8192) == 0) checkInterrupt();
        } 
        
        if (chr_id >= frags.size()) {
            frags.resize(chr_id+1);
        }
        frags[chr_id] = out;
    } 

    // Get cell and chromosome names
    for (int i = 0; ; i++) {
        const char* cell_name = fragments.cellNames(i);
        if (cell_name == NULL) break;
        cell_names.push_back(cell_name);
    }
    for (int i = 0; ; i++) {
        const char* chr_name = fragments.chrNames(i);
        if (chr_name == NULL) break;
        chr_names.push_back(chr_name);
    }

    return true;
}



} // end namespace BPCells