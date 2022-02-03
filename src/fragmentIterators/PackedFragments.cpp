#include "PackedFragments.h"

namespace BPCells {

PackedFragments::PackedFragments(const std::vector<packed_frags> frags, 
        const std::vector<std::string> cell_names, const std::vector<std::string> chr_names) : 
        frags(frags), chr_names(chr_names), cell_names(cell_names),
        current_chr(-1), current_block(0)
        {}

bool PackedFragments::isSeekable() const {return true;};

void PackedFragments::seek(uint32_t chr_id, uint32_t base) {
    current_chr = chr_id;

    current_block = std::lower_bound(
            &frags[current_chr].end_max[0], 
            &frags[current_chr].end_max[frags[current_chr].end_max.size],
            base) - &frags[current_chr].end_max[0];
}

void PackedFragments::restart() {
    current_chr = 0;
    current_block = 0;
};

int PackedFragments::chrCount() const { return chr_names.size(); };
int PackedFragments::cellCount() const { return cell_names.size(); };

const char* PackedFragments::chrNames(uint32_t chr_id) const {
    if (chr_id >= chr_names.size()) return NULL;
    return chr_names[chr_id].c_str(); 
};
const char* PackedFragments::cellNames(uint32_t cell_id) const {
    if (cell_id >= cell_names.size()) return NULL;
    return cell_names[cell_id].c_str(); 
};

bool PackedFragments::nextChr() {
    current_chr++;
    current_block = 0;
    if (current_chr >= frags.size()) return false;
    return true;
};

uint32_t PackedFragments::currentChr() const {return current_chr;}

int32_t PackedFragments::load(uint32_t count, FragmentArray &buf) {
    if (current_block >= frags[current_chr].n_chunks()) return 0;
    if (count < 128) return -1;
    count = count & ~(128-1); // Round count to nearest 128
    
    for (int i = 0; i < count; i += 128) {
        unpack_128_frags(
            frags[current_chr], 
            &buf.start[i], 
            &buf.end[i], 
            &buf.cell[i], 
            current_block
        );
        current_block += 1;
        if (current_block == frags[current_chr].n_chunks()) {
            return i + 128 - (128*frags[current_chr].n_chunks() - frags[current_chr].n_frags());
        }
    }
    return count;
}


std::vector<packed_frags> PackedFragmentsWriter::getPackedFrags() {return packed_frags_vec; }
std::vector<std::string> PackedFragmentsWriter::getChrNames() {return chr_names; }
std::vector<std::string> PackedFragmentsWriter::getCellNames() {return cell_names; }

vec_uint32_t PackedFragmentsWriter::convert_to_vec(std::vector<uint32_t> &v) {
    vec_uint32_t out;
    out.data = v.data();
    out.capacity = v.capacity();
    out.size = v.size();
    return out;
}

bool PackedFragmentsWriter::write(FragmentsIterator &fragments, void (*checkInterrupt)(void)) {
    uint32_t start_bits, end_bits, cell_id_bits;
    uint32_t start_idx, end_idx, cell_id_idx;
    uint32_t start_buffer[128];
    uint32_t end_buffer[128];
    uint32_t cell_buffer[128];
    uint32_t chr_id, i;
    bool chr_done;
    
    while (fragments.nextChr()) {
        chr_id = fragments.currentChr();

        vec_packed_frags out;

        start_idx = 0;
        end_idx = 0;
        cell_id_idx = 0;
        out.start_idx.push_back(start_idx);
        out.end_idx.push_back(end_idx);
        out.cell_idx.push_back(cell_id_idx);

        out.len = 0;

        chr_done = false;
        while (!chr_done) {
            for (i = 0; i < 128; i++) {
                if (!fragments.nextFrag()) {
                    chr_done = true;
                    break;
                };
                start_buffer[i] = fragments.start();
                end_buffer[i] = fragments.end();
                cell_buffer[i] = fragments.cell();
            }
            if (i == 0) break; //Handle chr with multiple of 128 frags
            
            out.len += i;    

            if (checkInterrupt != NULL && (out.len / 128) % 1024 == 0) checkInterrupt();

            // Cleanup loop in case we're at end of chromosome and
            // couldn't read all 128
            for (;i < 128; i++) {
                start_buffer[i] = start_buffer[i-1];
                end_buffer[i] = end_buffer[i-1];
                cell_buffer[i] = cell_buffer[i-1];
            }
            // Pack starts
            start_bits = simdmaxbitsd1(start_buffer[0], start_buffer);
            out.start_start.push_back(start_buffer[0]);
            out.start_data.resize(start_idx + start_bits*4);
            out.start_idx.push_back(start_idx + start_bits*4);
            simdpackd1(start_buffer[0], start_buffer, &out.start_data.data()[start_idx], start_bits);
            start_idx += start_bits*4;

            // Pack ends
            out.end_max.push_back(std::max(
                simdmax(&end_buffer[0]),
                out.end_max.size() == 0 ? 0 : out.end_max[out.end_max.size() - 1]
            ));
            simdsubtract(end_buffer, start_buffer);
            end_bits = simdmaxbits(end_buffer);
            out.end_data.resize(end_idx + end_bits*4);
            out.end_idx.push_back(end_idx + end_bits*4);
            simdpack(end_buffer, &out.end_data.data()[end_idx], end_bits);
            end_idx += end_bits*4;
            

            // Pack cell_ids
            cell_id_bits = simdmaxbits(cell_buffer);
        
            out.cell_data.resize(cell_id_idx + cell_id_bits*4);
            out.cell_idx.push_back(cell_id_idx + cell_id_bits*4);

            simdpack(cell_buffer, &out.cell_data.data()[cell_id_idx], cell_id_bits);
            cell_id_idx += cell_id_bits*4;
        } // End loop over fragments in chromosome

        if (chr_id >= frags.size()) {
            frags.resize(chr_id+1);
        }

        frags[chr_id] = out;
    } // End loop over chromosomes

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

    // Copy references needed to packed_frags
    for (int i = 0; i < frags.size(); i++) {
        packed_frags f;
        f.starts.data = convert_to_vec(frags[i].start_data);
        f.starts.starts = convert_to_vec(frags[i].start_start);
        f.starts.idx = convert_to_vec(frags[i].start_idx);
        f.starts.len = frags[i].len;

        f.ends.data = convert_to_vec(frags[i].end_data);
        f.ends.idx = convert_to_vec(frags[i].end_idx);
        f.end_max = convert_to_vec(frags[i].end_max);
        f.ends.len = frags[i].len;

        f.cell_ids.data = convert_to_vec(frags[i].cell_data);
        f.cell_ids.idx = convert_to_vec(frags[i].cell_idx);
        f.cell_ids.len = frags[i].len;
        packed_frags_vec.push_back(f);
    }

    return true;
}

} // end namespace BPCells