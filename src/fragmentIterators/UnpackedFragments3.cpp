#include "UnpackedFragments3.h"

namespace BPCells {

UnpackedFragments3::UnpackedFragments3(ReaderPtr &&cell, ReaderPtr &&start, ReaderPtr &&end, 
        ReaderPtr &&end_max, ReaderPtr && chr_ptr, 
        std::unique_ptr<StringReader> &&chr_names,
        std::unique_ptr<StringReader> &&cell_names) :
    cell(std::move(cell)), start(std::move(start)), end(std::move(end)), end_max(std::move(end_max)),
    chr_ptr(std::move(chr_ptr)), chr_names(std::move(chr_names)), cell_names(std::move(cell_names)) {}

// Read end_max_buf from end_max iterator, making it equal to the values between
// start_idx and end_idx
void UnpackedFragments3::readEndMaxBuf(uint32_t start_idx, uint32_t end_idx) {
    // Read the end_max buffer
    end_max_buf.resize(end_idx/128 - start_idx/128 + 1);
    end_max->seek(start_idx/128);
    uint32_t i = 0;
    while (true) {
        uint32_t load_amount = std::min(end_max->capacity(), (uint32_t) end_max_buf.size() - i);
        if (load_amount == 0)
            throw std::runtime_error("Failed to load end_max data");
        std::memmove(&end_max_buf[i], end_max->data(), load_amount * sizeof(uint32_t));
        i += load_amount;
        if (i >= end_max_buf.size()) break;
        end_max->next();
    }
    // Handle starting a chromosome on a non-multiple of 128, such end_max[0]
    // may be the max of the prior chromosome rather than max of the new one
    if (start_idx % 128 != 0 && end_max_buf.size() > 1) {
        end_max_buf[0] = std::min(end_max_buf[0], end_max_buf[1]);
    }
}

bool UnpackedFragments3::isSeekable() const {return true;}
void UnpackedFragments3::seek(uint32_t chr_id, uint32_t base) {
    if (chr_id != current_chr) {
        current_chr = chr_id;
        chr_ptr->seek(chr_id*2);
        chr_start_ptr = chr_ptr->read_one();
        chr_end_ptr = chr_ptr->read_one();
        
        readEndMaxBuf(chr_start_ptr, chr_end_ptr);
    }

    // Binary search for base in end_max
    uint32_t current_block = chr_start_ptr/128 + 
        std::lower_bound(end_max_buf.begin(), end_max_buf.end(), base) - 
        end_max_buf.begin();

    // Add this max in case we're seeking to the first block in a chromosome that
    // doesn't start on a multiple of 128
    current_idx = std::max(chr_start_ptr, current_block*128);

    start->seek(current_idx);
    end->seek(current_idx);
    cell->seek(current_idx);
    
}

void UnpackedFragments3::restart() {
    current_chr = UINT32_MAX;
    current_idx = UINT32_MAX;

}

int UnpackedFragments3::chrCount() const {return chr_names->size();}
int UnpackedFragments3::cellCount() const {return cell_names->size();}

const char* UnpackedFragments3::chrNames(uint32_t chr_id) const {
    return chr_names->get(chr_id);
}
const char* UnpackedFragments3::cellNames(uint32_t cell_id) const {
    return cell_names->get(cell_id); 
}

bool UnpackedFragments3::nextChr() {
    current_chr += 1;
    if (current_chr >= chrCount()) {
        current_chr -= 1;
        return false;
    }

    chr_start_ptr = chr_ptr->read_one();
    chr_end_ptr = chr_ptr->read_one();
    // Check if we need to perform any seeks, or if we're all good
    // since we just finished a chromosome
    if (current_idx != chr_start_ptr) {
        // We need to perform seeks to get to the right data
        // reading location
        start->seek(chr_start_ptr);
        end->seek(chr_start_ptr);
        cell->seek(chr_start_ptr);
    }
    current_idx = chr_start_ptr;
    readEndMaxBuf(chr_start_ptr, chr_end_ptr);
    return true;
}
uint32_t UnpackedFragments3::currentChr() const {return current_chr;}

int32_t UnpackedFragments3::load(uint32_t count, FragmentArray buffer) {
    if (current_idx >= chr_end_ptr) {
        return 0;
    }
    // Load cell, start, or end if necessary
    if (cell->capacity() == 0 && !cell->next()) {throw std::runtime_error("Failed to load next cell data");}
    if (start->capacity() == 0 && !start->next()) throw std::runtime_error("Failed to load next start data");
    if (end->capacity() == 0 && !end->next()) throw std::runtime_error("Failed to load next end data");
    
    uint32_t load_a = std::min(cell->capacity(), start->capacity());
    uint32_t load_b = std::min(count, chr_end_ptr - current_idx);
    uint32_t load = std::min(load_a, end->capacity());
    load = std::min(load, load_b);

    current_idx += load;
    std::memmove(buffer.cell, cell->data(), load*sizeof(uint32_t));
    std::memmove(buffer.start, start->data(), load*sizeof(uint32_t));
    std::memmove(buffer.end, end->data(), load*sizeof(uint32_t));

    cell->advance(load);
    start->advance(load);
    end->advance(load);
    return load;
}


UnpackedFragmentsWriter3::UnpackedFragmentsWriter3(WriterPtr &&cell, WriterPtr &&start,
     WriterPtr &&end, WriterPtr &&end_max, WriterPtr && chr_ptr,
     std::unique_ptr<StringWriter> &&chr_names, std::unique_ptr<StringWriter> &&cell_names) :
     cell(std::move(cell)), start(std::move(start)), end(std::move(end)), 
     end_max(std::move(end_max)), chr_ptr(std::move(chr_ptr)),
     chr_names(std::move(chr_names)), cell_names(std::move(cell_names)) {}

bool UnpackedFragmentsWriter3::write(FragmentsIterator &fragments, void (*checkInterrupt)(void)) {
    uint32_t cur_end_max = 0;
    uint32_t prev_end_max = 0;
    uint32_t idx = 0;

    std::vector<uint32_t> chr_ptr_buf;

    while (fragments.nextChr()) {
        uint32_t chr_id = fragments.currentChr();
        if (chr_id*2 + 2 >= chr_ptr_buf.size()) {
            chr_ptr_buf.resize(chr_id*2 + 2);
        }
        chr_ptr_buf[chr_id*2] = idx;

        // This helps handle the case of multiple chromosomes ending within
        // the same block of 128 fragments
        prev_end_max = std::max(prev_end_max, cur_end_max); 
        cur_end_max = 0;
        bool chr_done = false;
        while (!chr_done) {
            const uint32_t read_batch = 128 - idx % 128;
            cell->ensureCapacity(read_batch);
            start->ensureCapacity(read_batch);
            end->ensureCapacity(read_batch);
            uint32_t i;
            for (i = 0; i < read_batch; i++) {
                if (!fragments.nextFrag()) {
                    chr_done = true;
                    break;
                }
                cell->data()[i] = fragments.cell();
                start->data()[i] = fragments.start();
                end->data()[i] = fragments.end();
                cur_end_max = std::max(cur_end_max, fragments.end());
            }
            cell->advance(i);
            start->advance(i);
            end->advance(i);
            idx += i;

            // This handles overhang if a block of 128 crosses a chromosome boundary
            if (idx % 128 == 0) {
                end_max->write_one(std::max(cur_end_max, prev_end_max));
                prev_end_max = 0;
            }

            if (idx % 8192 == 0) checkInterrupt();
        }
        chr_ptr_buf[chr_id*2+1] = idx;
    }
    if (idx % 128 != 0) {
        end_max->write_one(std::max(cur_end_max, prev_end_max));
    }

    chr_ptr->ensureCapacity(chr_ptr_buf.size());
    std::memmove(chr_ptr->data(), chr_ptr_buf.data(), chr_ptr_buf.size() * sizeof(uint32_t));
    chr_ptr->advance(chr_ptr_buf.size());

    cell->backup(cell->capacity()); cell.reset(nullptr);
    start->backup(start->capacity()); start.reset(nullptr);
    end->backup(end->capacity()); end.reset(nullptr);
    end_max->backup(end_max->capacity()); end_max.reset(nullptr);
    chr_ptr->backup(chr_ptr->capacity()); chr_ptr.reset(nullptr);

    // Get cell and chromosome names. This probably incurs a few extra copies,
    // but it shouldn't matter since writing the actual fragments should dominate cost
    std::vector<std::string> cell_names;
    for (int i = 0; ; i++) {
        const char* cell_name = fragments.cellNames(i);
        if (cell_name == NULL) break;
        cell_names.push_back(std::string(cell_name));
    }
    this->cell_names->write(VecStringReader(cell_names));
    
    std::vector<std::string> chr_names;
    for (int i = 0; ; i++) {
        const char* chr_name = fragments.chrNames(i);
        if (chr_name == NULL) break;
        chr_names.push_back(std::string(chr_name));
    }
    this->chr_names->write(VecStringReader(chr_names));
    return true;
}

} // end namespace BPCells