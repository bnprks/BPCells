#include "PackedFragments3.h"

namespace BPCells {

PackedFragments3::PackedFragments3(ReaderPtr &&cell_data, ReaderPtr &&cell_idx, ReaderPtr &&start_data, ReaderPtr &&start_idx, ReaderPtr &&start_starts, 
        ReaderPtr &&end_data, ReaderPtr &&end_idx, ReaderPtr &&end_max, ReaderPtr &&chr_ptr, 
        std::unique_ptr<StringReader> &&chr_names,
        std::unique_ptr<StringReader> &&cell_names) :
            cell_data(std::move(cell_data)),
            cell_idx(std::move(cell_idx)),
            start_data(std::move(start_data)),
            start_idx(std::move(start_idx)),
            start_starts(std::move(start_starts)),
            end_data(std::move(end_data)),
            end_idx(std::move(end_idx)),
            end_max(std::move(end_max)),
            chr_ptr(std::move(chr_ptr)),
            chr_names(std::move(chr_names)), 
            cell_names(std::move(cell_names)) {}

// Read end_max_buf from end_max iterator, making it equal to the values between
// start_idx and end_idx
void PackedFragments3::readEndMaxBuf(uint32_t start_idx, uint32_t end_idx) {
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

void PackedFragments3::load128(uint32_t *cell_out, uint32_t *start_out, uint32_t *end_out) {
    uint32_t next_cell_idx = cell_idx->read_one();
    uint32_t next_start_idx = start_idx->read_one();
    uint32_t start_start = start_starts->read_one();
    uint32_t next_end_idx = end_idx->read_one();

    // Unpack cell
    uint32_t bits_cell = (next_cell_idx - prev_cell_idx) / 4;
    cell_data->ensureCapacity(bits_cell * 4);
    simdunpack(cell_data->data(), cell_out, bits_cell);
    cell_data->advance(bits_cell * 4);
    prev_cell_idx = next_cell_idx;

    // Unpack start
    uint32_t bits_start = (next_start_idx - prev_start_idx) / 4;
    start_data->ensureCapacity(bits_start * 4);
    simdunpackd1(start_start, start_data->data(), start_out, bits_start);
    start_data->advance(bits_start * 4);
    prev_start_idx = next_start_idx;

    // Unpack end
    uint32_t bits_end = (next_end_idx - prev_end_idx) / 4;
    end_data->ensureCapacity(bits_end * 4);
    simdunpack(end_data->data(), end_out, bits_end);
    end_data->advance(bits_end * 4);
    simdadd(end_out, start_out);
    prev_end_idx = next_end_idx;
}

bool PackedFragments3::isSeekable() const {return true;}
void PackedFragments3::seek(uint32_t chr_id, uint32_t base) {
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

    cell_idx->seek(current_block);
    start_idx->seek(current_block);
    end_idx->seek(current_block);

    prev_cell_idx = cell_idx->read_one();
    prev_start_idx = start_idx->read_one();
    prev_end_idx = end_idx->read_one();

    cell_data->seek(prev_cell_idx);
    start_data->seek(prev_start_idx);
    end_data->seek(prev_end_idx);

    if (current_idx % 128 != 0) {
        load128(cell_buf, start_buf, end_buf);
    }
}

void PackedFragments3::restart() {
    current_chr = UINT32_MAX;
    current_idx = UINT32_MAX;

}

int PackedFragments3::chrCount() const {return chr_names->size();}
int PackedFragments3::cellCount() const {return cell_names->size();}

const char* PackedFragments3::chrNames(uint32_t chr_id) const {
    return chr_names->get(chr_id);
}
const char* PackedFragments3::cellNames(uint32_t cell_id) const {
    return cell_names->get(cell_id); 
}

bool PackedFragments3::nextChr() {
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
        cell_idx->seek(chr_start_ptr / 128);
        start_idx->seek(chr_start_ptr / 128);
        end_idx->seek(chr_start_ptr / 128);

        prev_cell_idx = cell_idx->read_one();
        prev_start_idx = start_idx->read_one();
        prev_end_idx = end_idx->read_one();

        cell_data->seek(prev_cell_idx);
        start_data->seek(prev_start_idx);
        end_data->seek(prev_end_idx);

        if (chr_start_ptr % 128 != 0) {
            load128(cell_buf, start_buf, end_buf);
        }
    }
    current_idx = chr_start_ptr;
    readEndMaxBuf(chr_start_ptr, chr_end_ptr);
    return true;
}
uint32_t PackedFragments3::currentChr() const {return current_chr;}

int32_t PackedFragments3::load(uint32_t count, FragmentArray buffer) {
    if (current_idx >= chr_end_ptr) {
        return 0;
    }
    if (count < 128) throw std::runtime_error("Must load >128 entries at a time from PackedFragments");    

    uint32_t i = 0;
    // Check if we are starting a column and need to copy from the buffer
    if (current_idx % 128 != 0) {
        uint32_t copy_num = std::min(128 - current_idx % 128, chr_end_ptr - current_idx);
        std::memmove(buffer.cell, &cell_buf[current_idx % 128], copy_num * sizeof(uint32_t));
        std::memmove(buffer.start, &start_buf[current_idx % 128], copy_num * sizeof(uint32_t));
        std::memmove(buffer.end, &end_buf[current_idx % 128], copy_num * sizeof(uint32_t));
        i += copy_num;
        current_idx += copy_num;
    }

    for (; i + 128 <= count && chr_end_ptr - current_idx >= 128; i+= 128) {
        load128(&buffer.cell[i], &buffer.start[i], &buffer.end[i]);
        current_idx += 128;
    }

    if (i < count && current_idx < chr_end_ptr) {
        // Save any straglers into the buffer so we can read the rest of the
        // block of 128 on the subsequent load call
        uint32_t remaining = std::min(count-i, chr_end_ptr-current_idx);
        load128(cell_buf, start_buf, end_buf);
        std::memmove(&buffer.cell[i], cell_buf, sizeof(uint32_t) * remaining);
        std::memmove(&buffer.start[i], start_buf, sizeof(uint32_t) * remaining);
        std::memmove(&buffer.end[i], end_buf, sizeof(uint32_t) * remaining);
        i += remaining;
        current_idx += remaining;
    }
    return i;
}


PackedFragmentsWriter3::PackedFragmentsWriter3(WriterPtr &&cell_data, WriterPtr &&cell_idx, WriterPtr &&start_data, WriterPtr &&start_idx, WriterPtr &&start_starts, 
        WriterPtr &&end_data, WriterPtr &&end_idx, WriterPtr &&end_max, WriterPtr &&chr_ptr, 
        std::unique_ptr<StringWriter> &&chr_names,
        std::unique_ptr<StringWriter> &&cell_names) :
            cell_data(std::move(cell_data)),
            cell_idx(std::move(cell_idx)),
            start_data(std::move(start_data)),
            start_idx(std::move(start_idx)),
            start_starts(std::move(start_starts)),
            end_data(std::move(end_data)),
            end_idx(std::move(end_idx)),
            end_max(std::move(end_max)),
            chr_ptr(std::move(chr_ptr)),
            chr_names(std::move(chr_names)),
            cell_names(std::move(cell_names)) {}

void PackedFragmentsWriter3::pack128(const uint32_t* cell_in, const uint32_t* start_in, uint32_t* end_in,
                uint32_t &cur_cell_idx, uint32_t &cur_start_idx, uint32_t &cur_end_idx) {

    // Note: assumes cell_data, start_data and end_data have sufficient capacity
    // Note: Will modify the input end data, transforming it into fragment widths

    // Pack cell
    uint32_t cell_bits = simdmaxbits(cell_in);
    simdpack(cell_in, cell_data->data(), cell_bits);
    cell_data->advance(cell_bits*4);
    cur_cell_idx += cell_bits*4;
    cell_idx->write_one(cur_cell_idx);

    // Pack end (must do this before packing start, since the current setup has
    // the packed start data written to the begin of start_out, thus corrupting it)
    simdsubtract(end_in, start_in);
    uint32_t end_bits = simdmaxbits(end_in);
    simdpack(end_in, end_data->data(), end_bits);
    end_data->advance(end_bits*4);
    cur_end_idx += end_bits*4;
    end_idx->write_one(cur_end_idx);

    // Pack start
    uint32_t start_bits = simdmaxbitsd1(start_in[0], start_in);
    start_starts->write_one(start_in[0]);
    simdpackd1(start_in[0], start_in, start_data->data(), start_bits);
    start_data->advance(start_bits*4);
    cur_start_idx += start_bits*4;
    start_idx->write_one(cur_start_idx);

    

}
bool PackedFragmentsWriter3::write(FragmentsIterator &fragments, void (*checkInterrupt)(void)) {
    uint32_t cur_cell_idx = 0;
    uint32_t cur_start_idx = 0;
    uint32_t cur_end_idx = 0;
    
    uint32_t cur_end_max = 0;
    uint32_t prev_end_max = 0;
    uint32_t idx = 0;

    uint32_t loaded = 0; // Number of values loaded most recently
    std::vector<uint32_t> chr_ptr_buf;

    uint32_t cell_buffer[128], start_buffer[128], end_buffer[128];

    cell_idx->write_one(0);
    start_idx->write_one(0);
    end_idx->write_one(0);

    while (fragments.nextChr()) {
        if (checkInterrupt != NULL) checkInterrupt();

        uint32_t chr_id = fragments.currentChr();
        if (chr_id*2 + 2 >= chr_ptr_buf.size()) {
            chr_ptr_buf.resize(chr_id*2 + 2);
        }
        chr_ptr_buf[chr_id*2] = idx;

        // This helps handle the case of multiple chromosomes ending within
        // the same block of 128 fragments
        prev_end_max = std::max(prev_end_max, cur_end_max); 
        cur_end_max = 0;

        while (true) {
            // Strategy: Load up to 256 entries into buffers. 
            // If we end up getting a non-multiple of 128 values to compress, then try to make our next
            // load continue where the last one left off.
            // However, if we get to a point where there is space for less than 128 more values,
            // then copy the leftovers to the beginning of the buffer and continue.

            // Writing strategy: Use the output data buffers
            // in order to load >= 256 matrix entries prior to actually compressing + writing them.
            // We will need our own buffer space for 128 values to store stragglers
            uint32_t capacity = std::min(cell_data->capacity(), start_data->capacity());
            capacity = std::min(capacity, end_data->capacity());
            
            uint32_t *cell_data_buf = cell_data->data();
            uint32_t *start_data_buf = start_data->data();
            uint32_t *end_data_buf = end_data->data();
            uint32_t i;

            for (i = 0; i + 128 <= capacity; i += 128) {
                pack128(cell_data_buf + i, start_data_buf + i, end_data_buf + i,
                        cur_cell_idx, cur_start_idx, cur_end_idx);
            }
            uint32_t hangover = capacity % 128;
            // If our remaining capacity is not a multiple of 128, then stash the data before calling ensureCapacity
            if (hangover != 0) {
                std::memmove(cell_buffer, cell_data_buf + i, sizeof(uint32_t)*hangover);
                std::memmove(start_buffer, start_data_buf + i, sizeof(uint32_t)*hangover);
                std::memmove(end_buffer, end_data_buf + i, sizeof(uint32_t)*hangover);
            }

            cell_data->ensureCapacity(256);
            start_data->ensureCapacity(256);
            end_data->ensureCapacity(256);
            capacity = std::min(cell_data->capacity(), start_data->capacity());
            capacity = std::min(capacity, end_data->capacity());

            if (hangover != 0) {
                std::memmove(cell_data->data(), cell_buffer, sizeof(uint32_t)*hangover);
                std::memmove(start_data->data(), start_buffer, sizeof(uint32_t)*hangover);
                std::memmove(end_data->data(), end_buffer, sizeof(uint32_t)*hangover);
            }

            if (checkInterrupt != NULL) checkInterrupt();

            loaded = fragments.load(
                capacity-hangover, 
                {start_data->data()+hangover, 
                end_data->data()+hangover,
                cell_data->data()+hangover, 
                (uint32_t) capacity-hangover}
            );

            // Calculate & output end_max for all of the loaded values
            // Only count max for the actually loaded values
            for (i = 0; i < loaded; i++) {
                cur_end_max = std::max(
                        cur_end_max, 
                        end_data->data()[hangover + i]);
                if ((idx + i) % 128 == 127) {
                    end_max->write_one(std::max(cur_end_max, prev_end_max));
                    prev_end_max = 0;
                }
            }
            
            idx += loaded;

            cell_data->backup(cell_data->capacity() - loaded - hangover);
            start_data->backup(start_data->capacity() - loaded - hangover);
            end_data->backup(end_data->capacity() - loaded - hangover);
            if (loaded == 0) break;
            // END PACKED MATRIX LOOP COPY
        }
        chr_ptr_buf[chr_id*2+1] = idx;
    }
    if (idx % 128 != 0) {
        end_max->write_one(std::max(cur_end_max, prev_end_max));
    }
    // Finish compression by taking care of any leftovers
    // capacity should == hangover from above
    uint32_t capacity = std::min(cell_data->capacity(), start_data->capacity());
    capacity = std::min(capacity, end_data->capacity());
    if (capacity != 0) {
        std::memmove(cell_buffer, cell_data->data(), sizeof(uint32_t)*capacity);
        std::memmove(start_buffer, start_data->data(), sizeof(uint32_t)*capacity);
        std::memmove(end_buffer, end_data->data(), sizeof(uint32_t)*capacity);
        for (uint32_t i = capacity; i % 128 != 0; i++) {
            cell_buffer[i] = cell_buffer[i-1];
            start_buffer[i] = start_buffer[i-1];
            end_buffer[i] = end_buffer[i-1];
        }
        cell_data->ensureCapacity(128);
        start_data->ensureCapacity(128);
        end_data->ensureCapacity(128);
        pack128(cell_buffer, start_buffer, end_buffer,
                        cur_cell_idx, cur_start_idx, cur_end_idx);
    }

    chr_ptr->ensureCapacity(chr_ptr_buf.size());
    std::memmove(chr_ptr->data(), chr_ptr_buf.data(), chr_ptr_buf.size() * sizeof(uint32_t));
    chr_ptr->advance(chr_ptr_buf.size());

    cell_data->backup(cell_data->capacity()); cell_data.reset(nullptr);
    cell_idx->backup(cell_idx->capacity()); cell_idx.reset(nullptr);
    start_data->backup(start_data->capacity()); start_data.reset(nullptr);
    start_idx->backup(start_idx->capacity()); start_idx.reset(nullptr);
    start_starts->backup(start_starts->capacity()); start_starts.reset(nullptr);
    end_data->backup(end_data->capacity()); end_data.reset(nullptr);
    end_idx->backup(end_idx->capacity()); end_idx.reset(nullptr);
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