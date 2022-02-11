#pragma once

#include <algorithm>

#include "../bitpacking/packing_utils.h"
#include "FragmentsIterator.h"

namespace BPCells {



template<class Loader>
class PackedFragments2 : public FragmentsLoader {
public:
    PackedFragments2(Loader&& loader, 
        const std::vector<std::string> &cell_names, const std::vector<std::string> &chr_names) :
        loader(std::move(loader)), chr_names(chr_names), cell_names(cell_names), current_chr(UINT32_MAX) {}
    PackedFragments2(Loader&& loader) :
        loader(std::move(loader)), 
        chr_names(this->loader.readChrNames()), cell_names(this->loader.readCellNames()), 
        current_chr(UINT32_MAX) {}
    
    PackedFragments2() = delete;

    bool isSeekable() const override {return true;}
    
    void seek(uint32_t chr_id, uint32_t base) override {
        // Reload end_max buffer if needed
        if (chr_id != current_chr) {
            current_chr = chr_id;
            
            initChr();
        }
        // Binary search for base in end_max
        current_chunk = std::lower_bound(end_max_buf.begin(), end_max_buf.end(), base) - 
            end_max_buf.begin();
        
        // Invariants:
        // start_idx, end_idx, cell_idx holds the index of the
        // *next* data block to read.
        // chr_frags.start_data, end_data, and cell_data should
        // be seeked to that next block
        
        chr_frags.start_idx.seek(current_chunk);
        chr_frags.start_starts.seek(current_chunk);
        chr_frags.end_idx.seek(current_chunk);
        chr_frags.cell_idx.seek(current_chunk);

        chr_frags.start_idx.read(&start_idx, 1);
        chr_frags.end_idx.read(&end_idx, 1);
        chr_frags.cell_idx.read(&cell_idx, 1);

        chr_frags.start_data.seek(start_idx);
        chr_frags.end_data.seek(end_idx);
        chr_frags.cell_data.seek(cell_idx);
    };

    void restart() override {
        current_chr = UINT32_MAX;
    }

    int chrCount() const override { return chr_names.size(); }
    int cellCount() const override { return cell_names.size(); }

    const char* chrNames(uint32_t chr_id) const override { 
        if (chr_id >= chr_names.size()) return NULL;
        return chr_names[chr_id].c_str(); 
    }
    const char* cellNames(uint32_t cell_id) const override { 
        if (cell_id >= cell_names.size()) return NULL;
        return cell_names[cell_id].c_str(); 
    }

    // Advance the iterator to the next chromosome. Return false if there are no more chromosomes
    bool nextChr() override {
        current_chr++;
        if (current_chr >= chr_names.size()) {
            current_chr -= 1;
            return false;
        }
        initChr();

        chr_frags.start_idx.seek(0);
        chr_frags.start_starts.seek(0);
        chr_frags.end_idx.seek(0);
        chr_frags.cell_idx.seek(0);

        chr_frags.start_idx.read(&start_idx, 1);
        chr_frags.end_idx.read(&end_idx, 1);
        chr_frags.cell_idx.read(&cell_idx, 1);

        chr_frags.start_data.seek(0);
        chr_frags.end_data.seek(0);
        chr_frags.cell_data.seek(0);

        current_chunk = 0;

        return true;
    }

    // Return chromosome ID of current fragment
    uint32_t currentChr() const override { return current_chr; }
private:
    Loader loader;
    const std::vector<std::string> chr_names, cell_names;
    std::vector<uint32_t> end_max_buf;
    uint32_t chr_frags_count;
    uint32_t current_chr;
    uint32_t current_chunk;
    PackedFrags<typename Loader::UIntReader> chr_frags;
    // Invariants:
    // start_idx, end_idx, cell_idx hold the index of the
    // *next* data block to read.
    // frags[current_chr].start_data, end_data, and cell_data should
    // be seeked to that next block
    uint32_t start_idx, end_idx, cell_idx;

    void initChr() {
        chr_frags = std::move(loader.chrReaderPacked(current_chr));
        size_t size = chr_frags.end_max.size();
        end_max_buf.resize(size);

        chr_frags.end_max.seek(0);
        chr_frags.end_max.read(end_max_buf.data(), size);

        chr_frags.count.seek(0);
        chr_frags.count.read(&chr_frags_count, 1);
    }

    int32_t load(uint32_t count, FragmentArray buffer) override {
        if (current_chunk == end_max_buf.size()) return 0;
        if (count < 128) throw std::runtime_error("Must load >128 fragments at a time from PackedFragments");
        count = count & ~(128-1); // Round count down to nearest 128
        
        uint32_t next_start_idx, next_end_idx, next_cell_idx;
        uint32_t start_start, data_buf[128];

        for (int i = 0; i < count; i += 128) {                        
            // Read the next data offsets for calculating bit widths
            chr_frags.start_idx.read(&next_start_idx, 1);
            chr_frags.end_idx.read(&next_end_idx, 1);
            chr_frags.cell_idx.read(&next_cell_idx, 1);
            current_chunk += 1;

            // Unpack start
            uint32_t bits_start = (next_start_idx - start_idx) / 4;
            chr_frags.start_data.read(data_buf, bits_start * 4);
            chr_frags.start_starts.read(&start_start, 1);
            simdunpackd1(start_start, data_buf, &buffer.start[i], bits_start);

            // Unpack end
            uint32_t bits_end = (next_end_idx - end_idx) / 4;
            chr_frags.end_data.read(data_buf, bits_end * 4);
            simdunpack(data_buf, &buffer.end[i], bits_end);
            simdadd(&buffer.end[i], &buffer.start[i]);

            // Unpack cell
            uint32_t bits_cell = (next_cell_idx - cell_idx) / 4;
            chr_frags.cell_data.read(data_buf, bits_cell * 4);
            // Original design considered using FOR packing for cell_ids.
            // Since we are currently just storing pos-sorted without
            // much cell chunking I'll disable the FOR encoding
            simdunpack(data_buf, &buffer.cell[i], bits_cell);

            start_idx = next_start_idx;
            end_idx = next_end_idx;
            cell_idx = next_cell_idx;

            if (current_chunk == end_max_buf.size()) {
                uint32_t overhang = (128-chr_frags_count%128) % 128;
                return i + 128 - overhang;
            };
        }
        return count;
    }
};

// Design of PackedFragmentsWriter2
// The template helper is *just* for writing the integer data, with out worrying
// about storage for cell + chr names yet.
// One pointy bit is figuring out how many fragments are actually written per chr
// Another pointy bit is being able to construct new PackedFrags of the appropriate type
// (i.e. initialize the correct hdf5 groups, or create the right files, or construct in-memory vectors)

template<class Saver>
class PackedFragmentsWriter2 : public FragmentsWriter {
    Saver saver;
public:
    PackedFragmentsWriter2(Saver&& saver) : saver(std::move(saver)) {}
    PackedFragmentsWriter2(const PackedFragmentsWriter2&) = delete;
    PackedFragmentsWriter2& operator=(const PackedFragmentsWriter2& other) = delete;

    bool write(FragmentsIterator &fragments, void (*checkInterrupt)(void) = NULL) override {
        uint32_t start_bits, end_bits, cell_id_bits;
        uint32_t start_idx, end_idx, cell_id_idx;
        uint32_t start_buffer[128];
        uint32_t end_buffer[128];
        uint32_t tmp_buffer[128];
        uint32_t cell_buffer[128];
        uint32_t chr_id, i;
        uint32_t chr_frags_count;
        uint32_t end_max;
        bool chr_done;

        while (fragments.nextChr()) {
            chr_id = fragments.currentChr();
            // TODO: Assert that we aren't re-creating the same chromosome we've already seen?
            auto chr_frags = saver.chrWriterPacked(chr_id);

            start_idx = 0;
            end_idx = 0;
            cell_id_idx = 0;
            chr_frags.start_idx.write(&start_idx, 1);
            chr_frags.end_idx.write(&end_idx, 1);
            chr_frags.cell_idx.write(&cell_id_idx, 1);

            chr_frags_count = 0;
            end_max = 0;

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
                
                chr_frags_count += i;    

                if (checkInterrupt != NULL && (chr_frags_count / 128) % 1024 == 0) checkInterrupt();

                // Cleanup loop in case we're at end of chromosome and
                // couldn't read all 128
                for (;i < 128; i++) {
                    start_buffer[i] = start_buffer[i-1];
                    end_buffer[i] = end_buffer[i-1];
                    cell_buffer[i] = cell_buffer[i-1];
                }
                
                // Pack starts
                start_bits = simdmaxbitsd1(start_buffer[0], start_buffer);
                chr_frags.start_starts.write(&start_buffer[0], 1);
                simdpackd1(start_buffer[0], start_buffer, tmp_buffer, start_bits);
                chr_frags.start_data.write(tmp_buffer, start_bits * 4);
                start_idx += start_bits*4;
                chr_frags.start_idx.write(&start_idx, 1);

                // Pack ends
                end_max = std::max(simdmax(end_buffer), end_max);
                chr_frags.end_max.write(&end_max, 1);
                simdsubtract(end_buffer, start_buffer);
                end_bits = simdmaxbits(end_buffer);

                simdpack(end_buffer, tmp_buffer, end_bits);
                chr_frags.end_data.write(tmp_buffer, end_bits*4);
                end_idx += end_bits*4;
                chr_frags.end_idx.write(&end_idx, 1);

                // Pack cell_ids
                cell_id_bits = simdmaxbits(cell_buffer);
                simdpack(cell_buffer, tmp_buffer, cell_id_bits);
                chr_frags.cell_data.write(tmp_buffer, cell_id_bits * 4);
                cell_id_idx += cell_id_bits*4;
                chr_frags.cell_idx.write(&cell_id_idx, 1);
            } // End loop over fragments in chromosome
            chr_frags.count.write(&chr_frags_count, 1);
        } // End loop over chromosomes

        // Get cell and chromosome names
        std::vector<std::string> cell_names;
        for (int i = 0; ; i++) {
            const char* cell_name = fragments.cellNames(i);
            if (cell_name == NULL) break;
            cell_names.push_back(cell_name);
        }
        saver.writeCellNames(cell_names);
        std::vector<std::string> chr_names;
        for (int i = 0; ; i++) {
            const char* chr_name = fragments.chrNames(i);
            if (chr_name == NULL) break;
            chr_names.push_back(chr_name);
        }
        saver.writeChrNames(chr_names);
        return true;
    }
private:
    std::vector<PackedFrags<typename Saver::UIntWriter> > frags;
};

} // end namespace BPCells