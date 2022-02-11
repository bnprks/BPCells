#pragma once

#include <algorithm>

#include "../arrayIO/frags.h"
#include "../bitpacking/packing_utils.h"
#include "FragmentsIterator.h"


namespace BPCells {



template<class Loader>
class UnpackedFragments2 : public FragmentsLoader {
public:
    // Construct by passing the arguments to construct the underlying Loader class
    UnpackedFragments2(Loader&& loader) :
        loader(std::move(loader)), 
        chr_names(this->loader.readChrNames()), cell_names(this->loader.readCellNames()), 
        current_chr(UINT32_MAX) {}
    
    UnpackedFragments2() = delete;

    bool isSeekable() const override {return true;}
    // Move iterator to just before fragments which end after "base".
    // It's possible that fragments returned after seek will end before "base",
    // but it's guaranteed that it won't skip over any fragments ending before "base"
    void seek(uint32_t chr_id, uint32_t base) override {
        // Reload end_max buffer if needed
        if (chr_id != current_chr) {
            current_chr = chr_id;
            
            initChr();
        }
        // Binary search for base in end_max
        uint32_t current_block = std::lower_bound(end_max_buf.begin(), end_max_buf.end(), base) - 
            end_max_buf.begin();

        chr_frags.start.seek(current_block * 128);
        chr_frags.end.seek(current_block * 128);
        chr_frags.cell.seek(current_block * 128);
    };

    // Reset the iterator to start from the beginning
    void restart() override {
        current_chr = UINT32_MAX;
    }

    // Return the number of cells/chromosomes, or return -1 if this number is 
    // not known ahead of time
    int chrCount() const override { return chr_names.size(); }
    int cellCount() const override { return cell_names.size(); }

    // Return name for a given chr_id or cell_id. Only valid to call
    // for chromosme or cell_ids that have been actually returned by the iterator
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
        chr_frags.start.seek(0);
        chr_frags.end.seek(0);
        chr_frags.cell.seek(0);
        return true;
    }

    // Return chromosome ID of current fragment
    uint32_t currentChr() const override {return current_chr; }

    // Return number of items loaded. Should repeatedly return 0 at the end of a chromosome.
    // Return -1 for error
    int32_t load(uint32_t count, FragmentArray buffer) override {
        uint32_t ret = chr_frags.start.read(buffer.start, count);
        if (ret != chr_frags.end.read(buffer.end, count))
            return -1;
        if (ret != chr_frags.cell.read(buffer.cell, count))
            return -1;
        return ret;
    }
private:
    Loader loader;
    const std::vector<std::string> chr_names, cell_names;
    std::vector<uint32_t> end_max_buf;
    uint32_t chr_frags_count;
    uint32_t current_chr;
    UnpackedFrags<typename Loader::UIntReader> chr_frags;
    
    void initChr() {
        // Load the end_max buf into memory to enable seeking, and load the number of
        // fragments for this chromosome
        chr_frags = std::move(loader.chrReaderUnpacked(current_chr));
        size_t size = chr_frags.end_max.size();
        end_max_buf.resize(size);

        chr_frags.end_max.seek(0);
        chr_frags.end_max.read(end_max_buf.data(), size);
        chr_frags_count = chr_frags.start.size();
    }
};


template<class Saver>
class UnpackedFragmentsWriter2 : public FragmentsWriter {
public:
    using UIntWriter = typename Saver::UIntWriter;

    // Construct by passing the arguments to construct the underlying Saver class
    UnpackedFragmentsWriter2(Saver &&saver) : saver(std::move(saver)) {}
    UnpackedFragmentsWriter2(const UnpackedFragmentsWriter2&) = delete;
    UnpackedFragmentsWriter2& operator=(const UnpackedFragmentsWriter2& other) = delete;

    bool write(FragmentsIterator &fragments, void (*checkInterrupt)(void) = NULL) override {
        uint32_t start_buffer[128];
        uint32_t end_buffer[128];
        uint32_t cell_buffer[128];

        while (fragments.nextChr()) {
            uint32_t i;
            uint32_t j;

            uint32_t chr_id = fragments.currentChr();

            // TODO: Assert that we aren't re-creating the same chromosome we've already seen?
            UnpackedFrags<UIntWriter> chr_frags = saver.chrWriterUnpacked(chr_id);

            bool chr_done = false;
            uint32_t end_max = 0;
            while (!chr_done) {
                // Write in batches of 128 ints, and check for an interrupt
                // after every 100 batches
                for (j = 0; j < 100 && !chr_done; j++) {
                    if (checkInterrupt != NULL) checkInterrupt();
                    for (i = 0; i < 128; i++) {
                        if (!fragments.nextFrag()) {
                            chr_done = true;
                            break;
                        }
                        start_buffer[i] = fragments.start();
                        end_buffer[i] = fragments.end();
                        cell_buffer[i] = fragments.cell();
                        end_max = std::max(end_max, fragments.end());
                    }

                    chr_frags.start.write(start_buffer, i);
                    chr_frags.end.write(end_buffer, i);
                    chr_frags.cell.write(cell_buffer, i);
                    chr_frags.end_max.write(&end_max, 1);
                }
            }   
            chr_frags.start.finalize();
            chr_frags.end.finalize();
            chr_frags.cell.finalize();
            chr_frags.end_max.finalize();
        }
        

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
    Saver saver;
};



} // end namespace BPCells