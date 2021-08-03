#pragma once
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace BPCells {

struct FragmentArray {
    uint32_t *start, *end, *cell;
    uint32_t capacity;
};

// Interface for loading sorted fragments from files or memory.
// Main workhorse is the load function, which copies fragments into a buffer.
// Users will call load repeatedly, until it hits the end of a chromosome or the end
// of a file. 
// Optionally, a FragmentLoader can support seeking, which skips the current
// reading position to a given genomic coordinate.
// Requirements for building a FragmentLoader:
// 1. Coordinates are assumed to be 0-based and half-open, same as the bed file format. 
//    This means the first base in a chromosome is at position 0, and the end coordinate
//    of a fragment is one basepair after the last mapped base.
//    For more discussion, see here: http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/
// 2. Assign chromosme and cell IDs that count up from zero with no skipping. It's
//    okay if a given ID is never emitted, but there must be a string ID which can
//    be returned from chrNames or cellNames
// 3. The number of valid chromosomes and cell IDs must be known in advance for
//    seekable FragmentIterators, along with their corresponding names
// 4. Fragments are grouped by chromosme, although they need not be covered in
//    order of chromosome ID.
// 5. Implement load function, which will always load >= 1 fragment unless there is
//    an error or we have reached the end of a chromosome. (return negative for error,
//    return 0 repeatedly at the end of a chromosome)
// 6. Implement nextChr function, which will skip the rest of the fragments on the  
//    current chromosome (if any remain), and go to the next chromosome if it exists
// 7. chrNames and cellNames should return NULL if either the chromosome/cell ID doesn't
//    exist, (or perhaps it hasn't been read from disk yet for 10x fragments)
class FragmentsLoader {
public:
    virtual ~FragmentsLoader() = default;

    virtual bool isSeekable() const = 0;
    // Move loader to just before fragments which end after "base".
    // It's possible that fragments returned after seek will end before "base",
    // but it's guaranteed that it won't skip over any fragments ending before "base"
    virtual void seek(uint32_t chr_id, uint32_t base) = 0;
    
    // Reset the loader to start from the beginning
    virtual void restart() = 0;

    // Return the number of cells/chromosomes, or return -1 if this number is 
    // not known ahead of time
    virtual int chrCount() const = 0;
    virtual int cellCount() const = 0;

    // Return name for a given chr_id or cell_id. Only valid to call
    // for chromosme or cell_ids that have been actually returned by the loader
    virtual const char* chrNames(uint32_t chr_id) const = 0;
    virtual const char* cellNames(uint32_t cell_id) const = 0;
    
    // Advance the loader to the next chromosome. Return false if there are no more chromosomes
    virtual bool nextChr() = 0;
    // Return chromosome ID of current fragment
    virtual uint32_t currentChr() const = 0;

    // Return number of items loaded. Should repeatedly return 0 at the end of a chromosome.
    // Return -1 for error
    virtual int32_t load(uint32_t count, FragmentArray &buffer) = 0;
};

// Wrapper for a FragmentsLoader, forwarding all methods but load to the
// inner loader object. Designed to allow easier writing of
// loaders that perform transformations and filters on other loaders
class FragmentsLoaderWrapper : public FragmentsLoader {
protected:
    FragmentsLoader &loader;
public:
    FragmentsLoaderWrapper(FragmentsLoader &loader);
    bool isSeekable() const override;
    void seek(uint32_t chr_id, uint32_t base) override;
    
    void restart() override;

    int chrCount() const override;
    int cellCount() const override;

    const char* chrNames(uint32_t chr_id) const override;
    const char* cellNames(uint32_t cell_id) const override;

    bool nextChr() override;
    uint32_t currentChr() const override;

};


// Class to conveniently iterate over fragments from a FragmentsLoader
class FragmentsIterator : public FragmentsLoaderWrapper {
private:
    const uint32_t chunk_capacity;
    int32_t chunk_size;
    uint32_t idx;
    uint32_t current_chr;
    std::vector<uint32_t> start_buf, end_buf, cell_buf;
    FragmentArray fragments_buf;
public:
    // Construct iterator with a given internal buffer size (must be a power of 2 >= 128)
    FragmentsIterator(FragmentsLoader &loader, uint32_t buffer_size = 1024);

    virtual ~FragmentsIterator() = default;
    
    // Return false if there isn't a nextFragment in the current chromosome
    inline bool nextFrag() {
        idx += 1;
        if (idx >= chunk_size) {
            chunk_size = loader.load(chunk_capacity, fragments_buf);
            idx = 0;
        }
        return chunk_size > 0;
    }
    // Advance the iterator to the next chromosome. Return false if there are no more chromosomes
    inline bool nextChr() override {
        bool res = loader.nextChr();
        if (res) current_chr = loader.currentChr();
        return res;
    }
    // Access chr, start, end, cell from current fragment
    inline uint32_t chr() const {return current_chr; };
    inline uint32_t start() const {return fragments_buf.start[idx]; };
    inline uint32_t end() const {return fragments_buf.end[idx]; };
    inline uint32_t cell() const {return fragments_buf.cell[idx]; };   
    
    inline uint32_t get_chunk_capacity() {return chunk_capacity;}; 

    // Move iterator to just before fragments which end after "base".
    // It's possible that fragments returned after seek will end before "base",
    // but it's guaranteed that it won't skip over any fragments ending before "base"
    inline void seek(uint32_t chr_id, uint32_t base) override {
        loader.seek(chr_id, base);
        idx = chunk_capacity;
    }   

    int32_t load(uint32_t count, FragmentArray &buffer) override;
};

class FragmentsWriter {
public:
    // Return false on failure, true on success. 
    // During progress, call checkInterrupt which will raise an exception if there's 
    // a user-requested interrupt
    virtual bool write(FragmentsIterator &fragments, void (*checkInterrupt)(void) = NULL) = 0;
};

} // end namespace BPCells