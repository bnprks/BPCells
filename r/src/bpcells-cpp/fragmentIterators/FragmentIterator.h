// Copyright 2021 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once
#include <atomic>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace BPCells {

// Interface for loading sorted fragments from files or memory.
// Main workhorse is the `load` function, which will load data into an internal buffer
// that users can access through `cellData()`, `startData()`, and `endData()`. There will be
// `capacity()` elements loaded starting at each of those pointers, and the loaded
// data can be modified by the caller as desired.
// Users will call load repeatedly, until it hits the end of a chromosome or the end
// of a file.
// Optionally, a FragmentLoader can support seeking, which skips the current
// reading position to a given genomic coordinate.
// Requirements for building a FragmentLoader:
// 1. Coordinates are assumed to be 0-based and half-open, same as the bed file format.
//    This means the first base in a chromosome is at position 0, and the end coordinate
//    of a fragment is one basepair after the last mapped base.
//    For more discussion, see here:
//    http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/
// 2. Assign chromosme and cell IDs that count up from zero with no skipping. It's
//    okay if a given ID is never emitted, but there must be a string ID which can
//    be returned from chrNames or cellNames
// 3. The number of valid chromosomes and cell IDs must be known in advance for
//    seekable FragmentIterators, along with their corresponding names
// 4. Fragments are grouped by chromosme, although they need not be covered in
//    order of chromosome ID.
// 5. Implement load function, which will always load >= 1 fragment, returning false
//    if there are no more fragments in the current chromosome.
// 6. Implement nextChr function, which will skip the rest of the fragments on the
//    current chromosome (if any remain), and go to the next chromosome if it exists
// 7. chrNames and cellNames should return NULL if either the chromosome/cell ID doesn't
//    exist, (or perhaps it hasn't been read from disk yet for 10x fragments)
class FragmentLoader {
  public:
    virtual ~FragmentLoader() = default;

    virtual bool isSeekable() const = 0;
    // Skip the read location to fragments which overlap a coord >= `base`, without
    // loading any actual data.
    // Note that seek does not guarantee the next loaded fragment will have an end > base,
    // it's just a performance optimization that will get close without overshooting
    virtual void seek(uint32_t chr_id, uint32_t base) = 0;

    // Reset the loader to start from the beginning
    virtual void restart() = 0;

    // Return the number of cells/chromosomes, or return -1 if this number is
    // not known ahead of time
    virtual int chrCount() const = 0;
    virtual int cellCount() const = 0;

    // Return name for a given chr_id or cell_id. Only valid to call
    // for chromosme or cell_ids that have been actually returned by the loader
    virtual const char *chrNames(uint32_t chr_id) = 0;
    virtual const char *cellNames(uint32_t cell_id) = 0;

    // Advance the loader to the next chromosome. Return false if there are no more chromosomes
    virtual bool nextChr() = 0;
    // Return chromosome ID of current fragment
    virtual uint32_t currentChr() const = 0;

    // Return false if there are no more fragments to load on the current chromosome
    virtual bool load() = 0;

    // Number of loaded fragments available in the cell, start, and end pointers
    virtual uint32_t capacity() const = 0;

    // Pointers to the loaded data
    virtual uint32_t *cellData() = 0;
    virtual uint32_t *startData() = 0;
    virtual uint32_t *endData() = 0;
};

// Wrapper for a FragmentLoader, forwarding all to the inner loader object.
// Typically, a child class will override load and/or other methods
// Designed to allow easier writing of
// loaders that perform transformations and filters on other loaders
class FragmentLoaderWrapper : public FragmentLoader {
  protected:
    std::unique_ptr<FragmentLoader> loader;
    bool take_ownership = true;
  public:
    FragmentLoaderWrapper(std::unique_ptr<FragmentLoader> &&loader);
    
    ~FragmentLoaderWrapper() {
      if (!take_ownership) loader.release();
    }

    FragmentLoaderWrapper() = default;
    FragmentLoaderWrapper(FragmentLoaderWrapper&&) = default;
    FragmentLoaderWrapper& operator=(FragmentLoaderWrapper&&) = default;

    // Set the object so that the inner loader will be preserved
    // rather than calling the destructor when this loader is destructed
    void preserve_input_loader() {take_ownership = false;}
    
    bool isSeekable() const override;
    void seek(uint32_t chr_id, uint32_t base) override;

    void restart() override;

    int chrCount() const override;
    int cellCount() const override;

    const char *chrNames(uint32_t chr_id) override;
    const char *cellNames(uint32_t cell_id) override;

    bool nextChr() override;
    uint32_t currentChr() const override;

    bool load() override;
    uint32_t capacity() const override;

    uint32_t *cellData() override;
    uint32_t *startData() override;
    uint32_t *endData() override;
};

// Class to conveniently iterate over fragments from a FragmentLoader
class FragmentIterator : public FragmentLoaderWrapper {
  private:
    uint32_t idx = UINT32_MAX;
    uint32_t current_chr;
    uint32_t current_capacity = 0;
    uint32_t *current_cell, *current_start, *current_end;

  public:
    FragmentIterator(std::unique_ptr<FragmentLoader> &&loader);

    inline void restart() override {
        loader->restart();
        idx = UINT32_MAX;
        current_capacity = 0;
    }

    // Return false if there isn't a nextFragment in the current chromosome
    inline bool nextFrag() {
        idx += 1;
        if (idx >= current_capacity) {
            return load();
        }
        return true;
    }
    // Advance the iterator to the next chromosome. Return false if there are no more chromosomes
    inline bool nextChr() override {
        bool res = loader->nextChr();
        if (res) current_chr = loader->currentChr();
        idx = UINT32_MAX;
        current_capacity = 0;
        return res;
    }
    // Access chr, start, end, cell from current fragment
    inline uint32_t chr() const { return current_chr; };
    inline uint32_t start() const { return current_start[idx]; };
    inline uint32_t end() const { return current_end[idx]; };
    inline uint32_t cell() const { return current_cell[idx]; };

    // Move iterator to just before fragments which end after "base".
    // It's possible that fragments returned after seek will end before "base",
    // but it's guaranteed that it won't skip over any fragments ending before "base"
    inline void seek(uint32_t chr_id, uint32_t base) override {
        loader->seek(chr_id, base);
        idx = UINT32_MAX;
        current_capacity = 0;
    }

    bool load() override;
    uint32_t capacity() const override;
};

class FragmentWriter {
  public:
    // Write fragments
    // During progress, quit eraly if user_interrupt becomes true
    virtual void write(FragmentLoader &fragments, std::atomic<bool> *user_interrupt = NULL) = 0;
};

} // end namespace BPCells