#pragma once

#include <array>
#include <string>
#include <unordered_map>

#include <zlib.h>

#include "FragmentsIterator.h"

namespace BPCells {

// Read a fragment TSV with columns chr, start, end, cell_id [optional others] in that order.
// cell and chromosme IDs are assigned in sequential order from the order they're seen
class BedFragments : public FragmentsLoader {
public:
    BedFragments(const char *path, const char *comment_prefix = "");

    ~BedFragments();
    
    BedFragments() = delete;
    BedFragments(const BedFragments&) = delete;
    BedFragments& operator=(const BedFragments& other) = delete;

    // Reset the iterator to start from the beginning
    void restart() override;

    // Return the number of cells/chromosomes, or return -1 if this number is 
    // not known ahead of time
    int chrCount() const override;
    int cellCount() const override;

    const char* chrNames(uint32_t chr_id) const override;
    const char* cellNames(uint32_t cell_id) const override;
    
    bool nextChr() override;
    uint32_t currentChr() const override;

    bool isSeekable() const override;
    void seek(uint32_t chr_id, uint32_t base) override;

private:
    std::string(path);
    gzFile f;
    std::array<char, 1<<20 > line_buf;
    std::vector<std::string> chr_names, cell_names;
    std::unordered_map<std::string, uint32_t> chr_lookup, cell_id_lookup;
    uint32_t next_chr_id, next_cell_id;
    bool eof = false;
    std::string current_chr;
    std::string comment;
    uint32_t last_start;

    const char* nextField(const char * c) ;

    // Read the next line, returning false if we tried reading past the end of
    // the file
    bool read_line();

    // Parse the line in line_buf, returning the chromosome name as the actual
    // return value, with output parameters for start, end, cell_id.
    // Will assign a cell_id if it sees a new cell name.
    // Returns empty string at eof
    std::string parse_line(uint32_t &start, uint32_t &end, uint32_t &cell_id);

    bool validInt(const char* c);

    int32_t load(uint32_t count, FragmentArray buffer) override;
};


class BedFragmentsWriter : public FragmentsWriter {
public:
    BedFragmentsWriter(const char *path, bool append_5th_column=false,
                    uint32_t buffer_size = 1 << 20);
    ~BedFragmentsWriter();
    bool write(FragmentsIterator &fragments, void (*checkInterrupt)(void) = NULL) override;
private:
    gzFile f;
    bool append_5th_column;
};

} // end namespace BPCells