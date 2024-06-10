// Copyright 2021 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <array>
#include <atomic>
#include <string>
#include <unordered_map>

#include <zlib.h>

#include "../utils/gzfile_wrapper.h"

#include "FragmentIterator.h"

namespace BPCells {

// Read a fragment TSV with columns chr, start, end, cell_id [optional others] in that order.
// cell and chromosme IDs are assigned in sequential order from the order they're seen
class BedFragments : public FragmentLoader {
  public:
    BedFragments(const char *path, const char *comment_prefix = "");

    BedFragments() = delete;
    BedFragments(const BedFragments &) = delete;
    BedFragments &operator=(const BedFragments &other) = delete;

    // Reset the iterator to start from the beginning
    void restart() override;

    // Return the number of cells/chromosomes, or return -1 if this number is
    // not known ahead of time
    int chrCount() const override;
    int cellCount() const override;

    const char *chrNames(uint32_t chr_id) override;
    const char *cellNames(uint32_t cell_id) override;

    bool nextChr() override;
    uint32_t currentChr() const override;

    bool isSeekable() const override;
    void seek(uint32_t chr_id, uint32_t base) override;

    bool load() override;
    uint32_t capacity() const override;

    uint32_t *cellData() override;
    uint32_t *startData() override;
    uint32_t *endData() override;

  private:
    std::string path;
    gzFileWrapper f;
    std::array<char, 1 << 10> line_buf;
    std::vector<std::string> chr_names, cell_names;
    std::unordered_map<std::string, uint32_t> chr_lookup, cell_id_lookup;
    uint32_t next_chr_id, next_cell_id;
    bool eof = false;
    std::string current_chr;
    std::string comment;
    uint32_t last_start = 0;
    std::vector<uint32_t> cell, start, end;

    const char *nextField(const char *c);

    // Read the next line, returning false if we tried reading past the end of
    // the file
    bool read_line();

    // Parse the line in line_buf, returning the chromosome name as the actual
    // return value, with output parameters for start, end, cell_id.
    // Will assign a cell_id if it sees a new cell name.
    // Returns empty string at eof
    std::string_view parse_line(uint32_t &start, uint32_t &end, uint32_t &cell_id);

    bool validInt(const char *c);
};

class BedFragmentsWriter : public FragmentWriter {
  public:
    BedFragmentsWriter(
        const char *path, bool append_5th_column = false, uint32_t buffer_size = 1 << 20
    );
    void write(FragmentLoader &fragments, std::atomic<bool> *user_interrupt = NULL) override;

  private:
    gzFileWrapper f;
    bool append_5th_column;
};

} // end namespace BPCells