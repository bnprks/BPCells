// Copyright 2021 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include <atomic>
#include <charconv>
#include "BedFragments.h"
#include "../utils/filesystem_compat.h"

namespace BPCells {

BedFragments::BedFragments(const char *path, const char *comment_prefix)
    : path(path)
    , comment(comment_prefix) {
    restart();
}

// Return the number of cells/chromosomes, or return -1 if this number is
// not known ahead of time
int BedFragments::chrCount() const { return -1; }
int BedFragments::cellCount() const { return -1; }

const char *BedFragments::chrNames(uint32_t chr_id) {
    if (chr_id >= chr_names.size()) return NULL;
    return chr_names[chr_id].c_str();
}

const char *BedFragments::cellNames(uint32_t cell_id) {
    if (cell_id >= cell_names.size()) return NULL;
    return cell_names[cell_id].c_str();
}

uint32_t BedFragments::currentChr() const { return chr_lookup.at(current_chr); }

bool BedFragments::isSeekable() const { return false; }

void BedFragments::seek(uint32_t chr_id, uint32_t base) {
    throw std::invalid_argument("Cannot seek BedFragments");
}

const char *BedFragments::nextField(const char *c) {
    while (*c != '\0' && *c != '\t' && *c != '\n') {
        c++;
    }
    return c;
}

// Read the next line, returning false if we tried reading past the end of
// the file
bool BedFragments::read_line() {
    if (gzgets(*f, &line_buf[0], line_buf.size()) == NULL) {
        if (eof) {
            line_buf[0] = '\0';
            return false;
        } else if (!gzeof(*f)) {
            throw std::runtime_error("Error reading from gzfile");
        }
        eof = true;
    }

    return true;
}

// Parse the line in line_buf, returning the chromosome name as the actual
// return value, with output parameters for start, end, cell_id.
// Will assign a cell_id if it sees a new cell name.
// Returns empty string at eof
std::string_view BedFragments::parse_line(uint32_t &start, uint32_t &end, uint32_t &cell_id) {
    const char *cur_field, *next_field;

    cur_field = &line_buf[0];
    next_field = nextField(cur_field);
    std::string_view chr(cur_field, next_field - cur_field);

    if (next_field == cur_field) return chr;
    if (*next_field != '\t') throw std::runtime_error("Invalid TSV file");

    next_field = parseInt(next_field + 1, start);
    next_field = parseInt(next_field + 1, end);

    cur_field = next_field + 1;
    next_field = nextField(cur_field);
    cell_id = lookupCellID(std::string_view(cur_field, next_field - cur_field));

    return chr;
}

uint32_t BedFragments::lookupCellID(std::string_view name) {
    size_t mask = cell_id_slots.size() - 1;
    size_t slot = std::hash<std::string_view>()(name) & mask;
    for (; cell_id_slots[slot] != 0; slot = (slot + 1) & mask) {
        if (cell_names[cell_id_slots[slot] - 1] == name) return cell_id_slots[slot] - 1;
    }
    cell_names.emplace_back(name);
    cell_id_slots[slot] = cell_names.size();
    // Keep the table at most half full by doubling its size and re-inserting all names
    if (2 * cell_names.size() > cell_id_slots.size()) {
        cell_id_slots.assign(2 * cell_id_slots.size(), 0);
        mask = cell_id_slots.size() - 1;
        for (uint32_t i = 0; i < cell_names.size(); i++) {
            slot = std::hash<std::string_view>()(cell_names[i]) & mask;
            while (cell_id_slots[slot] != 0) slot = (slot + 1) & mask;
            cell_id_slots[slot] = i + 1;
        }
    }
    return cell_names.size() - 1;
}

const char *BedFragments::parseInt(const char *c, uint32_t &value) {
    auto [field_end, error] = std::from_chars(c, line_buf.data() + line_buf.size(), value);
    if (error != std::errc() || *field_end != '\t') throw std::runtime_error("Invalid TSV file");
    return field_end;
}

void BedFragments::restart() {
    // This will automatically handle closing the old copy of f if needed
    f = gzFileWrapper(path.c_str(), "r");


    // Reset the instance variables
    last_start = 0;
    current_chr = "";
    chr_lookup.clear();
    chr_names.clear();
    next_chr_id = 0;

    cell_id_slots.assign(1024, 0);
    cell_names.clear();

    cell.resize(0);
    start.resize(0);
    end.resize(0);
    eof = false;

    read_line();
    if (comment.size() == 0) return;

    // Loop through comment lines
    while (true) {
        if (line_buf[0] == '\0') break;
        bool matches_comment = true;
        for (uint32_t i = 0; i < comment.size(); i++) {
            if (line_buf[i] != comment[i]) {
                matches_comment = false;
                break;
            }
        }
        if (!matches_comment) break;
        read_line();
    }
}

bool BedFragments::nextChr() {
    if (line_buf[0] == '\0' || line_buf[0] == '\n') return false;

    // Keep reading fragments until we get to the next chromosome
    uint32_t dummy_start, dummy_end, dummy_cell;
    while (true) {
        std::string_view chr = parse_line(dummy_start, dummy_end, dummy_cell);
        if (chr == "" || chr != current_chr) {
            current_chr = chr;
            break;
        }
        if (dummy_start < last_start)
            throw std::runtime_error("TSV not in sorted order by chr, start");
        last_start = dummy_start;
        if (!read_line()) return false;
    }

    auto chr_id_res = chr_lookup.emplace(current_chr, next_chr_id);
    if (chr_id_res.second) {
        chr_names.push_back(current_chr);
        next_chr_id++;
    } else {
        throw std::runtime_error("TSV not in sorted order by chr, start");
    }
    last_start = 0;
    return true;
}

bool BedFragments::load() {
    std::string_view chr;
    uint32_t i;

    cell.resize(1024);
    start.resize(1024);
    end.resize(1024);
    for (i = 0; i < 1024; i++) {
        // line_buf will contain the next line in file before start of loop
        chr = parse_line(start[i], end[i], cell[i]);
        if (chr == "" || chr != current_chr) {
            break;
        }
        if (start[i] < last_start)
            throw std::runtime_error("TSV not in sorted order by chr, start");
        last_start = start[i];

        if (!read_line()) break;
    }
    cell.resize(i);
    start.resize(i);
    end.resize(i);
    return i;
}

uint32_t BedFragments::capacity() const { return cell.size(); }

uint32_t *BedFragments::cellData() { return cell.data(); }
uint32_t *BedFragments::startData() { return start.data(); }
uint32_t *BedFragments::endData() { return end.data(); }

BedFragmentsWriter::BedFragmentsWriter(
    const char *path, bool append_5th_column, uint32_t buffer_size
)
    : f(std::string(path), "w", buffer_size)
    , append_5th_column(append_5th_column) {}

void BedFragmentsWriter::write(FragmentLoader &loader, std::atomic<bool> *user_interrupt) {
    FragmentIterator fragments((std::unique_ptr<FragmentLoader>(&loader)));
    // Don't take ownership of the loader object
    fragments.preserve_input_loader();
    uint32_t bytes_written;

    size_t total_fragments = 0;

    const char *output_format;
    if (append_5th_column) {
        output_format = "%s\t%d\t%d\t%s\t0\n";
    } else {
        output_format = "%s\t%d\t%d\t%s\n";
    }

    fragments.restart();
    while (fragments.nextChr()) {
        const char *chr_name = fragments.chrNames(fragments.currentChr());
        while (fragments.nextFrag()) {
            bytes_written = gzprintf(
                *f,
                output_format,
                chr_name,
                fragments.start(),
                fragments.end(),
                fragments.cellNames(fragments.cell())
            );

            if (bytes_written <= 0) {
                throw std::runtime_error("Failed to write data in BedFragmentsWriter");
            }
            total_fragments += 1;
            if (total_fragments++ % 1024 == 0 && user_interrupt != NULL && *user_interrupt) return;
        }
    }
    // Force f to close
    f = gzFileWrapper();
}

} // end namespace BPCells
