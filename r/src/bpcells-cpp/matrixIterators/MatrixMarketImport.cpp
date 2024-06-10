// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "MatrixMarketImport.h"

namespace BPCells {

// Import mtx file and create output file of correct type
void importMtx(
    std::string input_path,
    std::vector<std::string> &&row_names,
    std::vector<std::string> &&col_names,
    WriterBuilder &output,
    const char *tmpdir,
    uint64_t load_bytes,
    uint64_t sort_buffer_bytes,
    bool row_major,
    std::atomic<bool> *user_interrupt
) {
    MatrixMarketHeader h;
    h = MatrixMarketImport<uint32_t>::parse_header(input_path);
    if (h.rows != row_names.size() && row_names.size() != 0) throw std::runtime_error("importMtx: row_names not same length as row count");
    if (h.cols != col_names.size() && col_names.size() != 0) throw std::runtime_error("importMtx: col_names not same length as col count");
    if (row_major) {
        std::swap(h.rows, h.cols);
        std::swap(row_names, col_names);
    }
    if (h.type == MatrixMarketHeader::Type::Real) {
        MatrixMarketImport<double> importer(
            input_path, output, tmpdir, load_bytes, sort_buffer_bytes, row_major
        );
        importer.writeValues(
            std::move(row_names), std::move(col_names), h.rows, h.cols, user_interrupt
        );
        if (importer.remaining_entries != 0) {
            throw std::runtime_error("importMtx: Detected truncated mtx input");
        }
    } else {
        MatrixMarketImport<uint32_t> importer(
            input_path, output, tmpdir, load_bytes, sort_buffer_bytes, row_major
        );
        importer.writeValues(
            std::move(row_names), std::move(col_names), h.rows, h.cols, user_interrupt
        );
        if (importer.remaining_entries != 0) {
            throw std::runtime_error("importMtx: Detected truncated mtx input");
        }
    }
}

template <typename T> const char *MatrixMarketImport<T>::nextField(const char *c) {
    // Go past a set of spaces, then a set of non-spaces
    while (*c != '\0' && std::isspace(*c)) {
        c++;
    }
    while (*c != '\0' && !std::isspace(*c)) {
        c++;
    }
    return c;
}

// Read the next line into out vectr, returning false if we tried reading past the end of the file.
// Resize out as needed to make the data fit
template <typename T> bool MatrixMarketImport<T>::read_line(gzFile f, std::vector<char> &out) {
    size_t offset = 0;
    while (true) {
        while (offset + 1 >= out.size()) {
            out.resize(out.size() * 2 + 3);
        }
        out.back() = 'X';
        if (gzgets(f, out.data() + offset, out.size() - offset) == NULL) {
            if (!gzeof(f)) {
                int err;
                throw std::runtime_error(
                    "Error reading from gzfile:" + std::string(gzerror(f, &err))
                );
            }
            return false;
        }
        if (out.back() != '\0') break;
        if (out[out.size() - 1] == '\n') return true;
        offset = out.size() - 1;
    }

    return true;
}

template <typename T>
MatrixMarketHeader MatrixMarketImport<T>::parse_header(std::string path) {
    gzFileWrapper f(path, "r");
    std::vector<char> line_buf;

    if (!read_line(*f, line_buf))
        throw std::runtime_error("MatrixMarketImport: error parsing header line");
    char type_ptr[128], symmetry_ptr[128];
    if (2 !=
        sscanf(&line_buf[0], "%%%%MatrixMarket matrix coordinate %127s %127s", type_ptr, symmetry_ptr)) {
        throw std::runtime_error("MatrixMarketImport: error parsing header line");
    }
    std::string type(type_ptr), symmetry(symmetry_ptr);
    MatrixMarketHeader ret;

    if (type == "real") {
        ret.type = MatrixMarketHeader::Type::Real;
    } else if (type == "integer") {
        ret.type = MatrixMarketHeader::Type::Int;
    } else if (type == "pattern") {
        ret.type = MatrixMarketHeader::Type::Pattern;
    } else {
        throw std::runtime_error(
            "MatrixMarketImport: error parsing header line: unsupported data type: " + type
        );
    }

    if (symmetry != "general") {
        throw std::runtime_error(
            "MatrixMarketImport: error parsing header line: unsupported symmetry option: " +
            symmetry
        );
    }

    // Read the rest of the header
    while (true) {
        if (!read_line(*f, line_buf))
            throw std::runtime_error("MatrixMarketImport: error reading header");
        if (line_buf[0] != '%') break;
    }

    // Parse the dimensions
    const char *cur_field, *next_field;
    cur_field = &line_buf[0];
    next_field = nextField(cur_field);
    if (!std::isspace(*next_field))
        throw std::runtime_error(
            "MatrixMarketImport: dimensions line expected space or tab between numbers"
        );
    ret.rows = atoi(cur_field);

    cur_field = next_field + 1;
    next_field = nextField(cur_field);
    if (!std::isspace(*next_field))
        throw std::runtime_error(
            "MatrixMarketImport: dimensions line expected space or tab between numbers"
        );
    ret.cols = atoi(cur_field);

    cur_field = next_field + 1;
    next_field = nextField(cur_field);
    if (!std::isspace(*next_field))
        throw std::runtime_error(
            "MatrixMarketImport: dimensions line expected space or tab between numbers"
        );
    ret.entries = atoll(cur_field);

    return ret;
}

template class MatrixMarketImport<uint32_t>;
template class MatrixMarketImport<double>;
} // namespace BPCells