// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <atomic>
#include <cstdlib>
#include <string>

#include "StoredMatrixSorter.h"
#include "../utils/gzfile_wrapper.h"

namespace BPCells {

class MatrixMarketHeader {
  public:
    enum class Type { Int, Real, Pattern } type;
    uint32_t rows, cols;
    uint64_t entries;
};

// Class to load data from a MatrixMarket format
// For now, skip support on symmetric and pattern matrices
template <typename T> class MatrixMarketImport : public StoredMatrixSorter<T> {
  private:
    gzFileWrapper f;
    std::vector<char> line_buf;
    int64_t remaining_entries;
    bool row_major;
    bool pattern;

    static const char *nextField(const char *c);

    static bool read_line(gzFile f, std::vector<char> &out);

    inline T parse_value(const char *str) {
        if constexpr (std::is_same_v<T, uint32_t>) {
            return std::atoi(str);
        }
        if constexpr (std::is_same_v<T, uint64_t>) {
            return std::atoll(str);
        }
        if constexpr (std::is_same_v<T, float>) {
            return std::strtof(str, NULL);
        }
        if constexpr (std::is_same_v<T, double>) {
            return std::strtod(str, NULL);
        }
    }

    inline bool parse_line(uint32_t &row, uint32_t &col, T &val) {
        if (!read_line(*f, line_buf)) return false;

        const char *cur_field, *next_field;
        // Parse the inputs
        cur_field = &line_buf[0];
        next_field = nextField(cur_field);
        if (!std::isspace(*next_field))
            throw std::runtime_error("MatrixMarketImport: expected space or tab between numbers");
        row = atoi(cur_field) - 1;

        cur_field = next_field + 1;
        next_field = nextField(cur_field);
        if (!std::isspace(*next_field))
            throw std::runtime_error("MatrixMarketImport: expected space or tab between numbers");
        col = atoi(cur_field) - 1;

        cur_field = next_field + 1;
        next_field = nextField(cur_field);
        if (pattern) {
            val = 1;
        } else {
            val = parse_value(cur_field);
        }

        return true;
    }

    static MatrixMarketHeader parse_header(std::string path);

  protected:
    size_t load_entries(
        std::vector<uint32_t> &row,
        std::vector<uint32_t> &col,
        std::vector<T> &val,
        std::atomic<bool> *user_interrupt
    ) override {
        if (row_major) std::swap(row, col);

        size_t i;
        for (i = 0; i < row.size(); i++) {
            if (i % 512 == 0 && user_interrupt != NULL && *user_interrupt) break;
            if (!parse_line(row[i], col[i], val[i])) break;
        }
        if (row_major) std::swap(row, col);
        remaining_entries -= i;
        return i;
    }
    bool row_sorted() const override { return true; }

  public:
    MatrixMarketImport(
        std::string input_path,
        WriterBuilder &output,
        const char *tmpdir,
        uint64_t load_bytes,
        uint64_t sort_buffer_bytes,
        bool row_major = false
    )
        : StoredMatrixSorter<T>(output, tmpdir, load_bytes, sort_buffer_bytes, row_major)
        , f(input_path, "r", 1<<20)
        , row_major(row_major) {

        MatrixMarketHeader h = parse_header(input_path);
        remaining_entries = h.entries;
        if (h.type == MatrixMarketHeader::Type::Real && std::is_integral_v<T>) {
            throw std::runtime_error(
                "MatrixMarketImport: Requested parsing real data type as integers"
            );
        }
        pattern = h.type == MatrixMarketHeader::Type::Pattern;
        // Read through the header lines
        while (read_line(*f, line_buf) && line_buf[0] == '%') {
            ;
        }
    }

    friend void importMtx(
        std::string input_path,
        std::vector<std::string> &&row_names,
        std::vector<std::string> &&col_names,
        WriterBuilder &output,
        const char *tmpdir,
        uint64_t load_bytes,
        uint64_t sort_buffer_bytes,
        bool row_major,
        std::atomic<bool> *user_interrupt
    );
};

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
);

} // end namespace BPCells