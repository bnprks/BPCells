// Copyright 2023 BPCells contributors
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include <vector>

#include "fragments.hpp"
#include "py_interrupts.hpp"

#include "bpcells-cpp/arrayIO/binaryfile.h"
#include "bpcells-cpp/arrayIO/vector.h"
#include "bpcells-cpp/fragmentIterators/BedFragments.h"
#include "bpcells-cpp/fragmentIterators/CellSelect.h"
#include "bpcells-cpp/fragmentIterators/ShiftCoords.h"
#include "bpcells-cpp/fragmentIterators/StoredFragments.h"
#include "bpcells-cpp/matrixIterators/ConcatenateMatrix.h"
#include "bpcells-cpp/matrixIterators/MatrixIndexSelect.h"
#include "bpcells-cpp/matrixIterators/StoredMatrixSparseColumn.h"
#include "bpcells-cpp/matrixIterators/RenameDims.h"
#include "bpcells-cpp/matrixIterators/TileMatrix.h"
#include "bpcells-cpp/utils/filesystem_compat.h"


namespace BPCells::py {

void import_10x_fragments(
    std::string input_10x,
    std::string output_bpcells,
    int shift_start,
    int shift_end,
    std::optional<std::vector<std::string>> keeper_cells
) {
    std::unique_ptr<FragmentLoader> frags = std::make_unique<BedFragments>(input_10x.c_str(), "#");
    FileWriterBuilder wb(output_bpcells);

    if (shift_start != 0 || shift_end != 0) {
        frags = std::make_unique<ShiftCoords>(std::move(frags), shift_start, shift_end);
    }
    if (keeper_cells) {
        frags = std::make_unique<CellNameSelect>(std::move(frags), keeper_cells.value());
    }

    run_with_py_interrupt_check(
        &StoredFragmentsWriter::write, StoredFragmentsWriter::createPacked(wb), std::ref(*frags)
    );
}

std::vector<std::string> cell_names_fragments_dir(std::string input_bpcells) {
    FileReaderBuilder rb(input_bpcells);
    StoredFragmentsPacked frags = StoredFragmentsPacked::openPacked(rb);
    std::vector<std::string> ret;
    for (int i = 0;; i++) {
        const char *name = frags.cellNames(i);
        if (name == NULL) break;
        ret.push_back(name);
    }
    return ret;
}

std::vector<std::string> chr_names_fragments_dir(std::string input_bpcells) {
    FileReaderBuilder rb(input_bpcells);
    StoredFragmentsPacked frags = StoredFragmentsPacked::openPacked(rb);
    std::vector<std::string> ret;
    for (int i = 0;; i++) {
        const char *name = frags.chrNames(i);
        if (name == NULL) break;
        ret.push_back(name);
    }
    return ret;
}

static Eigen::MatrixXi matrix_loader_to_eigen_helper(
    std::unique_ptr<MatrixLoader<uint32_t>> &&loader, std::atomic<bool> *user_interrupt
) {
    // Load matrix to Eigen dense while throwing out the last row
    MatrixIterator<uint32_t> it(std::move(loader));

    Eigen::MatrixXi res(it.rows() - 1, it.cols());
    res.setZero();

    while (it.nextCol()) {
        const uint32_t col = it.currentCol();
        if (user_interrupt != NULL && *user_interrupt) return res;
        while (it.nextValue()) {
            if (it.row() < res.rows()) {
                res(it.row(), col) += it.val();
            }
        }
    }
    return res;
}

Eigen::MatrixXi pseudobulk_coverage(
    std::string fragments_path,
    std::vector<std::string> chr,
    std::vector<uint32_t> start,
    std::vector<uint32_t> end,
    std::vector<int32_t> cell_groups,
    int bin_size
) {
    FileReaderBuilder rb(fragments_path);
    std::unique_ptr<FragmentLoader> frags =
        std::make_unique<StoredFragmentsPacked>(StoredFragmentsPacked::openPacked(rb));

    // Merge cells from the same groups
    uint32_t num_groups = 1 + *std::max_element(cell_groups.cbegin(), cell_groups.cend());

    std::vector<std::string> dummy_names;
    for (uint32_t i = 0; i < num_groups + 1; i++) {
        dummy_names.push_back(std::to_string(i));
    }

    std::vector<uint32_t> cell_groups_uint;
    for (const auto &x : cell_groups) {
        cell_groups_uint.push_back(x >= 0 ? x : num_groups + 1);
    }

    frags = std::make_unique<CellMerge>(
        std::move(frags), cell_groups_uint, std::make_unique<VecStringReader>(dummy_names)
    );

    // Construct the tile matrix
    if (chr.size() != start.size() || chr.size() != end.size()) {
        throw std::runtime_error("pseudbulk_coverage: chr, start, and end must be matching lengths"
        );
    }

    std::vector<uint32_t> width;
    for (size_t i = 0; i < start.size(); i++) {
        width.push_back(bin_size);
    }

    for (size_t i = 0; i < start.size(); i++) {
        if (end[i] - start[i] != end[0] - start[0]) {
            throw std::runtime_error(
                "pseudobulk_coverage: start - end must be identical for all input regions"
            );
        }
    }

    std::vector<uint32_t> chr_id;
    std::unordered_map<std::string, uint32_t> chr_name_lookup;
    std::vector<std::string> chr_levels;
    for (int32_t i = 0; i < frags->chrCount(); i++) {
        if (frags->chrNames(i) == NULL) {
            throw std::runtime_error("pseudobulk_coverage: missing chr names in input fragments");
        }
        chr_name_lookup[std::string(frags->chrNames(i))] = i;
        chr_levels.push_back(std::string(frags->chrNames(i)));
    }
    for (auto &c : chr) {
        chr_id.push_back(chr_name_lookup[c]);
    }

    std::unique_ptr<MatrixLoader<uint32_t>> tile_mat = std::make_unique<TileMatrix>(
        std::move(frags), chr_id, start, end, width, std::make_unique<VecStringReader>(chr_levels)
    );

    return run_with_py_interrupt_check(&matrix_loader_to_eigen_helper, std::move(tile_mat));
}

template <typename T>
void parallel_map_helper(std::vector<std::future<T>> &futures, size_t threads) {
    // Non-threaded fallback
    if (threads == 0) {
        for (size_t i = 0; i < futures.size(); i++) {
            futures[i].get();
        }
        return;
    }

    // Very basic threading, designed for small numbers of futures
    std::atomic<size_t> task_id(0);
    std::atomic<bool> has_error = false;
    std::exception_ptr exception;
    std::vector<std::thread> thread_vec;
    for (size_t i = 0; i < threads; i++) {
        thread_vec.push_back(std::thread([&futures, &task_id, &has_error, &exception] {
            while (true) {
                size_t cur_task = task_id.fetch_add(1);
                if (cur_task >= futures.size()) break;
                try {
                    futures[cur_task].get();
                } catch (...) {
                    if (!has_error) {
                        has_error = true;
                        exception = std::current_exception();
                    }
                    break;
                }
            }
        }));
    }
    for (auto &th : thread_vec) {
        if (th.joinable()) {
            th.join();
        }
    }
    if (has_error) {
        std::rethrow_exception(exception);
    }
}

// Write a chunk of the tile matrix columns to the given output path
static void precalculate_pseudobulk_coverage_helper(
    std::string fragments_path,
    std::string chunk_output_path,
    std::pair<uint32_t, uint32_t> chunk_col_range,

    const std::vector<uint32_t> &group_ids,
    const std::vector<std::string> &group_names,

    const std::vector<uint32_t> &chr_id,
    const std::vector<uint32_t> &start,
    const std::vector<uint32_t> &end,
    const std::vector<uint32_t> &width,
    const std::vector<std::string> &chr_levels,

    std::atomic<bool> *user_interrupt
) {
    FileReaderBuilder rb(fragments_path);
    std::unique_ptr<FragmentLoader> frags =
        std::make_unique<StoredFragmentsPacked>(StoredFragmentsPacked::openPacked(rb));

    // Merge cells
    frags = std::make_unique<CellMerge>(
        std::move(frags), group_ids, std::make_unique<VecStringReader>(group_names)
    );

    // Construct tile matrix
    std::unique_ptr<MatrixLoader<uint32_t>> tile_mat = std::make_unique<TileMatrix>(
        std::move(frags), chr_id, start, end, width, std::make_unique<VecStringReader>(chr_levels)
    );

    // Subset to the desired columns
    tile_mat = std::make_unique<MatrixColSlice<uint32_t>>(
        std::move(tile_mat), chunk_col_range.first, chunk_col_range.second
    );

    // Clear the row/col names
    std::vector<std::string> empty;
    tile_mat = std::make_unique<RenameDims<uint32_t>>(
        std::move(tile_mat), empty, empty, true, true
    );

    // Write to output
    FileWriterBuilder wb(chunk_output_path);
    EXPERIMENTAL_createPackedSparseColumn<uint32_t>(wb).write(*tile_mat, user_interrupt);
}

void precalculate_pseudobulk_coverage(
    std::string fragments_path,
    std::string output_path,
    std::string tmp_path,
    std::vector<std::string> chr,
    std::vector<uint32_t> chr_len,
    std::vector<int32_t> cell_groups,
    int bin_size,
    int threads
) {
    // Create the arguments needed for TileMatrix: start, tile_width, chr_id, chr_levels (end =
    // chr_len)
    if (chr.size() != chr_len.size()) {
        throw std::runtime_error(
            "precalculate_pseudobulk_coverage: chr must be same length as cell_len"
        );
    }

    std::vector<uint32_t> start(chr.size(), 0);
    std::vector<uint32_t> tile_width(chr.size(), bin_size);

    FileReaderBuilder rb(fragments_path);
    auto frags = StoredFragmentsPacked::openPacked(rb);

    std::vector<uint32_t> chr_id;
    std::unordered_map<std::string, uint32_t> chr_name_lookup;
    std::vector<std::string> chr_levels;
    for (int32_t i = 0; i < frags.chrCount(); i++) {
        if (frags.chrNames(i) == NULL) {
            throw std::runtime_error("pseudobulk_coverage: missing chr names in input fragments");
        }
        chr_name_lookup[std::string(frags.chrNames(i))] = i;
        chr_levels.push_back(std::string(frags.chrNames(i)));
    }
    for (auto &c : chr) {
        chr_id.push_back(chr_name_lookup[c]);
    }

    // Create the arguments needed for MergeCellls (group_id, group_names)
    uint32_t num_groups = 1 + *std::max_element(cell_groups.cbegin(), cell_groups.cend());

    std::vector<std::string> dummy_names;
    for (uint32_t i = 0; i < num_groups + 1; i++) {
        dummy_names.push_back(std::to_string(i));
    }

    std::vector<uint32_t> cell_groups_uint;
    for (const auto &x : cell_groups) {
        cell_groups_uint.push_back(x >= 0 ? x : num_groups + 1);
    }

    // Split columns into chunks
    size_t total_columns = 0;
    for (const auto &x : chr_len) {
        total_columns += x;
    }
    uint32_t chunks = std::max<uint32_t>(1, threads * 4);
    std::vector<std::pair<uint32_t,uint32_t>> chunk_col_splits;
    uint32_t idx = 0;
    for (uint32_t i = 0; i < chunks; i++) {
        uint32_t col_count = (total_columns - idx) / (chunks - i);
        chunk_col_splits.push_back({idx, idx + col_count});
        idx += col_count;
    }

    std::vector<std::string> chunk_output_paths;
    for (uint32_t i = 0; i < chunks; i++) {
        chunk_output_paths.push_back((std_fs::path(tmp_path) / std::to_string(i)).string());
    }

    // Make all the matrix chunks
    run_with_py_interrupt_check([&fragments_path,
                                        &chunk_output_paths,
                                        &chunk_col_splits,
                                        &cell_groups_uint,
                                        &dummy_names,
                                        &chr_id,
                                        &start,
                                        &chr_len,
                                        &tile_width,
                                        &chr_levels,
                                        threads,
                                        chunks](std::atomic<bool> *user_interrupt) {
        std::vector<std::future<void>> task_vec;
        for (size_t i = 0; i < chunks; i++) {
            task_vec.push_back(std::async(
                std::launch::deferred,
                &precalculate_pseudobulk_coverage_helper,
                fragments_path,
                chunk_output_paths[i],
                chunk_col_splits[i],

                std::cref(cell_groups_uint),
                std::cref(dummy_names),

                std::cref(chr_id),
                std::cref(start),
                std::cref(chr_len),
                std::cref(tile_width),
                std::cref(chr_levels),

                user_interrupt
            ));
        }

        parallel_map_helper(task_vec, threads);
    });

    std::vector<std::unique_ptr<MatrixLoader<uint32_t>>> matrix_chunks;
    for (size_t i = 0; i < chunks; i++) {
        FileReaderBuilder rb(chunk_output_paths[i]);
        matrix_chunks.push_back(std::make_unique<StoredMatrix<uint32_t>>(EXPERIMENTAL_openPackedSparseColumn<uint32_t>(rb)));
    }

    std::unique_ptr<MatrixLoader<uint32_t>> full_mat;
    if (chunks > 1) {
        full_mat = std::make_unique<ConcatCols<uint32_t>>(std::move(matrix_chunks), 0);
    } else {
        full_mat = std::move(matrix_chunks[0]);
    }
    
    // Unselect the last row, as thos are the cells we meant to discard
    std::vector<uint32_t> row_selection;
    for (uint32_t i = 0; i < num_groups; i++) {
        row_selection.push_back(i);
    }
    full_mat = std::make_unique<MatrixRowSelect<uint32_t>>(std::move(full_mat), row_selection);
    
    FileWriterBuilder wb(output_path);
    
    run_with_py_interrupt_check(
        &StoredMatrixWriter<uint32_t>::write,
        EXPERIMENTAL_createPackedSparseColumn<uint32_t>(wb),
        std::ref(*full_mat)
    );

    for (const auto &x : chunk_output_paths) {
        std_fs::remove_all(std_fs::path(x));
    }
}

Eigen::MatrixXi query_precalculated_pseudobulk_coverage(
    std::string mat_path,
    std::vector<uint32_t> range_starts,
    uint32_t range_len
) {
    
    FileReaderBuilder rb(mat_path);
    MatrixIterator<uint32_t> mat(std::make_unique<StoredMatrix<uint32_t>>(EXPERIMENTAL_openPackedSparseColumn<uint32_t>(rb)));
    
    Eigen::MatrixXi ret(range_starts.size() * range_len, mat.rows());
    ret.setZero();

    for (const auto &x : range_starts) {
        if (x + range_len > mat.cols()) {
            throw std::runtime_error("query_precalculated_pseudobulk_coverage: A query range exceeds number of columns");
        }
    }

    for (uint32_t i = 0; i < range_starts.size(); i++) {
        mat.seekCol(range_starts[i]);
        do {
            uint32_t col = i * range_len + (mat.col() - range_starts[i]);
            while (mat.nextValue()) {
                // if (col >= ret.rows() || mat.row() >= ret.cols()) {
                //     printf("i=%d\trange_len=%d\tmat.col()=%d\trange_starts[i]=%d\nret.rows()=%d\tret.cols()=%d\n",
                //         i, range_len, mat.col(), range_starts[i], ret.rows(), ret.cols()
                //     );
                //     throw std::runtime_error("Tried to calculate bad index");
                // }
                ret(col, mat.row()) = mat.val();
            }
        } while (mat.nextCol() && mat.col() < range_starts[i] + range_len);
    }
    return ret;
}

} // namespace BPCells::py