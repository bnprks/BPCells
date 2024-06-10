// Copyright 2023 BPCells contributors
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "matrix.hpp"
#include "py_interrupts.hpp"

#include <algorithm>
#include <future>
#include <memory>
#include <optional>
#include <string>
#include <thread>
#include <utility>

#include <Eigen/SparseCore>

#include "bpcells-cpp/arrayIO/binaryfile.h"
#include "bpcells-cpp/arrayIO/vector.h"

#include "bpcells-cpp/matrixIterators/CSparseMatrix.h"
#include "bpcells-cpp/matrixIterators/MatrixIndexSelect.h"
#include "bpcells-cpp/matrixIterators/StoredMatrix.h"
#include "bpcells-cpp/matrixIterators/StoredMatrixWriter.h"
#include "bpcells-cpp/matrixIterators/ImportMatrixHDF5.h"
#include "bpcells-cpp/matrixIterators/ConcatenateMatrix.h"


namespace BPCells::py {

void write_matrix_dir_from_memory(
    const Eigen::SparseMatrix<uint32_t> in, std::string out_path, bool row_major
) {
    const Eigen::Map<Eigen::SparseMatrix<uint32_t>> in_map(
        in.rows(),
        in.cols(),
        in.nonZeros(),
        (int *)in.outerIndexPtr(),
        (int *)in.innerIndexPtr(),
        (uint32_t *)in.valuePtr()
    );

    auto mat = std::make_unique<CSparseMatrix<uint32_t>>(in_map);

    FileWriterBuilder wb(out_path);

    run_with_py_interrupt_check(
        &StoredMatrixWriter<uint32_t>::write,
        StoredMatrixWriter<uint32_t>::createPacked(wb, row_major),
        std::ref(*mat)
    );
}

static bool is_row_major_matrix_dir(std::string path) {
    FileReaderBuilder rb(path);
    auto storage_order_reader = rb.openStringReader("storage_order");
    auto storage_order = storage_order_reader->get(0);
    bool row_major = false;
    if (std::string_view("row") == storage_order) row_major = true;
    else if (std::string("col") == storage_order) row_major = false;
    else
        throw std::runtime_error(
            std::string("storage_order must be either \"row\" or \"col\", found: \"") +
            storage_order + "\""
        );
    return row_major;
}

void write_matrix_dir_from_concat(std::vector<std::string> in_paths, std::string out_path, bool concat_cols)  {
    std::vector<std::unique_ptr<MatrixLoader<uint32_t>>> mats;

    if (in_paths.size() == 0) {
        throw std::runtime_error("write_matrix_dir_from_hstack: Zero matrices given as input.");
    }

    bool row_major = is_row_major_matrix_dir(in_paths[0]);
    
    for (const std::string &path : in_paths) {
        FileReaderBuilder rb(path);
        mats.push_back(std::make_unique<StoredMatrix<uint32_t>>(StoredMatrix<uint32_t>::openPacked(rb)));
    }

    std::unique_ptr<MatrixLoader<uint32_t>> mat;
    if (row_major == concat_cols) {
        mat = std::make_unique<ConcatRows<uint32_t>>(std::move(mats), 0);
    } else {
        mat = std::make_unique<ConcatCols<uint32_t>>(std::move(mats), 0);
    }
    
    FileWriterBuilder wb(out_path);

    run_with_py_interrupt_check(
        &StoredMatrixWriter<uint32_t>::write,
        StoredMatrixWriter<uint32_t>::createPacked(wb, row_major),
        std::ref(*mat)
    );
}

void write_matrix_dir_from_h5ad(std::string h5ad_path, std::string out_path, std::string group) {
    std::vector<std::string> empty_names;
    auto row_names = std::make_unique<VecStringReader>(empty_names);
    auto col_names = std::make_unique<VecStringReader>(empty_names);

    auto mat_float = std::make_unique<StoredMatrix<float>>(openAnnDataMatrix(
        h5ad_path,
        group,
        16384L, // Just provide a default buffer size matching what R uses
        std::move(row_names),
        std::move(col_names)
    ));
    auto mat_int = std::make_unique<MatrixConverterLoader<float, uint32_t>>(std::move(mat_float));
    
    bool row_major = isRowOrientedAnnDataMatrix(h5ad_path, group);
    
    FileWriterBuilder wb(out_path);

    run_with_py_interrupt_check(
        &StoredMatrixWriter<uint32_t>::write,
        StoredMatrixWriter<uint32_t>::createPacked(wb, row_major),
        std::ref(*mat_int)
    );
}

template <typename T>
std::vector<T> parallel_map_helper(std::vector<std::future<T>> &futures, size_t threads) {
    std::vector<T> result(futures.size());

    // Non-threaded fallback
    if (threads == 0) {
        for (size_t i = 0; i < futures.size(); i++) {
            result[i] = futures[i].get();
        }
        return result;
    }

    // Very basic threading, designed for small numbers of futures
    std::atomic<size_t> task_id(0);
    std::vector<std::thread> thread_vec;
    for (size_t i = 0; i < threads; i++) {
        thread_vec.push_back(std::thread([&futures, &result, &task_id] {
            while (true) {
                size_t cur_task = task_id.fetch_add(1);
                if (cur_task >= futures.size()) break;
                result[cur_task] = futures[cur_task].get();
            }
        }));
    }
    for (auto &th : thread_vec) {
        if (th.joinable()) {
            th.join();
        }
    }
    return result;
}

VecReaderWriterBuilder load_matrix_dir_to_memory(std::string matrix_path) {
    FileReaderBuilder rb(matrix_path);
    std::unique_ptr<MatrixLoader<uint32_t>> mat =
        std::make_unique<StoredMatrix<uint32_t>>(StoredMatrix<uint32_t>::openPacked(rb));

    VecReaderWriterBuilder wb;
    run_with_py_interrupt_check(
        &StoredMatrixWriter<uint32_t>::write,
        StoredMatrixWriter<uint32_t>::createPacked(wb),
        std::ref(*mat)
    );
    return wb;
}

Eigen::SparseMatrix<uint32_t> load_matrix_subset_helper(
    ReaderBuilder &rb,
    std::optional<std::vector<uint32_t>> rows,
    std::vector<uint32_t> columns,
    std::atomic<bool> *user_interrupt
) {
    std::unique_ptr<MatrixLoader<uint32_t>> mat =
        std::make_unique<StoredMatrix<uint32_t>>(StoredMatrix<uint32_t>::openPacked(rb));

    mat = std::make_unique<MatrixColSelect<uint32_t>>(std::move(mat), columns);
    if (rows) {
        mat = std::make_unique<MatrixRowSelect<uint32_t>>(std::move(mat), rows.value());
    }

    CSparseMatrixWriter<uint32_t> writer;
    run_with_py_interrupt_check(
        &CSparseMatrixWriter<uint32_t>::write, std::ref(writer), std::ref(*mat)
    );
    return writer.getMat();
}

std::vector<Eigen::SparseMatrix<uint32_t>> load_matrix_subset(
    ReaderBuilder &rb,
    std::optional<std::vector<uint32_t>> rows,
    std::vector<uint32_t> columns,
    uint32_t threads
) {
    // Split columns into chunks
    std::vector<std::vector<uint32_t>> col_splits;
    uint32_t chunks = std::max<uint32_t>(1, threads);
    uint32_t idx = 0;
    for (uint32_t i = 0; i < chunks; i++) {
        std::vector<uint32_t> c;
        uint32_t col_count = (columns.size() - idx) / (chunks - i);
        for (uint32_t j = 0; j < col_count; j++) {
            c.push_back(columns[idx]);
            idx++;
        }
        col_splits.push_back(c);
    }

    // Perform multi-threaded reads
    return run_with_py_interrupt_check(
        [&rb, &rows, &col_splits, threads](std::atomic<bool> *user_interrupt) {
            std::vector<std::future<Eigen::SparseMatrix<uint32_t>>> task_vec;
            for (size_t i = 0; i < col_splits.size(); i++) {
                task_vec.push_back(std::async(
                    std::launch::deferred,
                    &load_matrix_subset_helper,
                    std::ref(rb),
                    rows,
                    col_splits[i],
                    user_interrupt
                ));
            }

            return parallel_map_helper(task_vec, threads);
        }
    );
}

std::vector<Eigen::SparseMatrix<uint32_t>> load_matrix_dir_subset(
    std::string matrix_path,
    std::optional<std::vector<uint32_t>> rows,
    std::vector<uint32_t> columns,
    uint32_t threads
) {
    FileReaderBuilder rb(matrix_path);
    return load_matrix_subset(rb, rows, columns, threads);
}

std::vector<Eigen::SparseMatrix<uint32_t>> load_matrix_memory_subset(
    VecReaderWriterBuilder &rb,
    std::optional<std::vector<uint32_t>> rows,
    std::vector<uint32_t> columns,
    uint32_t threads
) {
    return load_matrix_subset(rb, rows, columns, threads);
}

std::tuple<uint32_t, uint32_t> dims_matrix_dir(std::string matrix_path) {
    FileReaderBuilder rb(matrix_path);
    std::unique_ptr<MatrixLoader<uint32_t>> mat =
        std::make_unique<StoredMatrix<uint32_t>>(StoredMatrix<uint32_t>::openPacked(rb));

    return {mat->rows(), mat->cols()};
}

} // namespace BPCells::py