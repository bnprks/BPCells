// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once
#include <atomic>
#include <future>
#include <vector>

#include "MatrixIterator.h"

namespace BPCells {

namespace {

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

const Eigen::Map<Eigen::MatrixXd>
subset_map_cols(const Eigen::Map<Eigen::MatrixXd> mat, int start_col, int n_cols) {
    return Eigen::Map<Eigen::MatrixXd>((double *)&mat(0, start_col), mat.rows(), n_cols);
}

const Eigen::Map<Eigen::VectorXd>
subset_map_vec(const Eigen::Map<Eigen::VectorXd> vec, int start, int count) {
    return Eigen::Map<Eigen::VectorXd>((double *)&vec(start), count, 1);
}

} // namespace

// Concatenate a list of MatrixLoaders by row.
// Column names will be taken from the first matrix in the list
template <typename T> class ConcatRows : public MatrixLoader<T> {
  protected:
    std::vector<std::unique_ptr<MatrixLoader<T>>> mats;
    std::vector<uint32_t> row_offset;
    uint32_t cur_mat = 0;
    uint32_t threads = 0;

  public:
    ConcatRows(std::vector<std::unique_ptr<MatrixLoader<T>>> &&mats, uint32_t threads) : mats(std::move(mats)), threads(threads) {

        if (this->mats.size() < 2) throw std::runtime_error("Must have >= 2 matrices to merge");

        row_offset.push_back(0);
        uint32_t cols = this->mats.front()->cols();
        for (const auto &m : this->mats) {
            row_offset.push_back(row_offset.back() + m->rows());
            if (m->cols() != cols)
                throw std::runtime_error("ConcatRows: Matrices must have equal numbers of columns");
        }

        this->threads = std::min(this->threads, (uint32_t) this->mats.size());
    }

    uint32_t rows() const override { return row_offset.back(); }
    uint32_t cols() const override { return mats.front()->cols(); }

    const char *rowNames(uint32_t row) override {
        auto it = std::upper_bound(row_offset.begin(), row_offset.end(), row);
        uint32_t idx = it - row_offset.begin() - 1;

        if (idx == mats.size()) return NULL;

        return mats[idx]->rowNames(row - row_offset[idx]);
    }
    const char *colNames(uint32_t col) override { return mats.front()->colNames(col); }

    // Reset the iterator to start from the beginning
    void restart() override {
        cur_mat = 0;
        for (auto &m : mats)
            m->restart();
    }

    // Seek to a specific column without reading data
    // Next call should be to load(). col must be < cols()
    void seekCol(uint32_t col) override {
        cur_mat = 0;
        for (auto &m : mats)
            m->seekCol(col);
    }

    // Advance to the next column, or return false if there
    // are no more columns
    bool nextCol() override {
        cur_mat = 0;
        bool any_fail = false;
        bool all_fail = true;
        for (auto &m : mats) {
            if (m->nextCol()) all_fail = false;
            else any_fail = true;
        }
        if (all_fail) return false;
        if (any_fail)
            throw std::runtime_error(
                "ConcatRows: Some matrices reached nextCol while others did not"
            );
        return true;
    }

    // Return the index of the current column
    uint32_t currentCol() const override { return mats.front()->currentCol(); }

    // Return false if there are no more entries to load
    bool load() override {
        while (!mats[cur_mat]->load()) {
            if (cur_mat + 1 == mats.size()) return false;
            cur_mat++;
        }

        uint32_t *row_data = mats[cur_mat]->rowData();
        uint32_t cap = mats[cur_mat]->capacity();
        for (uint32_t i = 0; i < cap; i++)
            row_data[i] += row_offset[cur_mat];
        return true;
    }

    // Number of loaded entries available
    uint32_t capacity() const override { return mats[cur_mat]->capacity(); }

    // Pointers to the loaded entries
    uint32_t *rowData() override { return mats[cur_mat]->rowData(); }
    T *valData() override { return mats[cur_mat]->valData(); }

    Eigen::MatrixXd
    denseMultiplyRight(const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt) override {
        if (cols() != B.rows()) throw std::runtime_error("Incompatible dimensions for matrix multiply");
        std::vector<std::future<Eigen::MatrixXd>> task_vec;
        // Multiply separately, and concatenate chunks
        for (size_t i = 0; i < mats.size(); i++) {
            task_vec.push_back(std::async(
                std::launch::deferred,
                &MatrixLoader<T>::denseMultiplyRight,
                mats[i].get(),
                B,
                user_interrupt
            ));
        }
        std::vector<Eigen::MatrixXd> sub_results = parallel_map_helper(task_vec, threads);
        Eigen::MatrixXd res(B.cols(), rows());
        if (user_interrupt != NULL && *user_interrupt) return res.transpose();
        for (size_t i = 0; i < mats.size(); i++) {
            res.middleCols(row_offset[i], mats[i]->rows()) = sub_results[i].transpose();
        }
        return res.transpose();
    }

    Eigen::VectorXd
    vecMultiplyRight(const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt) override {
        if (cols() != v.rows()) throw std::runtime_error("Incompatible dimensions for vector multiply");
        std::vector<std::future<Eigen::VectorXd>> task_vec;
        // Multiply separately, and concatenate chunks
        for (size_t i = 0; i < mats.size(); i++) {
            task_vec.push_back(std::async(
                std::launch::deferred,
                &MatrixLoader<T>::vecMultiplyRight,
                mats[i].get(),
                v,
                user_interrupt
            ));
        }
        std::vector<Eigen::VectorXd> sub_results = parallel_map_helper(task_vec, threads);
        Eigen::VectorXd res(rows());
        if (user_interrupt != NULL && *user_interrupt) return res;
        for (size_t i = 0; i < mats.size(); i++) {
            res.middleRows(row_offset[i], mats[i]->rows()) = sub_results[i];
        }
        return res;
    }

    Eigen::MatrixXd
    denseMultiplyLeft(const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt) override {
        if (rows() != B.cols()) throw std::runtime_error("Incompatible dimensions for matrix multiply");
        std::vector<std::future<Eigen::MatrixXd>> task_vec;
        // Multiply chunks, and add outputs
        for (size_t i = 0; i < mats.size(); i++) {
            task_vec.push_back(std::async(
                std::launch::deferred,
                &MatrixLoader<T>::denseMultiplyLeft,
                mats[i].get(),
                subset_map_cols(B, row_offset[i], mats[i]->rows()),
                user_interrupt
            ));
        }
        std::vector<Eigen::MatrixXd> sub_results = parallel_map_helper(task_vec, threads);
        Eigen::MatrixXd res(B.rows(), cols());
        if (user_interrupt != NULL && *user_interrupt) return res;
        res.setZero();
        for (size_t i = 0; i < mats.size(); i++) {
            res += sub_results[i];
        }
        return res;
    }

    Eigen::VectorXd
    vecMultiplyLeft(const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt) override {
        if (rows() != v.rows()) throw std::runtime_error("Incompatible dimensions for vector multiply");
        std::vector<std::future<Eigen::VectorXd>> task_vec;
        // Multiply chunks, and add outputs
        for (size_t i = 0; i < mats.size(); i++) {
            task_vec.push_back(std::async(
                std::launch::deferred,
                &MatrixLoader<T>::vecMultiplyLeft,
                mats[i].get(),
                subset_map_vec(v, row_offset[i], mats[i]->rows()),
                user_interrupt
            ));
        }
        std::vector<Eigen::VectorXd> sub_results = parallel_map_helper(task_vec, threads);
        Eigen::VectorXd res(cols());
        if (user_interrupt != NULL && *user_interrupt) return res;
        res.setZero();
        for (size_t i = 0; i < mats.size(); i++) {
            res += sub_results[i];
        }
        return res;
    }

    std::vector<T> colSums(std::atomic<bool> *user_interrupt) override {
        std::vector<std::future<std::vector<T>>> task_vec;
        for (size_t i = 0; i < mats.size(); i++) {
            task_vec.push_back(std::async(
                std::launch::deferred, &MatrixLoader<T>::colSums, mats[i].get(), user_interrupt
            ));
        }
        std::vector<std::vector<T>> sub_results = parallel_map_helper(task_vec, threads);
        std::vector<T> res(cols(), 0);
        if (user_interrupt != NULL && *user_interrupt) return res;
        for (size_t col = 0; col < cols(); col++) {
            for (size_t i = 0; i < mats.size(); i++) {
                res[col] += sub_results[i][col];
            }
        }
        return res;
    }

    std::vector<T> rowSums(std::atomic<bool> *user_interrupt) override {
        std::vector<std::future<std::vector<T>>> task_vec;
        for (size_t i = 0; i < mats.size(); i++) {
            task_vec.push_back(std::async(
                std::launch::deferred, &MatrixLoader<T>::rowSums, mats[i].get(), user_interrupt
            ));
        }
        std::vector<std::vector<T>> sub_results = parallel_map_helper(task_vec, threads);
        std::vector<T> res;
        if (user_interrupt != NULL && *user_interrupt) return res;
        for (size_t i = 0; i < mats.size(); i++) {
            res.insert(res.end(), sub_results[i].begin(), sub_results[i].end());
        }
        return res;
    }

    StatsResult
    computeMatrixStats(Stats row_stats, Stats col_stats, std::atomic<bool> *user_interrupt) override {
        std::vector<std::future<StatsResult>> task_vec;
        // Multiply chunks, and add outputs
        for (size_t i = 0; i < mats.size(); i++) {
            task_vec.push_back(std::async(
                std::launch::deferred,
                &MatrixLoader<T>::computeMatrixStats,
                mats[i].get(),
                row_stats,
                col_stats,
                user_interrupt
            ));
        }
        std::vector<StatsResult> sub_results = parallel_map_helper(task_vec, threads);
        StatsResult res{
            Eigen::ArrayXXd((int)row_stats, rows()), Eigen::ArrayXXd((int)col_stats, cols())};
        res.row_stats.setZero();
        res.col_stats.setZero();
        // Concatenate the row results
        for (size_t i = 0; i < mats.size(); i++) {
            res.row_stats.middleCols(row_offset[i], mats[i]->rows()) = sub_results[i].row_stats;
        }
        // Smarter math for the col results
        for (size_t i = 0; i < mats.size(); i++) {
            // Nonzeros simply add
            if (res.col_stats.rows() > 0) {
                res.col_stats.row(0) += sub_results[i].col_stats.row(0);
            }
            // Means get scaled by the fraction of rows, then added
            if (res.col_stats.rows() > 1) {
                double scale_factor = mats[i]->rows() / ((double)rows());
                res.col_stats.row(1) += scale_factor * sub_results[i].col_stats.row(1);
            }
        }
        if (res.col_stats.rows() <= 2) return res;

        // For variances, we want to have already calculated the correct mean
        for (size_t i = 0; i < mats.size(); i++) {
            // 1. Get sum((val - orig_mean)^2)
            auto m2_orig = sub_results[i].col_stats.row(2) * (mats[i]->rows() - 1);
            // 2. Get sum((val - new_mean)^2) = sum((val - orig_mean)^2) + N * (new_mean -
            // orig_mean)^2
            auto m2_new =
                m2_orig +
                mats[i]->rows() * (res.col_stats.row(1) - sub_results[i].col_stats.row(1)).square();
            // 3. Add in to our running sum((val - new_mean)^2)
            res.col_stats.row(2) += m2_new;
        }
        // Convert sum((val - new_mean)^2) into sample variance
        res.col_stats.row(2) /= (rows() - 1);
        return res;
    }
};

// Concatenate a list of MatrixLoaders by col.
// Row names will be taken from the first matrix in the list
template <typename T> class ConcatCols : public MatrixLoader<T> {
  protected:
    std::vector<std::unique_ptr<MatrixLoader<T>>> mats;
    std::vector<uint32_t> col_offset;
    uint32_t cur_mat = 0;
    uint32_t threads = 0;
  public:
    ConcatCols(std::vector<std::unique_ptr<MatrixLoader<T>>> &&mats, uint32_t threads) : mats(std::move(mats)), threads(threads) {

        if (this->mats.size() < 2) throw std::runtime_error("Must have >= 2 matrices to merge");

        col_offset.push_back(0);
        uint32_t rows = this->mats.front()->rows();
        for (const auto &m : this->mats) {
            col_offset.push_back(col_offset.back() + m->cols());
            if (m->rows() != rows)
                throw std::runtime_error("ConcatCols: Matrices must have equal numbers of rows");
        }

        this->threads = std::min(this->threads, (uint32_t) this->mats.size());
    }

    uint32_t rows() const override { return mats.front()->rows(); }
    uint32_t cols() const override { return col_offset.back(); }

    const char *rowNames(uint32_t row) override { return mats.front()->rowNames(row); }

    const char *colNames(uint32_t col) override {
        auto it = std::upper_bound(col_offset.begin(), col_offset.end(), col);
        uint32_t idx = it - col_offset.begin() - 1;

        if (idx == mats.size()) return NULL;

        return mats[idx]->colNames(col - col_offset[idx]);
    }

    // Reset the iterator to start from the beginning
    void restart() override {
        cur_mat = 0;
        for (auto &m : mats)
            m->restart();
    }

    // Seek to a specific column without reading data
    // Next call should be to load(). col must be < cols()
    void seekCol(uint32_t col) override {
        auto it = std::upper_bound(col_offset.begin(), col_offset.end(), col);
        uint32_t idx = it - col_offset.begin() - 1;

        if (idx == mats.size())
            throw std::runtime_error(
                "ConcatCols: Cannot seek to a column larger than number of columns"
            );
        mats[idx]->seekCol(col - col_offset[idx]);
        cur_mat = idx;
    }

    // Advance to the next column, or return false if there
    // are no more columns
    bool nextCol() override {
        if (!mats[cur_mat]->nextCol()) {
            if (cur_mat + 1 == mats.size()) return false;
            cur_mat++;
            mats[cur_mat]->seekCol(0);
        }

        return true;
    }

    // Return the index of the current column
    uint32_t currentCol() const override {
        return mats[cur_mat]->currentCol() + col_offset[cur_mat];
    }

    // Return false if there are no more entries to load
    bool load() override { return mats[cur_mat]->load(); }

    // Number of loaded entries available
    uint32_t capacity() const override { return mats[cur_mat]->capacity(); }

    // Pointers to the loaded entries
    uint32_t *rowData() override { return mats[cur_mat]->rowData(); }
    T *valData() override { return mats[cur_mat]->valData(); }

    Eigen::MatrixXd
    denseMultiplyLeft(const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt) override {
        if (rows() != B.cols()) throw std::runtime_error("Incompatible dimensions for matrix multiply");
        std::vector<std::future<Eigen::MatrixXd>> task_vec;
        // Multiply separately, and concatenate chunks
        for (size_t i = 0; i < mats.size(); i++) {
            task_vec.push_back(std::async(
                std::launch::deferred,
                &MatrixLoader<T>::denseMultiplyLeft,
                mats[i].get(),
                B,
                user_interrupt
            ));
        }
        std::vector<Eigen::MatrixXd> sub_results = parallel_map_helper(task_vec, threads);
        Eigen::MatrixXd res(B.rows(), cols());
        if (user_interrupt != NULL && *user_interrupt) return res;
        for (size_t i = 0; i < mats.size(); i++) {
            res.middleCols(col_offset[i], mats[i]->cols()) = sub_results[i];
        }
        return res;
    }

    Eigen::VectorXd
    vecMultiplyLeft(const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt) override {
        if (rows() != v.rows()) throw std::runtime_error("Incompatible dimensions for vector multiply");
        std::vector<std::future<Eigen::VectorXd>> task_vec;
        // Multiply separately, and concatenate chunks
        for (size_t i = 0; i < mats.size(); i++) {
            task_vec.push_back(std::async(
                std::launch::deferred,
                &MatrixLoader<T>::vecMultiplyLeft,
                mats[i].get(),
                v,
                user_interrupt
            ));
        }
        std::vector<Eigen::VectorXd> sub_results = parallel_map_helper(task_vec, threads);
        Eigen::VectorXd res(cols());
        if (user_interrupt != NULL && *user_interrupt) return res;
        for (size_t i = 0; i < mats.size(); i++) {
            res.middleRows(col_offset[i], mats[i]->cols()) = sub_results[i];
        }
        return res;
    }

    Eigen::MatrixXd
    denseMultiplyRight(const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt) override {
        if (cols() != B.rows()) throw std::runtime_error("Incompatible dimensions for matrix multiply");
        std::vector<std::future<Eigen::MatrixXd>> task_vec;
        std::vector<Eigen::MatrixXd> B_chunks;
        // Multiply chunks, and add outputs
        for (size_t i = 0; i < mats.size(); i++) {
            B_chunks.push_back(B.middleRows(col_offset[i], mats[i]->cols()));
            task_vec.push_back(std::async(
                std::launch::deferred,
                &MatrixLoader<T>::denseMultiplyRight,
                mats[i].get(),
                Eigen::Map<Eigen::MatrixXd>(B_chunks[i].data(), B_chunks[i].rows(), B_chunks[i].cols()),
                user_interrupt
            ));
        }
        std::vector<Eigen::MatrixXd> sub_results = parallel_map_helper(task_vec, threads);
        Eigen::MatrixXd res(rows(), B.cols());
        if (user_interrupt != NULL && *user_interrupt) return res;
        res.setZero();
        for (size_t i = 0; i < mats.size(); i++) {
            res += sub_results[i];
        }
        return res;
    }

    Eigen::VectorXd
    vecMultiplyRight(const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt) override {
        if (cols() != v.rows()) throw std::runtime_error("Incompatible dimensions for vector multiply");
        std::vector<std::future<Eigen::VectorXd>> task_vec;
        // Multiply chunks, and add outputs
        for (size_t i = 0; i < mats.size(); i++) {
            task_vec.push_back(std::async(
                std::launch::deferred,
                &MatrixLoader<T>::vecMultiplyRight,
                mats[i].get(),
                subset_map_vec(v, col_offset[i], mats[i]->cols()),
                user_interrupt
            ));
        }
        std::vector<Eigen::VectorXd> sub_results = parallel_map_helper(task_vec, threads);
        Eigen::VectorXd res(rows());
        if (user_interrupt != NULL && *user_interrupt) return res;
        res.setZero();
        for (size_t i = 0; i < mats.size(); i++) {
            res += sub_results[i];
        }
        return res;
    }

    std::vector<T> rowSums(std::atomic<bool> *user_interrupt) override {
        std::vector<std::future<std::vector<T>>> task_vec;
        for (size_t i = 0; i < mats.size(); i++) {
            task_vec.push_back(std::async(
                std::launch::deferred, &MatrixLoader<T>::rowSums, mats[i].get(), user_interrupt
            ));
        }
        std::vector<std::vector<T>> sub_results = parallel_map_helper(task_vec, threads);
        std::vector<T> res(rows(), 0);
        if (user_interrupt != NULL && *user_interrupt) return res;
        for (size_t row = 0; row < rows(); row++) {
            for (size_t i = 0; i < mats.size(); i++) {
                res[row] += sub_results[i][row];
            }
        }
        return res;
    }

    std::vector<T> colSums(std::atomic<bool> *user_interrupt) override {
        std::vector<std::future<std::vector<T>>> task_vec;
        for (size_t i = 0; i < mats.size(); i++) {
            task_vec.push_back(std::async(
                std::launch::deferred, &MatrixLoader<T>::colSums, mats[i].get(), user_interrupt
            ));
        }
        std::vector<std::vector<T>> sub_results = parallel_map_helper(task_vec, threads);
        std::vector<T> res;
        if (user_interrupt != NULL && *user_interrupt) return res;
        for (size_t i = 0; i < mats.size(); i++) {
            res.insert(res.end(), sub_results[i].begin(), sub_results[i].end());
        }
        return res;
    }

    StatsResult
    computeMatrixStats(Stats row_stats, Stats col_stats, std::atomic<bool> *user_interrupt) override {
        std::vector<std::future<StatsResult>> task_vec;
        // Multiply chunks, and add outputs
        for (size_t i = 0; i < mats.size(); i++) {
            task_vec.push_back(std::async(
                std::launch::deferred,
                &MatrixLoader<T>::computeMatrixStats,
                mats[i].get(),
                row_stats,
                col_stats,
                user_interrupt
            ));
        }
        std::vector<StatsResult> sub_results = parallel_map_helper(task_vec, threads);
        StatsResult res{
            Eigen::ArrayXXd((int)row_stats, rows()), Eigen::ArrayXXd((int)col_stats, cols())};
        res.row_stats.setZero();
        res.col_stats.setZero();
        // Concatenate the col results
        for (size_t i = 0; i < mats.size(); i++) {
            res.col_stats.middleCols(col_offset[i], mats[i]->cols()) = sub_results[i].col_stats;
        }
        // Smarter math for the row results
        for (size_t i = 0; i < mats.size(); i++) {
            // Nonzeros simply add
            if (res.row_stats.rows() > 0) {
                res.row_stats.row(0) += sub_results[i].row_stats.row(0);
            }
            // Means get scaled by the fraction of rows, then added
            if (res.row_stats.rows() > 1) {
                double scale_factor = mats[i]->cols() / ((double)cols());
                res.row_stats.row(1) += scale_factor * sub_results[i].row_stats.row(1);
            }
        }
        if (res.row_stats.rows() <= 2) return res;

        // For variances, we want to have already calculated the correct mean
        for (size_t i = 0; i < mats.size(); i++) {
            // 1. Get sum((val - orig_mean)^2)
            auto m2_orig = sub_results[i].row_stats.row(2) * (mats[i]->cols() - 1);
            // 2. Get sum((val - new_mean)^2) = sum((val - orig_mean)^2) + N * (new_mean -
            // orig_mean)^2
            auto m2_new =
                m2_orig +
                mats[i]->cols() * (res.row_stats.row(1) - sub_results[i].row_stats.row(1)).square();
            // 3. Add in to our running sum((val - new_mean)^2)
            res.row_stats.row(2) += m2_new;
        }
        // Convert sum((val - new_mean)^2) into sample variance
        res.row_stats.row(2) /= (cols() - 1);
        return res;
    }
};

} // end namespace BPCells