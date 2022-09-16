#pragma once

#include "MatrixIterator.h"

namespace BPCells {

// Concatenate a list of MatrixLoaders by row.
// Column names will be taken from the first matrix in the list
template <typename T> class ConcatRows : public MatrixLoader<T> {
  protected:
    std::vector<MatrixLoader<T> *> mats;
    std::vector<uint32_t> row_offset;
    uint32_t cur_mat = 0;

  public:
    ConcatRows(const std::vector<MatrixLoader<T> *> &mats) : mats(mats) {

        if (mats.size() < 2) throw std::runtime_error("Must have >= 2 matrices to merge");

        row_offset.push_back(0);
        uint32_t cols = mats.front()->cols();
        for (auto m : this->mats) {
            row_offset.push_back(row_offset.back() + m->rows());
            if (m->cols() != cols)
                throw std::runtime_error("ConcatRows: Matrices must have equal numbers of columns");
        }
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
};

// Concatenate a list of MatrixLoaders by col.
// Row names will be taken from the first matrix in the list
template <typename T> class ConcatCols : public MatrixLoader<T> {
  protected:
    std::vector<MatrixLoader<T> *> mats;
    std::vector<uint32_t> col_offset;
    uint32_t cur_mat = 0;

  public:
    ConcatCols(const std::vector<MatrixLoader<T> *> &mats) : mats(mats) {

        if (mats.size() < 2) throw std::runtime_error("Must have >= 2 matrices to merge");

        col_offset.push_back(0);
        uint32_t rows = mats.front()->rows();
        for (auto m : this->mats) {
            col_offset.push_back(col_offset.back() + m->cols());
            if (m->rows() != rows)
                throw std::runtime_error("ConcatCols: Matrices must have equal numbers of rows");
        }
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
        while (!mats[cur_mat]->nextCol()) {
            if (cur_mat + 1 == mats.size()) return false;
            cur_mat++;
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
};

} // end namespace BPCells