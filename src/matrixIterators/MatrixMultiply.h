#pragma once

#include "MatrixIterator.h"

namespace BPCells {

// Perform sparse matrix multiplication given two input matrices
template <typename T> class SparseMultiply : public MatrixLoader<T> {
  protected:
    // Running invariants:
    // - right is always about to read the next column for output
    // - row_idx is the next row to resume output from
    MatrixLoader<T> &left, &right;
    std::vector<T> col_buffer, val_buffer;
    std::vector<uint32_t> row_buffer;
    uint32_t row_idx = 0;
    uint32_t max_load_size, loaded;

  public:
    SparseMultiply(MatrixLoader<T> &left, MatrixLoader<T> &right, uint32_t load_size = 1024)
        : left(left)
        , right(right)
        , col_buffer(left.rows())
        , val_buffer(load_size)
        , row_buffer(load_size)
        , max_load_size(load_size) {

        if (left.cols() != right.rows())
            throw std::runtime_error("Matrices have incompatible dimensions for multiplication");
    }

    uint32_t rows() const override { return left.rows(); }
    uint32_t cols() const override { return right.cols(); }

    const char *rowNames(uint32_t row) override { return left.rowNames(row); }
    const char *colNames(uint32_t col) override { return right.colNames(col); }

    // Reset the iterator to start from the beginning
    void restart() override {
        right.restart();
        row_idx = 0;
    }

    // Seek to a specific column without reading data
    // Next call should be to load(). col must be < cols()
    void seekCol(uint32_t col) override {
        right.seekCol(col);
        row_idx = 0;
    }

    // Advance to the next column, or return false if there
    // are no more columns
    bool nextCol() override {
        if (!right.nextCol()) return false;
        row_idx = 0;
        return true;
    }

    // Return the index of the current column
    uint32_t currentCol() const override { return right.currentCol(); }

    // Return false if there are no more entries to load
    bool load() override {
        if (row_idx == 0) {
            // Load the next column of data into col_buffer
            for (auto &x : col_buffer) {
                x = 0;
            }
            // For each entry in the current column of right, add entry_val * left[:,entry_row] to
            // the col_buffer
            while (right.load()) {
                uint32_t r_cap = right.capacity();
                for (uint32_t i = 0; i < r_cap; i++) {
                    if (left.currentCol() + 1 == right.rowData()[i]) left.nextCol();
                    else left.seekCol(right.rowData()[i]);

                    T val = right.valData()[i];
                    while (left.load()) {
                        uint32_t l_cap = left.capacity();
                        for (uint32_t j = 0; j < l_cap; j++) {
                            col_buffer[left.rowData()[j]] += left.valData()[j] * val;
                        }
                    }
                }
            }
        }

        loaded = 0;
        for (; row_idx < col_buffer.size(); row_idx++) {
            if (col_buffer[row_idx] == 0) continue;
            val_buffer[loaded] = col_buffer[row_idx];
            row_buffer[loaded] = row_idx;
            loaded += 1;
            if (loaded >= max_load_size) break;
        }
        return loaded > 0;
    }

    // Number of loaded entries available
    uint32_t capacity() const override { return loaded; }

    // Pointers to the loaded entries
    uint32_t *rowData() override { return row_buffer.data(); }
    T *valData() override { return val_buffer.data(); }

    // MATH OPERATIONS: utilize associative property that A*B*C = (A*B)*C = A*(B*C)
    // Calculate matrix-matrix product A*B where A (this) is sparse and B is a dense matrix.
    Eigen::MatrixXd denseMultiplyRight(
        const Eigen::Map<Eigen::MatrixXd> B, void (*checkInterrupt)(void) = NULL
    ) override {
        auto tmp = right.denseMultiplyRight(B, checkInterrupt);
        Eigen::Map<Eigen::MatrixXd> map(tmp.data(), tmp.rows(), tmp.cols());
        return left.denseMultiplyRight(map, checkInterrupt);
    }
    Eigen::MatrixXd denseMultiplyLeft(
        const Eigen::Map<Eigen::MatrixXd> B, void (*checkInterrupt)(void) = NULL
    ) override {
        auto tmp = left.denseMultiplyLeft(B, checkInterrupt);
        Eigen::Map<Eigen::MatrixXd> map(tmp.data(), tmp.rows(), tmp.cols());
        return right.denseMultiplyLeft(map, checkInterrupt);
    }
    // Calculate matrix-vector product A*v where A (this) is sparse and B is a dense matrix.
    Eigen::VectorXd vecMultiplyRight(
        const Eigen::Map<Eigen::VectorXd> v, void (*checkInterrupt)(void) = NULL
    ) override {
        auto tmp = right.vecMultiplyRight(v, checkInterrupt);
        Eigen::Map<Eigen::VectorXd> map(tmp.data(), tmp.rows(), tmp.cols());
        return left.vecMultiplyRight(map, checkInterrupt);
    }
    Eigen::VectorXd vecMultiplyLeft(
        const Eigen::Map<Eigen::VectorXd> v, void (*checkInterrupt)(void) = NULL
    ) override {
        auto tmp = left.vecMultiplyLeft(v, checkInterrupt);
        Eigen::Map<Eigen::VectorXd> map(tmp.data(), tmp.rows(), tmp.cols());
        return right.vecMultiplyLeft(map, checkInterrupt);
    }
    // Calculate row/column sums of the matrix
    std::vector<T> colSums(void (*checkInterrupt)(void) = NULL) override {
        Eigen::VectorXd v;
        v.setOnes(rows());
        Eigen::Map<Eigen::VectorXd> map(v.data(), v.rows(), v.cols());
        auto res = vecMultiplyLeft(map, checkInterrupt);
        std::vector<T> ret(cols());
        for (uint32_t i = 0; i < cols(); i++) {
            ret[i] = res[i];
        }
        return ret;
    }
    std::vector<T> rowSums(void (*checkInterrupt)(void) = NULL) override {
        Eigen::VectorXd v;
        v.setOnes(cols());
        Eigen::Map<Eigen::VectorXd> map(v.data(), v.rows(), v.cols());
        auto res = vecMultiplyRight(map, checkInterrupt);
        std::vector<T> ret(rows());
        for (uint32_t i = 0; i < rows(); i++) {
            ret[i] = res[i];
        }
        return ret;
    }
};

} // end namespace BPCells