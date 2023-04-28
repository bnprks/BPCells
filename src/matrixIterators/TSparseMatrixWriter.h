#pragma once

#include <atomic>
#include "MatrixIterator.h"

namespace BPCells {

template <typename T> class TSparseMatrixWriter : public MatrixWriter<T> {
  public:
    bool write(MatrixLoader<T> &mat, std::atomic<bool> *user_interrupt = NULL) override {
        MatrixIterator<T> it(mat);
        uint32_t count = 0;
        while (it.nextCol()) {
            while (it.nextValue()) {
                rows.push_back(it.row());
                cols.push_back(it.col());
                vals.push_back(it.val());
                if (count++ % 8192 == 0 && user_interrupt != NULL && *user_interrupt) return false;
            }
        }
        return true;
    };

    std::vector<uint32_t> rows;
    std::vector<uint32_t> cols;
    std::vector<T> vals;
};

} // end namespace BPCells