#pragma once

#include "MatrixIterator.h"

namespace BPCells {

template <typename T>
class TSparseMatrixWriter : public MatrixWriter<T> {
public:
    bool write(MatrixIterator<T> &mat, void (*checkInterrupt)(void) = NULL) override {
        uint32_t count = 0;
        while (mat.nextCol()) {
            while (mat.nextValue()) {
                rows.push_back(mat.row());
                cols.push_back(mat.col());
                vals.push_back(mat.val());
                if(count++ % 10000 == 0 && checkInterrupt != NULL) checkInterrupt();
            }
        }
        return true;
    };

    std::vector<uint32_t> rows;
    std::vector<uint32_t> cols;
    std::vector<T> vals;
};

} // end namespace BPCells