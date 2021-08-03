#include "MatrixIterator.h"

#include "../packed_fragments.h"

namespace BPCells {

class PackedMatrixIterator : public MatrixIterator<uint32_t> {
public:
    PackedMatrixIterator(MatrixIterator<T> &iter, uint32_t buffer_size = 1024) : 
        MatrixIterator<T>(buffer_size), iter(iter) {}

    // Reset the iterator to start from the beginning
    void restart() override { iter.restart(); }

    // Return the count of rows and columns
    uint32_t rowCount() override { return iter.colCount(); }
    uint32_t colCount() override { return iter.rowCount(); }

    int32_t load(uint32_t count, MatrixArray<T> buffer) override {
        MatrixArray<T> inner_buffer;
        inner_buffer.capacity = buffer.capacity;
        inner_buffer.row = buffer.col;
        inner_buffer.col = buffer.row;
        inner_buffer.val = buffer.val;

        return iter.load(count, inner_buffer);
    };
};

class PackedMatrixWriter : public MatrixWriter<uint32_t> {

public:
    // Return false on failure, true on success
    bool write(MatrixIterator<T> &mat, void (*checkInterrupt)(void) = NULL) override {
        std::vector<uint32_t> col_buffer;
        while(mat.nextCol()) {
            while(mat.nextValue()) {

            }
        }
    }
};


} // end namespace BPCells