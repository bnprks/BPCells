#pragma once

#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>

#include "MatrixIterator.h"
#include "../arrayIO/array_interfaces.h"

namespace BPCells {
// Base class to load entries from (sparse) matrices
// Subclasses implement reading from various in-memory and disk formats,
// calculation of matrices from fragments, as well as transformations like
// normalization, and re-indexing.
//
// For simplicity, matrix loaders are purely in column-major order, similar to
// Eigen and R's default compressed sparse column layouts. 
// Transposition can be performed only with a MatrixWriter,
// since the re-ordering requires storing all entries in an intermediate matrix.
//
// For flexibility with indexing, matrix iterators need only be grouped by column,
// and neither column or row ordering matters, so long as all entries for a single column
// are consecutive
// To implement a new MatrixIterator:
// 1. Implement rows() and cols() methods to return dimension
// 2. Implement load() method to load the next chunk from the current column.
//    this should return 0 repeatedly at the end of a column until nextCol is called.
// 3. Implement nextCol() method to advance to the next available column.
// 4. Implement restart() method to restart the iterator from the beginning
class UnpackedMatrix: public MatrixLoader<uint32_t> {
private:
    using ReaderPtr = std::unique_ptr<UIntReader>;
    ReaderPtr val, row, col_ptr, row_count;
    uint32_t n_rows;
    uint32_t n_cols;
    uint32_t current_col = UINT32_MAX;
    uint32_t current_idx = UINT32_MAX;
    uint32_t next_col_ptr;

public:
    UnpackedMatrix(ReaderPtr &&val, ReaderPtr &&row, 
                ReaderPtr &&col_ptr, ReaderPtr &&row_count);
                
    // Return the count of rows and columns
    uint32_t rows() const override;
    uint32_t cols() const override;

    // Reset the iterator to start from the beginning
    void restart() override;

    // Advance to the next column, or return false if there
    // are no more columns
    bool nextCol() override;

    // Return the index of the current column
    uint32_t currentCol() const override;
    
    // Return number of matrix entries loaded. Should repeatedly return 0 
    // at the end of a column, until nextCol is called
    // Return -1 for error
    int32_t load(uint32_t count, SparseVector<uint32_t> buffer) override;
};

class UnpackedMatrixWriter: public MatrixWriter<uint32_t> {
private: 
    using WriterPtr = std::unique_ptr<UIntWriter>;
    WriterPtr val, row, col_ptr, row_count;
            
public:
    UnpackedMatrixWriter(WriterPtr &&val, WriterPtr &&row, 
                       WriterPtr &&col_ptr, WriterPtr &&row_count);
    // Return false on failure, true on success
    bool write(MatrixLoader<uint32_t> &mat, void (*checkInterrupt)(void) = NULL) override;
};

} // end namespace BPCells