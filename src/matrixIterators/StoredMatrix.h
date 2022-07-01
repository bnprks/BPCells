#pragma once

#include "MatrixIterator.h"
#include "../arrayIO/array_interfaces.h"
#include "../arrayIO/bp128.h"

namespace BPCells {

// Main class for accessing matrices stored on disk.
// Templated to help with compatibility reading 10x and AnnData matrix formats.
template<class T>
class StoredMatrix: public MatrixLoader<T> {
private:
    NumReader<uint32_t> row;
    NumReader<T> val;
    NumReader<uint32_t> col_ptr;
    std::unique_ptr<StringReader> row_names, col_names;
    uint32_t n_rows;
    uint32_t n_cols;
    uint32_t current_col = UINT32_MAX;
    uint32_t current_idx = 0;
    uint32_t next_col_ptr;
    uint32_t current_capacity = 0;

public:
    StoredMatrix() = default;
    StoredMatrix(StoredMatrix &&other) = default;
    StoredMatrix& operator=(StoredMatrix &&other) = default;
    StoredMatrix(NumReader<uint32_t> &&row, NumReader<T> &&val, NumReader<uint32_t> &&col_ptr, 
        NumReader<uint32_t> &&row_count, 
        std::unique_ptr<StringReader> &&row_names,
        std::unique_ptr<StringReader> &&col_names) :
        row(std::move(row)), val(std::move(val)), col_ptr(std::move(col_ptr)), 
        row_names(std::move(row_names)), col_names(std::move(col_names)),
        n_rows(row_count.read_one()), n_cols(this->col_ptr.size() - 1),
        next_col_ptr(this->col_ptr.read_one()) {}

    static StoredMatrix<uint32_t> openUnpacked(ReaderBuilder &rb, std::unique_ptr<StringReader> &&row_names, std::unique_ptr<StringReader> &&col_names) {
        if (rb.readVersion() != "unpacked-uint-matrix-v1") {
            throw std::runtime_error(std::string("Version does not match unpacked-uint-matrix-v1: ") + rb.readVersion());
        }

        return StoredMatrix<uint32_t>(
            rb.openUIntReader("row"),
            rb.openUIntReader("val"),
            rb.openUIntReader("col_ptr"),
            rb.openUIntReader("row_count"),
            std::move(row_names),
            std::move(col_names)
        );
    }
    static StoredMatrix<uint32_t> openUnpacked(ReaderBuilder &rb) {
        return StoredMatrix<uint32_t>::openUnpacked(rb, rb.openStringReader("row_names"), rb.openStringReader("col_names"));
    }


    static StoredMatrix<uint32_t> openPacked(ReaderBuilder &rb, uint32_t load_size, std::unique_ptr<StringReader> &&row_names, std::unique_ptr<StringReader> &&col_names) {
        if (rb.readVersion() != "packed-uint-matrix-v1") {
            throw std::runtime_error(std::string("Version does not match packed-uint-matrix-v1: ") + rb.readVersion());
        }

        UIntReader col_ptr = rb.openUIntReader("col_ptr");
        col_ptr.seek(col_ptr.size() - 1);
        uint32_t count = col_ptr.read_one();
        col_ptr.seek(0);

        return StoredMatrix(
            UIntReader(
                std::make_unique<BP128_D1Z_UIntReader>(
                    rb.openUIntReader("row_data"),
                    rb.openUIntReader("row_idx"),
                    rb.openUIntReader("row_starts"),
                    count
                ),
                load_size, load_size
            ),
            UIntReader(
                std::make_unique<BP128_FOR_UIntReader>(
                    rb.openUIntReader("val_data"),
                    rb.openUIntReader("val_idx"),
                    count
                ),
                load_size, load_size
            ),
            std::move(col_ptr),
            rb.openUIntReader("row_count"),
            std::move(row_names),
            std::move(col_names)
        );
    }
    static StoredMatrix<uint32_t> openPacked(ReaderBuilder &rb, uint32_t load_size=1024) {
        return StoredMatrix<uint32_t>::openPacked(rb, load_size, rb.openStringReader("row_names"), rb.openStringReader("col_names"));
    }

    // Return the count of rows and columns
    uint32_t rows() const override {return n_rows;}
    uint32_t cols() const override {return n_cols;}

    const char* rowNames(uint32_t row) const override {return row_names->get(row);}
    const char* colNames(uint32_t col) const override {return col_names->get(col);}

    // Reset the iterator to start from the beginning
    void restart() override {
        current_col = UINT32_MAX;
        // Don't change current_idx so we will correctly seek when nextCol iscalled
        current_capacity = 0;
        col_ptr.seek(0);
        next_col_ptr = col_ptr.read_one();
    }

    // Seek to a specific column without reading data
    void seekCol(uint32_t col) override {
        current_col = col-1;
        col_ptr.seek(col);
        next_col_ptr = col_ptr.read_one();
        nextCol();
    }

    // Advance to the next column, or return false if there
    // are no more columns
    bool nextCol() override {
        current_col += 1;
        if (current_col >= n_cols) {
            current_col -= 1;
            return false;
        }

        // Check if we need to perform any seeks, or if we're all good
        // since we just finished a column
        if (next_col_ptr != current_idx) {
            // We need to perform seeks to get to the right data
            // reading location
            current_idx = next_col_ptr;
            val.seek(current_idx);
            row.seek(current_idx);
        }
        
        next_col_ptr = col_ptr.read_one();
        current_capacity = 0;
        return true;
    }

    // Return the index of the current column
    uint32_t currentCol() const override {return current_col;}
    
    // Return false if there are no more entries to load
    bool load() override {
        val.advance(current_capacity);
        row.advance(current_capacity);

        if (current_idx >= next_col_ptr) {
            current_capacity = 0;
            return false;
        }

        // Load more data if necessary
        if (val.capacity() == 0) val.ensureCapacity(1);
        if (row.capacity() == 0) row.ensureCapacity(1);

        current_capacity = std::min({val.capacity(), row.capacity(), next_col_ptr-current_idx});
        current_idx += current_capacity;
        return true;
    }

    // Number of loaded entries available
    uint32_t capacity() const override {return current_capacity;}

    // Pointers to the loaded entries
    uint32_t* rowData() override {return row.data();}
    T* valData() override {return val.data();}
};

class StoredMatrixWriter: public MatrixWriter<uint32_t> {
private:
    UIntWriter row, val, col_ptr, row_count;
    std::unique_ptr<StringWriter> row_names, col_names;
public:
    static StoredMatrixWriter createUnpacked(WriterBuilder &wb);
    static StoredMatrixWriter createPacked(WriterBuilder &wb, uint32_t buffer_size=1024);
    StoredMatrixWriter(UIntWriter &&row, UIntWriter &&val, UIntWriter &&col_ptr, 
        UIntWriter &&row_count, 
        std::unique_ptr<StringWriter> &&row_names, std::unique_ptr<StringWriter> &&col_names);  
    void write(MatrixLoader<uint32_t> &loader, void (*checkInterrupt)(void) = NULL) override;
};

} // end namespace BPCells