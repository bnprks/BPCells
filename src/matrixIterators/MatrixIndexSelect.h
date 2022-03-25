#pragma once
#include <unordered_map>
#include "MatrixIterator.h"
#include "../radix_sort.h"

namespace BPCells {

// Select specific columns from a dataset 
template<class T>
class MatrixColSelect : public MatrixLoaderWrapper<T> {
private:
    uint32_t current_col = UINT32_MAX;
    const std::vector<uint32_t> col_indices;
public:
    // col_indices -- vector of columns to select
    MatrixColSelect(MatrixLoader<T> &loader, const std::vector<uint32_t> col_indices) :
        MatrixLoaderWrapper<T>(loader), col_indices(col_indices) {}

    ~MatrixColSelect() = default;

    uint32_t cols() const override {return col_indices.size();}
    const char* colNames(uint32_t col) const override {
        if (col < col_indices.size()) return this->loader.colNames(col_indices[col]);
        return NULL;
    }

    void restart() override {current_col = 0; this->loader.restart();}
    void seekCol(uint32_t col) override {
        if (col >= col_indices.size()) throw std::runtime_error("Requested column is greater than number of columns");
        this->loader.seekCol(col_indices[col]);
        current_col = col_indices[col];
    }

    bool nextCol() override {
        current_col += 1;
        if (current_col > 0 && col_indices[current_col - 1] == col_indices[current_col] - 1) {
            return this->loader.nextCol();
        } else if (current_col >= col_indices.size()) {
            current_col -= 1;
            return false;
        } else {
            this->loader.seekCol(col_indices[current_col]);
            return true;
        }
    }
    uint32_t currentCol() const override {return current_col;}
};

// Select specific rows from a dataset
// Note: unlike the MatrixColSelect, MatrixRowSelect does not support duplicating rows,
// only filtering and/or reordering
template<class T>
class MatrixRowSelect : public MatrixLoaderWrapper<T> {
private:
    std::vector<uint32_t> row_data, row_buf;
    std::vector<T> val_data, val_buf;
    
    uint32_t idx = 0;
    uint32_t cap = 0;
    uint32_t load_size;
    
    // Reverse lookup for row indices -- reverse_indices[i] gives the output row_id
    // for input row_id i
    std::vector<uint32_t> reverse_indices;
    const std::vector<uint32_t> row_indices;
    bool needs_reorder = false;

public:
    // cell_names -- vector with length <= the number of chromosomes in the input
    //     FragmentLoader. The output cell `i` will come from input cell
    //     `chr_assignments[i]`. The entries of cell_names must be unique
    MatrixRowSelect(MatrixLoader<T> &loader, const std::vector<uint32_t> row_indices, uint32_t load_size=1024) :
        MatrixLoaderWrapper<T>(loader), load_size(load_size), reverse_indices(loader.rows(), UINT32_MAX), row_indices(row_indices) {
        for (uint32_t i = 0; i < row_indices.size(); i++) {
            if (reverse_indices[row_indices[i]] != UINT32_MAX)
                throw std::runtime_error("Cannot duplicate rows using MatrixRowSelect");
            reverse_indices[row_indices[i]] = i;
        }
        for (uint32_t i = 1; i < row_indices.size(); i++) {
            if (row_indices[i] < row_indices[i-1]) {
                needs_reorder = true;
                break;
            }
        }
        if (needs_reorder) {
            row_data.resize(this->loader.rows());
            row_buf.resize(this->loader.rows());
            val_data.resize(this->loader.rows());
            val_buf.resize(this->loader.rows());
        }
    }

    ~MatrixRowSelect() = default;
    
    uint32_t rows() const override {return row_indices.size();}
    const char* rowNames(uint32_t row) const override {
        if (row < row_indices.size()) return this->loader.rowNames(row_indices[row]);
        return NULL;
    }

    void restart() override {idx=0; cap=0; this->loader.restart();}
    void seekCol(uint32_t col) override {idx=0; cap=0; this->loader.seekCol(col);}

    bool nextCol() override {idx = 0; cap = 0; return this->loader.nextCol();}

    bool load() override {
        if (!needs_reorder) {
            // Just perform a straight filter and load incrementally
            uint32_t loaded = 0;
            while (loaded == 0) {
                if (!this->loader.load()) return false;
                
                cap = this->loader.capacity();
                uint32_t *row_ptr = this->loader.rowData();
                T *val_ptr = this->loader.valData();
                
                for (uint32_t i = 0; i < cap; i++) {
                    uint32_t new_row = reverse_indices[row_ptr[i]];
                    row_ptr[loaded] = new_row;
                    val_ptr[loaded] = val_ptr[i];
                    loaded += new_row != UINT32_MAX;
                }
            }
            cap = loaded;
            return true;
        } else if (idx == 0 && cap == 0) {
            // Load all the values for the column and sort
            while (this->loader.load()) {
                uint32_t loaded = this->loader.capacity();
                
                uint32_t *row_ptr = this->loader.rowData();
                T *val_ptr = this->loader.valData();
                
                for (uint32_t i = 0; i < loaded; i++) {
                    row_data[cap] = reverse_indices[row_ptr[i]];
                    val_data[cap] = val_ptr[i];
                    cap += reverse_indices[row_ptr[i]] != UINT32_MAX;
                }
            }
            // Sort the entries by increasing row ID
            lsdRadixSortArrays<T>(cap, row_data, val_data, row_buf, val_buf);
        } else {
            idx += std::min(cap - idx, load_size);
        }
        return idx < cap;
    }

    uint32_t capacity() const override {return cap-idx;}
    uint32_t* rowData() override {return needs_reorder ? row_data.data() + idx : this->loader.rowData();}
    T* valData() override {return needs_reorder ? val_data.data() + idx : this->loader.valData();}
};


} // end namespace BPCells