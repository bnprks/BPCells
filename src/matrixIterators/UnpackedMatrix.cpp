#include "UnpackedMatrix.h"

namespace BPCells {

UnpackedMatrix::UnpackedMatrix(ReaderPtr &&val, ReaderPtr &&row,
                ReaderPtr &&col_ptr, ReaderPtr &&row_count):
            val(std::move(val)), row(std::move(row)), col_ptr(std::move(col_ptr)), row_count(std::move(row_count)) {
    n_rows = this->row_count->read_one();
    n_cols = this->col_ptr->size() - 1;
    next_col_ptr = this->col_ptr->read_one();
}

uint32_t UnpackedMatrix::rows() const {return n_rows;}
uint32_t UnpackedMatrix::cols() const {return n_cols;}

void UnpackedMatrix::restart() {
    current_col=UINT32_MAX; 
    current_idx=UINT32_MAX; 
    col_ptr->seek(0);
    next_col_ptr = col_ptr->read_one();
    // The other buffers will get reset properly as soon 
    // as nextCol is called
}

bool UnpackedMatrix::nextCol() {
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
        val->seek(current_idx);
        row->seek(current_idx);
    }
    
    next_col_ptr = col_ptr->read_one();
    return true;
}

uint32_t UnpackedMatrix::currentCol() const {return current_col;}

int32_t UnpackedMatrix::load(uint32_t count, SparseVector<uint32_t> buffer) {
    if (current_idx >= next_col_ptr) {
        return 0;
    }

    // Load row / val if necessary
    if (row->capacity() == 0 && !row->next()) {throw std::runtime_error("Failed to load next row data");}
    if (val->capacity() == 0 && !val->next()) throw std::runtime_error("Failed to load next val data");
    
    uint32_t capacity = std::min(row->capacity(), val->capacity());
    uint32_t load = std::min(next_col_ptr - current_idx, count);
    load = std::min(load, capacity);
    current_idx += load;

    std::memmove(buffer.idx, row->data(), load*sizeof(uint32_t));
    std::memmove(buffer.val, val->data(), load*sizeof(uint32_t));
    
    row->advance(load);
    val->advance(load);

    return load;
}

UnpackedMatrixWriter::UnpackedMatrixWriter(WriterPtr &&val, WriterPtr &&row, 
                WriterPtr &&col_ptr, WriterPtr &&row_count) :
                val(std::move(val)), row(std::move(row)), col_ptr(std::move(col_ptr)), row_count(std::move(row_count)) {}

bool UnpackedMatrixWriter::write(MatrixLoader<uint32_t> &mat, void (*checkInterrupt)(void)) {
    uint32_t col = 0;
    uint32_t idx = 0; // Index of for col_ptr array
    uint32_t loaded = 0; // Number of values loaded most recently

    col_ptr->write_one(idx);

    while (mat.nextCol()) {
        if (checkInterrupt != NULL) checkInterrupt();
        if (mat.currentCol() < col)
            throw std::runtime_error("UnpackedMatrixWriter encountered out-of-order columns");
        while (col < mat.currentCol()) {
            col_ptr->write_one(idx);
            col += 1;
        }
        while(true) {
            // This is a bit of a workaround, since the packed loaders require load requests
            // for at least 128 values at a time
            row->ensureCapacity(128);
            val->ensureCapacity(128);
            
            size_t capacity = std::min(val->capacity(), row->capacity());
            loaded = mat.load(capacity, {row->data(), val->data(), (uint32_t) capacity});
            row->advance(loaded);
            val->advance(loaded);
            idx += loaded;
            if (loaded == 0) break;
        } // End of reading column
    }
    row_count->write_one(mat.rows());
    col_ptr->write_one(idx);

    // Flush any data that we've advanced past already
    row->backup(row->capacity()); row.reset(nullptr);
    val->backup(val->capacity()); val.reset(nullptr);
    col_ptr->backup(col_ptr->capacity()); col_ptr.reset(nullptr);
    row_count->backup(row_count->capacity()); row_count.reset(nullptr);

    return true;
};

} // end namespace BPCells