#include "PackedMatrix.h"

namespace BPCells {

PackedMatrix::PackedMatrix(ReaderPtr &&val_data, ReaderPtr &&val_idx, ReaderPtr &&row_data, 
                ReaderPtr &&row_starts, ReaderPtr &&row_idx, ReaderPtr &&col_ptr, ReaderPtr &&row_count) :
            val_data(std::move(val_data)), val_idx(std::move(val_idx)), row_data(std::move(row_data)), row_starts(std::move(row_starts)), 
            row_idx(std::move(row_idx)), col_ptr(std::move(col_ptr)), row_count(std::move(row_count)) {
        
    n_rows = this->row_count->read_one();
    n_cols = this->col_ptr->size() - 1;
    next_col_ptr = this->col_ptr->read_one();
}

void PackedMatrix::load128(uint32_t *row_out, uint32_t *val_out) {
        uint32_t next_val_idx = val_idx->read_one(), 
                 next_row_idx = row_idx->read_one(),
                 row_start = row_starts->read_one();
    // Unpack start
    uint32_t bits_row = (next_row_idx - prev_row_idx) / 4;
    row_data->ensureCapacity(bits_row * 4);
    simdunpackd1z(row_start, row_data->data(), row_out, bits_row);
    row_data->advance(bits_row * 4);
    prev_row_idx = next_row_idx;

    // Unpack val
    uint32_t bits_val = (next_val_idx - prev_val_idx) / 4;
    val_data->ensureCapacity(bits_val * 4);
    simdunpackFOR(1, val_data->data(), val_out, bits_val);
    val_data->advance(bits_val * 4);
    prev_val_idx = next_val_idx;
}

uint32_t PackedMatrix::rows() const {return n_rows;}
uint32_t PackedMatrix::cols() const {return n_cols;}

void PackedMatrix::restart() {
    current_col=UINT32_MAX; 
    current_idx=UINT32_MAX; 
    col_ptr->seek(0);
    next_col_ptr = col_ptr->read_one();
    // The other buffers will get reset properly as soon 
    // as nextCol is called
}

bool PackedMatrix::nextCol() {
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
        val_idx->seek(current_idx / 128);
        row_idx->seek(current_idx / 128);
        row_starts->seek(current_idx / 128);

        prev_val_idx = val_idx->read_one();
        prev_row_idx = row_idx->read_one();

        val_data->seek(prev_val_idx);
        row_data->seek(prev_row_idx);
    }
    
    next_col_ptr = col_ptr->read_one();
    return true;
}

uint32_t PackedMatrix::currentCol() const {return current_col;}

int32_t PackedMatrix::load(uint32_t count, SparseVector<uint32_t> buffer) {
    if (current_idx >= next_col_ptr) {
        return 0;
    }
    if (count < 128) throw std::runtime_error("Must load >128 entries at a time from PackedMatrix");    
    
    uint32_t i = 0;
    // Check if we are starting a column and need to copy from the buffer
    if (current_idx % 128 != 0) {
        uint32_t copy_num = std::min(128 - current_idx % 128, next_col_ptr - current_idx);
        std::memmove(buffer.idx, &row_buf[current_idx % 128], copy_num * sizeof(uint32_t));
        std::memmove(buffer.val, &val_buf[current_idx % 128], copy_num * sizeof(uint32_t));
        i += copy_num;
        current_idx += copy_num;
    }
    
    for (; i + 128 <= count && next_col_ptr - current_idx >= 128; i+= 128) {
        load128(&buffer.idx[i], &buffer.val[i]);
        current_idx += 128;
    }

    if (i < count && current_idx < next_col_ptr) {
        // Save any straglers into the buffer so we can read the rest of the
        // block of 128 on the subsequent load call
        uint32_t remaining = std::min(count-i, next_col_ptr-current_idx);
        load128(row_buf, val_buf);
        std::memmove(&buffer.idx[i], row_buf, sizeof(uint32_t) * remaining);
        std::memmove(&buffer.val[i], val_buf, sizeof(uint32_t) * remaining);
        i += remaining;
        current_idx += remaining;
    }

    return i;
}

PackedMatrixWriter::PackedMatrixWriter(WriterPtr &&val_data, WriterPtr &&val_idx, WriterPtr &&row_data, 
                WriterPtr &&row_starts, WriterPtr &&row_idx, WriterPtr &&col_ptr, WriterPtr &&row_count) :
            val_data(std::move(val_data)), val_idx(std::move(val_idx)), row_data(std::move(row_data)), row_starts(std::move(row_starts)), 
            row_idx(std::move(row_idx)), col_ptr(std::move(col_ptr)), row_count(std::move(row_count)) {}
    
void PackedMatrixWriter::pack128(const uint32_t* idx_in, const uint32_t* val_in,
                uint32_t &cur_val_idx, uint32_t &cur_row_idx) {
    // Note: assumes row_data and val_data have sufficient capacity

    // Pack row
    uint32_t row_bits = simdmaxbitsd1z(idx_in[0], idx_in);
    row_starts->write_one(idx_in[0]);
    simdpackd1z(idx_in[0], idx_in, row_data->data(), row_bits);
    row_data->advance(row_bits*4);
    cur_row_idx += row_bits*4;
    row_idx->write_one(cur_row_idx);

    // Pack val
    uint32_t val_bits = simdmaxbitsFOR(1, val_in);
    simdpackFOR(1, val_in, val_data->data(), val_bits);
    val_data->advance(val_bits*4);
    cur_val_idx += val_bits*4;
    val_idx->write_one(cur_val_idx);
}

bool PackedMatrixWriter::write(MatrixLoader<uint32_t> &mat, void (*checkInterrupt)(void)) {
    uint32_t cur_val_idx = 0;
    uint32_t cur_row_idx = 0;

    uint32_t col = 0;
    uint32_t idx = 0; // Index of for col_ptr array
    uint32_t loaded = 0; // Number of values loaded most recently

    uint32_t val_buffer[128], row_buffer[128]; 

    col_ptr->write_one(idx);
    val_idx->write_one(cur_val_idx);
    row_idx->write_one(cur_row_idx);
    
    while (mat.nextCol()) {
        if (checkInterrupt != NULL) checkInterrupt();
        if (mat.currentCol() < col)
            throw std::runtime_error("PackedMatrixWriter encountered out-of-order columns");
        while (col < mat.currentCol()) {
            col_ptr->write_one(idx);
            col += 1;
        }
        while(true) {
            // Strategy: Load up to 256 entries into row_buffer and val_buffer. 
            // If we end up getting a non-multiple of 128 values to compress, then try to make our next
            // load continue where the last one left off.
            // However, if we get to a point where there is space for less than 128 more values,
            // then copy the leftovers to the beginning of row_buffer and val_buffer and continue.


            // Writing strategy: Use the output data buffer for val_data and row_data
            // in order to load >= 256 matrix entries prior to actually compressing + writing them.
            // We will need our own buffer space for 128 values to store stragglers
            size_t capacity = std::min(val_data->capacity(), row_data->capacity());
            
            uint32_t *row_data_buf = row_data->data();
            uint32_t *val_data_buf = val_data->data();
            uint32_t i;
            for (i = 0; i + 128 <= capacity; i += 128) {
                pack128(row_data_buf + i, val_data_buf + i, cur_val_idx, cur_row_idx);
            }
            uint32_t hangover = capacity % 128;
            // If our remaining capacity is not a multiple of 128, then stash the data before calling ensureCapacity
            if (hangover != 0) {
                std::memmove(row_buffer, row_data_buf + i, sizeof(uint32_t)*hangover);
                std::memmove(val_buffer, val_data_buf + i, sizeof(uint32_t)*hangover);
            }

            val_data->ensureCapacity(256);
            row_data->ensureCapacity(256);
            capacity = std::min(val_data->capacity(), row_data->capacity());

            if (hangover != 0) {
                std::memmove(row_data->data(), row_buffer, sizeof(uint32_t)*hangover);
                std::memmove(val_data->data(), val_buffer, sizeof(uint32_t)*hangover);
            }

            loaded = mat.load(capacity-hangover, {row_data->data()+hangover, val_data->data()+hangover, (uint32_t) capacity-hangover});
            idx += loaded;

            val_data->backup(val_data->capacity() - loaded - hangover);
            row_data->backup(row_data->capacity() - loaded - hangover);
            if (loaded == 0) break;
        } // End of reading column
    }
    // Finish compression by taking care of any leftovers
    size_t capacity = std::min(val_data->capacity(), row_data->capacity()); // should == hangover from above
    if (capacity != 0) {
        std::memmove(row_buffer, row_data->data(), sizeof(uint32_t)*capacity);
        std::memmove(val_buffer, val_data->data(), sizeof(uint32_t)*capacity);
        for (uint32_t i = capacity; i % 128 != 0; i++) {
            val_buffer[i] = val_buffer[i-1];
            row_buffer[i] = row_buffer[i-1];
        }
        row_data->ensureCapacity(128);
        val_data->ensureCapacity(128);
        pack128(row_buffer, val_buffer, cur_val_idx, cur_row_idx);
    }
    row_count->write_one(mat.rows());

    col_ptr->write_one(idx);

    // Flush any data that we've advanced past already
    val_data->backup(val_data->capacity()); val_data.reset(nullptr);
    val_idx->backup(val_idx->capacity()); val_idx.reset(nullptr);
    row_data->backup(row_data->capacity()); row_data.reset(nullptr);
    row_starts->backup(row_starts->capacity()); row_starts.reset(nullptr);
    row_idx->backup(row_idx->capacity()); row_idx.reset(nullptr);
    col_ptr->backup(col_ptr->capacity()); col_ptr.reset(nullptr);
    row_count->backup(row_count->capacity()); row_count.reset(nullptr);

    return true;
};

} // end namespace BPCells