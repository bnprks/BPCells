#pragma once

#include <cstdint>
#include <stdexcept>
#include <vector>

namespace BPCells {

template<typename T>
struct SparseVector {
    uint32_t *idx;
    T *val;
    uint32_t capacity;
};


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
template<typename T>
class MatrixLoader {
public:
    virtual ~MatrixLoader() {};

    // Return the count of rows and columns
    virtual uint32_t rows() const = 0;
    virtual uint32_t cols() const = 0;

    // Reset the iterator to start from the beginning
    virtual void restart() = 0;

    // Advance to the next column, or return false if there
    // are no more columns
    virtual bool nextCol() = 0;

    // Return the index of the current column
    virtual uint32_t currentCol() const = 0;
    
    // Return number of matrix entries loaded. Should repeatedly return 0 
    // at the end of a column, until nextCol is called
    // Return -1 for error
    virtual int32_t load(uint32_t count, SparseVector<T> buffer) = 0;
};

template<typename T>
class MatrixLoaderWrapper : public MatrixLoader<T> {
protected:
    MatrixLoader<T> &loader;
    uint32_t current_col;
public:
    MatrixLoaderWrapper(MatrixLoader<T> &loader) : loader(loader), current_col(UINT32_MAX-1) {};

    uint32_t rows() const override {return loader.rows(); }
    uint32_t cols() const override {return loader.cols(); }

    void restart() override {loader.restart(); }

    bool nextCol() override {
        if (!this->loader.nextCol()) return false;
        current_col = this->loader.currentCol();
        return true;
    }
    
    uint32_t currentCol() const override {return loader.currentCol();}
};

template<typename T> 
class MatrixIterator : public MatrixLoaderWrapper<T> {
private:
    const uint32_t chunk_capacity;
    int32_t chunk_size;
    uint32_t idx;
    std::vector<uint32_t> row_buf;
    std::vector<T> val_buf;
    SparseVector<T> buf;
public:
    // Construct iterator with a given internal buffer size (must be a power of 2 >= 128)
    MatrixIterator(MatrixLoader<T> &loader, uint32_t buffer_size = 1024) : 
        MatrixLoaderWrapper<T>(loader),
        chunk_capacity(buffer_size), chunk_size(buffer_size), idx(buffer_size), 
        row_buf(buffer_size), val_buf(buffer_size) {

        if (buffer_size < 128)
            throw std::invalid_argument("buffer_size must be >= 128");
        else if (buffer_size & (buffer_size - 1)) 
            throw std::invalid_argument("buffer_size must be a power of 2");

        buf.idx = &row_buf[0];
        buf.val = &val_buf[0];
        buf.capacity = buffer_size;
    }

    // Reset the iterator to start from the beginning
    void restart() override {
        idx = chunk_capacity;
        this->loader.restart();
    }
    
    // Return false if there isn't another column to access
    bool nextCol() override {
        if (!this->loader.nextCol()) return false;
        this->current_col = this->loader.currentCol();
        return true;
    }

    // Return false if there isn't another entry in the current column
    inline bool nextValue() {
        idx += 1;
        if (idx >= chunk_size) {
            chunk_size = load(chunk_capacity, buf);
            idx = 0;
        }
        return chunk_size > 0;
    }
    // Access current row, column, and value
    inline uint32_t row() const {return buf.idx[idx]; };
    inline uint32_t col() const {return this->current_col; };
    inline T val() const {return buf.val[idx]; };

    int32_t load(uint32_t count, SparseVector<T> buffer) override {return this->loader.load(count, buffer);};
};


template<typename T> 
class MatrixTransposeIterator : public MatrixLoaderWrapper<T> {
private:
    const uint32_t chunk_capacity;
    int32_t chunk_size;
    uint32_t idx;
    std::vector<uint32_t> row_buf;
    std::vector<T> val_buf;
    SparseVector<T> buf;
public:
    // Construct iterator with a given internal buffer size (must be a power of 2 >= 128)
    MatrixTransposeIterator(MatrixLoader<T> &loader, uint32_t buffer_size = 1024) : 
        MatrixLoaderWrapper<T>(loader),
        chunk_capacity(buffer_size), chunk_size(buffer_size), idx(buffer_size), 
        row_buf(buffer_size), val_buf(buffer_size) {

        if (buffer_size < 128)
            throw std::invalid_argument("buffer_size must be >= 128");
        else if (buffer_size & (buffer_size - 1)) 
            throw std::invalid_argument("buffer_size must be a power of 2");

        buf.idx = &row_buf[0];
        buf.val = &val_buf[0];
        buf.capacity = buffer_size;
    }

    // Reset the iterator to start from the beginning
    void restart() override {
        idx = chunk_capacity;
        this->loader.restart();
    }
    
    // Return false if there isn't another column to access
    bool nextCol() override {
        if (!this->loader.nextCol()) return false;
        this->current_col = this->loader.currentCol();
        return true;
    }

    // Return false if there isn't another entry in the current column
    inline bool nextValue();
    // Access current row, column, and value
    inline uint32_t col() const {return buf.idx[idx]; };
    inline uint32_t row() const {return this->current_col; };
    inline T val() const {return buf.val[idx]; };

    int32_t load(uint32_t count, SparseVector<T> buffer) override {return this->loader.load(count, buffer);};
};

template <typename T>
inline bool MatrixTransposeIterator<T>::nextValue() {
    idx += 1;
    if (idx >= chunk_size) {
        chunk_size = load(chunk_capacity, buf);
        idx = 0;
    }
    return chunk_size > 0;
}

template<typename T> 
class MatrixWriter {
public:
    // Return false on failure, true on success
    virtual bool write(MatrixLoader<T> &mat, void (*checkInterrupt)(void) = NULL) = 0;
};

template<typename Tin, typename Tout>
class MatrixConverterLoader : public MatrixLoader<Tout> {
private:
    MatrixLoader<Tin> &loader;
    std::vector<Tin> vals;
public:
    MatrixConverterLoader(MatrixLoader<Tin> &loader) : 
        loader(loader) {}

    uint32_t rows() const override {return loader.rows(); }
    uint32_t cols() const override {return loader.cols(); }

    void restart() override {loader.restart(); }

    bool nextCol() override {return loader.nextCol(); }
    uint32_t currentCol() const override {return loader.currentCol();}

    int32_t load(uint32_t count, SparseVector<Tout> buffer) override {
        SparseVector<Tin> inner_buffer;
        if (buffer.capacity > vals.size()) vals.resize(buffer.capacity);
        inner_buffer.capacity = buffer.capacity;
        inner_buffer.idx = buffer.idx;
        inner_buffer.val = &vals[0];

        int32_t ret = loader.load(count, inner_buffer);

        for (int i = 0; i < ret; i++) {
            buffer.val[i] = inner_buffer.val[i];
        }
        return ret;
    };
};

} // end namespace BPCells