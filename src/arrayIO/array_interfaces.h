#pragma once
#include <stdexcept>

namespace BPCells {

class UIntReadBuffer {
public:
    const uint32_t *data;
    uint32_t size;
};

class UIntWriteBuffer {
public:
    uint32_t *data;
    uint32_t size;
    inline bool isempty() {return size == 0;}
};


class UIntReader {
protected:
    size_t total_size;
    UIntReadBuffer read_buffer = {NULL, 0};

    // Helper function for subclasses to implement. This will get called only
    // when there's not enough loaded data, and we need to load more while copying
    // over the unread data
    virtual void _ensureCapacity(size_t capacity) = 0;
public:
    virtual ~UIntReader() = default;
    // Get a pointer the input data buffer, where input data can be read
    inline const uint32_t* data() const {return read_buffer.data;}
    // Get the number of elements ready to read from data()
    inline uint32_t capacity() const {return read_buffer.size;}
    // Advance the data buffer pointer by count without reading more of the underlying
    // stream
    inline void advance(uint32_t count) {
        if(count > read_buffer.size) {
            printf("Called advance with count > read_buffer.size, count=%d, size=%d\n",count, read_buffer.size);
            throw std::runtime_error("Invalid advance amount");
        }
        read_buffer.data += count;
        read_buffer.size -= count;
    }
    // Read one element of the input stream, throwing an exception if there are
    // no more entries to read
    inline uint32_t read_one() {
        ensureCapacity(1);
        uint32_t val = *read_buffer.data;
        advance(1);
        return val;
    }

    // Read from input and update data() and capacity() pointers.
    // Return false if no data could be read
    virtual bool next() = 0;
    // Seek to a different position in the stream (first integer is position 0),
    // and read from input, updating data() and capacity() pointers
    // (Tip: Don't call next directly after calling seek, as it would throw out
    // the data read during the seek() call)
    // Return false if no data could be read from the seek location
    virtual bool seek(size_t pos) = 0;
    // Return the total number of integers in this reader
    size_t size() const {return total_size;}

    // Ensure that the ReadBuffer has at least `capacity` items it can read.
    // May result in additional data copying + reading if not enough items are currently loaded.
    // throw an error if there isn't enough left in the input to read
    inline void ensureCapacity(size_t capacity) {
        if (this->capacity() == 0 && capacity != 0) {
            if (!next()) throw std::runtime_error("Not enough remaining data to ensure read capacity");
        }
        if (this->capacity() >= capacity) return;
        _ensureCapacity(capacity);
    }
};

class UIntWriter {
protected:
    size_t total_size;
    UIntWriteBuffer write_buffer = {NULL, 0};

    // Helper function for subclasses to implement. This will get called only
    // when there's not enough data buffer, and we need space to write more
    virtual void _ensureCapacity(size_t capacity) = 0;
public:
    // The destructor should flush any remaining data
    virtual ~UIntWriter() = default;

    // Get a pointer the output data buffer
    inline uint32_t* data() const {return write_buffer.data;}
    // Get the amount of elements which can be written to data()
    inline uint32_t capacity() const {return write_buffer.size;}

    // Advance the data buffer pointer by count without writing any of the
    // data buffer
    inline void advance(uint32_t count) {
        if(count > write_buffer.size) {
            printf("Called advance with count > write_buffer.size, count=%d, size=%d\n",count, write_buffer.size);
            throw std::runtime_error("Invalid advance amount");
        }
        write_buffer.size -= count;
        write_buffer.data += count;
    }
    inline void write_one(uint32_t val) {
        ensureCapacity(1);
        *write_buffer.data = val;
        advance(1);
    }

    // Write data as needed, and update data() and capacity() pointers.
    // buffer is only available for writing until the following next() or backup() call.
    virtual void next() = 0;

    // Instruct the writer to "back up", ignoring the last `count`
    // integers from the data buffer next time it flushes to disk
    void backup(size_t count) {
        write_buffer.size -= count;
    }

    // Ensure that the WriteBuffer has at least `capacity` items it can still write.
    // May result in additional data output if not enough items are currently loaded.
    // Note: Any data written beyond the current location of data() will be lost.
    inline void ensureCapacity(size_t capacity) {
        if (this->capacity() >= capacity) return;
        if (this->capacity() == 0) next();
        if (this->capacity() < capacity) _ensureCapacity(capacity);
    }

    
};


}// end namespace BPCells