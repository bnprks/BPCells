#pragma once
#include <cstdint>
#include <cstddef>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace BPCells {


class StringReader {
public:
    virtual ~StringReader() = default;
    virtual const char* get(uint32_t idx) const = 0;
    virtual uint32_t size() const = 0;
};

class StringWriter {
public:
    virtual ~StringWriter() = default;
    virtual void write(const StringReader &reader) = 0;
};

// Simple generic StringReader designed to allow for transparent reading
// of std::vector<std::string> 
class VecStringReader : public StringReader {
private:
    std::vector<std::string> data;
public:
    VecStringReader(std::vector<std::string> data);
    const char* get(uint32_t idx) const override;
    uint32_t size() const override;
};


class UIntBulkReader {
public:
    virtual ~UIntBulkReader() = default;

    // Return total number of integers in the reader
    virtual uint32_t size() const = 0;

    // Change the next load to start at index pos
    virtual void seek(uint32_t pos) = 0;

    // Copy up to `count` integers into `out`, returning the actual number copied.
    // Will always load >0 if count is >0
    // Note: It is the caller's responsibility to ensure there is no data overflow, i.e.
    // that a load does not try to read past size() total elements
    virtual uint32_t load(uint32_t *out, uint32_t count) = 0;
};

class UIntBulkWriter {
public:
    virtual ~UIntBulkWriter() = default;

    // Write up to `count` integers from `in`, returning the actual number written.
    // Will always write >0, otherwise throwing an exception for an error
    // Note: The writer is allowed to modify the input data, so it might be 
    // modified after calling write()
    virtual uint32_t write(uint32_t *in, uint32_t count) = 0;
    
    // Flush any remaining data to disk. No more write calls will be made after
    // calling fiinalize. By default a no-op.
    virtual void finalize() {}
};


// UIntReader -- Read a stream 32-bit unsigned integers, providing a conveinient
// interface to a UIntBulkReader.
// Usage:
// 1. At construction time, specify the input data source, how much data should be
//    buffered internally, and the maximum data to provide in 1 chunk to downstream users.
// 2. Call ensureCapacity() or requestCapacity() to load data
// 3. Loaded data can be accessed starting at the pointer returned by data(), continuing for
//    capacity() elements. If the data is being used piece-by-piece, call advance() to mark
//    the first few elements as consumed.
// Note: Users of the class are free to modify the memory between data() and capacity() as needed.
// This is to support zero-copy transformers, where the loaded data is transformed or filtered in-place.
class UIntReader {
protected:
    std::vector<uint32_t> buffer;
    uint32_t idx = 0; // Index of read data in buffer
    uint32_t available = 0; // Amount of read data in buffer
    uint32_t loaded = 0; // Amount of data loaded currently in buffer (note idx + available <= loaded as an invariant)
    uint32_t pos = 0; // Position of next integer that would be read

    std::unique_ptr<UIntBulkReader> reader;
    uint32_t total_size; // total size of reader

    uint32_t read_size; // Amount to provide users by default

public:
    UIntReader(std::unique_ptr<UIntBulkReader> &&reader, uint32_t buffer_size, uint32_t read_size=1024);

    // Pointer to data in buffer start
    inline uint32_t* data() {return buffer.data() + idx;}
    // Number of available entries in data() buffer
    inline uint32_t capacity() const {return available;};

    // Try to ensure there are at least `new_capacity` items available to read 
    // (i.e. capacity() > new_capacity). Return false if there was not enough
    // data left to fill out the requested capacity.
    // Any remaining data between data() and data()+capacity() will be readable
    // at the value of data() after the call.
    inline bool requestCapacity(uint32_t new_capacity) {
        if (new_capacity > read_size) return false;

        if (loaded - idx >= new_capacity) {
            // We already have the data loaded, so just expand capacity as required
            // We want -- if available >= new_capacity no change
            // otherwise, provide a bit extra with read_size
            available = available >= new_capacity ? available : std::min(read_size, loaded-idx);
            return true;
        } 

        if (idx != 0) {
            std::memmove(buffer.data(), data(), (loaded - idx)*sizeof(uint32_t));
            loaded = loaded - idx;
            idx = 0;
        }

        while (loaded < read_size) {
            uint32_t load_size = std::min((uint32_t) buffer.size() - loaded, total_size-pos);
            if (load_size == 0) break;

            uint32_t newly_loaded = reader->load(buffer.data() + loaded, load_size);
            loaded += newly_loaded;
            pos += newly_loaded;
        }
        available = std::min(loaded, read_size);
        return available >= new_capacity;
    }

    // Request a non-zero amount of read capacity, returning false if there is no more
    // data to load
    inline bool requestCapacity() {return requestCapacity(1);}

    // A variant of requestCapacity that throws an error if the requested capacity
    // could not be loaded.
    inline void ensureCapacity(uint32_t new_capacity) {
        if (requestCapacity(new_capacity)) return;
        if (new_capacity > read_size)
            throw std::invalid_argument("new_capacity can't be larger than load_size");
        throw std::runtime_error("Not enough remaining data to ensure read capacity");
    }

    // Total number of ints in the reader
    inline uint32_t size() const {return total_size;}

    // Seek to a different position in the stream (first integer is position 0),
    // resetting data() and capacity() pointers to have 0 capacity.
    inline void seek(uint32_t new_pos) {
        uint32_t buf_start_pos = pos - loaded;

        if (buf_start_pos <= new_pos && pos > new_pos) {
            // Don't seek, just let the next read come from the existing buffer
            idx = new_pos - buf_start_pos;
            available = 0;
        } else {
            new_pos = std::min(new_pos, total_size);
            reader->seek(new_pos);
            pos = new_pos;
            loaded = 0;
            idx = 0;
            available = 0;
        }
    }
    
    // Read one element of the input stream, throwing an exception if there are
    // no more entries to read
    inline uint32_t read_one() {
        ensureCapacity(1);
        uint32_t val = *data();
        advance(1);
        return val;
    }

    // Advance the data buffer by count without reading more of the underlying stream
    inline void advance(uint32_t count) {
        idx += count;
        available -= count;
    }
};


// Writer -- Write a stream 32-bit unsigned integers, providing a conveinient
// interface to a UIntBulkWriter.
// Usage:
// 1. At construction time, specify the input data source, how much data should be
//    buffered internally, and the maximum data to provide in 1 chunk to downstream users.
// 2. Call ensureCapacity() to request write-buffer space, flushing existing data in the buffer as needed
// 3. Data for output can be written starting at the pointer returned by data() continuing for 
//    capacity() elements. If the data is being written piece-by-piece, call advance() to mark
//    the first few elements as consumed.
// 4. Call finalize() after all data is done writing to flush any remainig data in the internal buffer
class UIntWriter {
protected:
    std::vector<uint32_t> buffer;
    uint32_t idx = 0; // Index of buffer for next write
    std::unique_ptr<UIntBulkWriter> writer;

private:
    // Write all data up to idx+available, and copy any stragglers to the beginning
    inline void flush() {
        uint32_t written = writer->write(buffer.data(), idx);
        if (written == 0) 
            throw std::runtime_error("No data written after write request");
        idx = idx - written;
        if (idx != 0) {
            // Copy leftover data to the front of the buffer
            std::memmove(buffer.data(), buffer.data() + written, sizeof(uint32_t) * idx);
        }
    }

public:
    UIntWriter(std::unique_ptr<UIntBulkWriter> &&writer, uint32_t buffer_size);
    UIntWriter(UIntWriter &&other) = default;
    
    // Pointer to data in buffer start
    inline uint32_t* data() {return buffer.data() + idx;}
    // Number of available entries in data() buffer
    inline uint32_t capacity() const {return buffer.size() - idx;};

    // Ensure that there are at least `capacity` items available to write.
    // throw an error if the requested capacity is larger than the write_size set at construction time.
    // May result in flushing data that has been passed using `advance`, but will not
    // flush any data from the existing entries betweeen `data()` and `data() + capacity()`
    inline void ensureCapacity(uint32_t new_capacity) {
        if (new_capacity > buffer.size())
            throw std::invalid_argument("new_capacity can't be larger than load_size");
        
        while (buffer.size() - idx < new_capacity) {
            flush();
        }
    }

    inline void ensureCapacity() {ensureCapacity(1);}
    
    // Read one element of the input stream, throwing an exception if there are
    // no more entries to read
    inline void write_one(uint32_t val) {
        ensureCapacity(1);
        *data() = val;
        advance(1);
    }

    // Advance the data buffer by count without reading more of the underlying stream
    inline void advance(uint32_t count) {
        idx += count;
    }

    inline void finalize() {
        while (idx != 0) flush();
        writer->finalize();
    }
};


class WriterBuilder {
public:
    virtual ~WriterBuilder() = default;
    virtual UIntWriter createUIntWriter(std::string name) = 0;
    virtual std::unique_ptr<StringWriter> createStringWriter(std::string name) = 0;
    virtual void writeVersion(std::string version) = 0;
};

class ReaderBuilder {
public:
    virtual ~ReaderBuilder() = default;
    virtual UIntReader openUIntReader(std::string name) = 0;
    virtual std::unique_ptr<StringReader> openStringReader(std::string name) = 0;
    virtual std::string readVersion() = 0;
};

}// end namespace BPCells