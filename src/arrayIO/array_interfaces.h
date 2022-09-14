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

class NullStringWriter : public StringWriter {
public:
    void write(const StringReader &reader) override {}
};

template<class T>
class BulkNumReader {
public:
    virtual ~BulkNumReader() = default;

    // Return total number of integers in the reader
    virtual uint32_t size() const = 0;

    // Change the next load to start at index pos
    virtual void seek(uint32_t pos) = 0;

    // Copy up to `count` integers into `out`, returning the actual number copied.
    // Will always load >0 if count is >0
    // Note: It is the caller's responsibility to ensure there is no data overflow, i.e.
    // that a load does not try to read past size() total elements
    virtual uint32_t load(T *out, uint32_t count) = 0;
};


template<class From, class To>
class BulkNumReaderConverter : public BulkNumReader<To> {
private:
    std::unique_ptr<BulkNumReader<From>> reader;
    std::vector<From> buffer;
public:
    BulkNumReaderConverter(std::unique_ptr<BulkNumReader<From>> &&reader) :
        reader(std::move(reader)) {}
    uint32_t size() const override {return reader->size();}
    void seek(uint32_t pos) override {return reader->seek(pos);}
    uint32_t load(To *out, uint32_t count) override {
        if (buffer.size() < count) buffer.resize(count);
        uint32_t loaded = reader->load(buffer.data(), count);
        for (uint32_t i = 0; i < loaded; i++) {
            out[i] = (To) buffer[i];
        }
        return loaded;
    }
};

using UIntBulkReader = BulkNumReader<uint32_t>;

template<class T>
class BulkNumWriter {
public:
    virtual ~BulkNumWriter() = default;

    // Write up to `count` integers from `in`, returning the actual number written.
    // Will always write >0, otherwise throwing an exception for an error
    // Note: The writer is allowed to modify the input data, so it might be 
    // modified after calling write()
    virtual uint32_t write(T *in, uint32_t count) = 0;
    
    // Flush any remaining data to disk. No more write calls will be made after
    // calling fiinalize. By default a no-op.
    virtual void finalize() {}
};

template<class From, class To>
class BulkNumWriterConverter : public BulkNumWriter<To> {
private:
    std::unique_ptr<BulkNumWriter<From>> writer;
    std::vector<From> buffer;
public:
    BulkNumWriterConverter(std::unique_ptr<BulkNumWriter<From>> writer) :
        writer(std::move(writer)) {}
    
    uint32_t write(To *in, uint32_t count) override {
        if (count > buffer.size()) buffer.resize(count);
        for (uint32_t i = 0; i < count; i++) {
            buffer[i] = (From) in[i];
        }
        return writer->write(buffer.data(), count);
    }

    void finalize() override {writer->finalize();}
};

using UIntBulkWriter = BulkNumWriter<uint32_t>;


// NumReader -- Read a stream of numbers, providing a conveinient
// interface to a BulkReader.
// Usage:
// 1. At construction time, specify the input data source, how much data should be
//    buffered internally, and the maximum data to provide in 1 chunk to downstream users.
// 2. Call ensureCapacity() or requestCapacity() to load data
// 3. Loaded data can be accessed starting at the pointer returned by data(), continuing for
//    capacity() elements. If the data is being used piece-by-piece, call advance() to mark
//    the first few elements as consumed.
// Note: Users of the class are free to modify the memory between data() and capacity() as needed.
// This is to support zero-copy transformers, where the loaded data is transformed or filtered in-place.
template<class T>
class NumReader {
protected:
    std::vector<T> buffer;
    uint32_t idx = 0; // Index of read data in buffer
    uint32_t available = 0; // Amount of read data in buffer
    uint32_t loaded = 0; // Amount of data loaded currently in buffer (note idx + available <= loaded as an invariant)
    uint32_t pos = 0; // Position of next integer that would be read

    std::unique_ptr<BulkNumReader<T>> reader;
    uint32_t total_size; // total size of reader

    uint32_t read_size; // Amount to provide users by default

public:
    NumReader() = default;
    NumReader(std::unique_ptr<BulkNumReader<T>> &&reader, uint32_t buffer_size, uint32_t read_size=1024) :
        buffer(buffer_size), reader(std::move(reader)), total_size(this->reader->size()), read_size(read_size) {}

    // Convert to a different output type, after which the current reader should
    // not be used
    template <class To>
    NumReader<To> convert() {
        return NumReader<To>(
            std::make_unique<BulkNumReaderConverter<T, To>>(std::move(reader)),
            buffer.size()
        );
    }

    // Pointer to data in buffer start
    inline T* data() {return buffer.data() + idx;}
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
            std::memmove(buffer.data(), data(), (loaded - idx)*sizeof(T));
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
        // Note: previous code allowed for re-using the existing loaded data buffer 
        // upon seek, but this could cause bugs when downstream code had already modified
        // the data buffer but expected clean data after seek+load. Look at file 
        // history to see the old version
        new_pos = std::min(new_pos, total_size);
        reader->seek(new_pos);
        pos = new_pos;
        loaded = 0;
        idx = 0;
        available = 0;  
    }
    
    // Read one element of the input stream, throwing an exception if there are
    // no more entries to read
    inline T read_one() {
        ensureCapacity(1);
        T val = *data();
        advance(1);
        return val;
    }

    // Advance the data buffer by count without reading more of the underlying stream
    inline void advance(uint32_t count) {
        idx += count;
        available -= count;
    }
};

using UIntReader = NumReader<uint32_t>;
using ULongReader = NumReader<uint64_t>;
using FloatReader = NumReader<float>;
using DoubleReader = NumReader<double>;

// Writer -- Write a stream of numbers, providing a conveinient
// interface to a BulkWriter.
// Usage:
// 1. At construction time, specify the input data source, how much data should be
//    buffered internally, and the maximum data to provide in 1 chunk to downstream users.
// 2. Call ensureCapacity() to request write-buffer space, flushing existing data in the buffer as needed
// 3. Data for output can be written starting at the pointer returned by data() continuing for 
//    capacity() elements. If the data is being written piece-by-piece, call advance() to mark
//    the first few elements as consumed.
// 4. Call finalize() after all data is done writing to flush any remainig data in the internal buffer
template <class T>
class NumWriter {
protected:
    std::vector<T> buffer;
    uint32_t idx = 0; // Index of buffer for next write
    std::unique_ptr<BulkNumWriter<T>> writer;
private:
    // Write all data up to idx+available, and copy any stragglers to the beginning
    inline void flush() {
        uint32_t written = writer->write(buffer.data(), idx);
        if (written == 0) 
            throw std::runtime_error("No data written after write request");
        idx = idx - written;
        if (idx != 0) {
            // Copy leftover data to the front of the buffer
            std::memmove(buffer.data(), buffer.data() + written, sizeof(T) * idx);
        }
    }

public:
    NumWriter(std::unique_ptr<BulkNumWriter<T>> &&writer, uint32_t buffer_size) :
        buffer(buffer_size), writer(std::move(writer)) {}
    NumWriter() = default;
    NumWriter(NumWriter<T> &&other) = default;
    NumWriter<T>& operator=(NumWriter<T>&& other) = default;
    
    // Convert to a different output type, after which the current writer should
    // not be used
    template <class To>
    NumWriter<To> convert() {
        return NumWriter<To>(
            std::make_unique<BulkNumWriterConverter<T, To>>(std::move(writer)),
            buffer.size()
        );
    }

    // Pointer to data in buffer start
    inline T* data() {return buffer.data() + idx;}
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

    inline uint32_t maxCapacity() const {return buffer.size();}
    
    // Read one element of the input stream, throwing an exception if there are
    // no more entries to read
    inline void write_one(T val) {
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
        writer.reset();
    }
};

using UIntWriter = NumWriter<uint32_t>;
using ULongWriter = NumWriter<uint64_t>;
using FloatWriter = NumWriter<float>;
using DoubleWriter = NumWriter<double>;

class WriterBuilder {
public:
    virtual ~WriterBuilder() = default;
    virtual UIntWriter createUIntWriter(std::string name) = 0;
    virtual ULongWriter createULongWriter(std::string name) = 0;
    virtual FloatWriter createFloatWriter(std::string name) = 0;
    virtual DoubleWriter createDoubleWriter(std::string name) = 0;

    template<class T>
    NumWriter<T> create(std::string name) {
        if constexpr (std::is_same_v<T, uint32_t>) {return createUIntWriter(name);}
        if constexpr (std::is_same_v<T, uint64_t>) {return createULongWriter(name);}
        if constexpr (std::is_same_v<T, float>) {return createFloatWriter(name);}
        if constexpr (std::is_same_v<T, double>) {return createDoubleWriter(name);}
    }

    virtual std::unique_ptr<StringWriter> createStringWriter(std::string name) = 0;
    virtual void writeVersion(std::string version) = 0;

    // Delete all writers with the given name, throwing an exception if this
    // implementation does not support deletion
    virtual void deleteWriter(std::string name) = 0;
};

class ReaderBuilder {
public:
    virtual ~ReaderBuilder() = default;
    virtual UIntReader openUIntReader(std::string name) = 0;
    virtual ULongReader openULongReader(std::string name) = 0;
    virtual FloatReader openFloatReader(std::string name) = 0;
    virtual DoubleReader openDoubleReader(std::string name) = 0;

    template<class T>
    NumReader<T> open(std::string name) {
        if constexpr (std::is_same_v<T, uint32_t>) {return openUIntReader(name);}
        if constexpr (std::is_same_v<T, uint64_t>) {return openULongReader(name);}
        if constexpr (std::is_same_v<T, float>) {return openFloatReader(name);}
        if constexpr (std::is_same_v<T, double>) {return openDoubleReader(name);}
    }

    virtual std::unique_ptr<StringReader> openStringReader(std::string name) = 0;
    virtual std::string readVersion() = 0;
};

}// end namespace BPCells