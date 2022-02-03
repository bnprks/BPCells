#include "../bitpacking/bp128.h"
#include <cstring>

namespace BPCells {

// Note: For now I'm not using these, since it's seeming easier
// to just write a dedicated PackedFragments and PackedMatrix class
// that are able to result in some slightly bigger efficiencies.

template<class UIntWriter>
class BP128UIntWriter {
private:
    uint32_t [128]in_buf;
    size_t total_written = 0;
    size_t total_compressed = 0;

    UIntWriter data, idx, len;

    inline void write128(uint32_t *compress_buf) {
        uint32_t bits = simdmaxbits(in_buf);
        simdpack(in_buf, compress_buf, bits);
        data.write(compress_buf, bits * 4);
        idx.write(&total_compresed, 1);
        total_compressed += bits*4;
    }
public:
    BP128UIntWriter(UIntWriter&& data, UIntWriter&& idx, UIntWriter&& len) : data(data), idx(idx), len(len) {}
    BP128UIntWriter(const BP128UIntWriter&) = delete;
    BP128UIntWriter& operator= (const BP128UIntWriter&) = delete;

    // Append integers to output stream; Throws exception on failure
    void write(const uint32_t *buffer, uint32_t count) {
        uint32_t [128]compress_buf;
        uint32_t written = 0;

        while (written < count) {
            uint32_t items = std::min(128 - total_written%128, count-written);
            std::memmove(&in_buf[total_written % 128], &buffer[written], items*sizeof(uint32_t));
            written += items;
            total_written += items;
            if (total_written % 128 == 0) write128(compress_buf);
        }
    }

    void finalize() {
        len.write(&total_written, 1);
        while(total_written % 128 != 0) {
            in_buf[total_written % 128] = in_buf[(total_written % 128) - 1];
            total_written += 1;
        }
        write128();
    }
};

template<class UIntReader>
class BP128UIntReader {
private:
    uint32_t [128]out_buf;
    UIntReader data, idx, len;
    size_t pos = 0;
    uint32_t idx_prev;

    inline void read128(uint32_t *data_buf) {
        uint32_t idx_new;
        if (idx.read(&idx_new, 1) != 1) 
            throw std::runtime_error("Failed to read next idx value");
        if (data.read(data_buf, idx_new-idx_prev) != idx_new-idx_prev)
            throw std::runtime_error("Failed to read compressed data chunk");
        
        uint32_t bits = (idx_new - idx_prev)/4;
        simdunpack(data_buf, out_buf, bits);
        
        idx_prev = idx_new;
    }
public:
    BP128UIntReader() = default;
    BP128UIntReader(UIntReader&& data, UIntReader&& idx, UIntReader&& len) : data(data), idx(idx), len(len) {
        if (idx.read(&idx_prev, 1) != 1) 
            throw std::runtime_error("Failed to read next idx value");
    }
    
    BP128UIntReader(const BP128UIntReader&) = delete;
    BP128UIntReader& operator= (const BP128UIntReader&) = delete;

    uint32_t read(uint32_t *buffer, uint32_t count) {
        uint32_t [128]data_buf;
        uint32_t read = 0;
        while (read < count) {
            if (pos % 128 == 0) read128(data_buf);
            uint32_t items = std::min(128 - pos%128, count-read);
            std::memmove(buffer, &out_buf[pos%128], items*sizeof(uint32_t));
            read += items;
            pos += items;
        }
        return read;
    }
    uint32_t size() {
        uint32_t s;
        len.read(&s, 1);
        return s;
    }
    void seek(const std::size_t pos) {
        this->pos = pos;
        idx.seek(pos / 128);
        if (idx.read(&idx_prev, 1) != 1) 
            throw std::runtime_error("Failed to read next idx value");
        data.seek(idx_prev);
    };
};



} // end namespace BPCells