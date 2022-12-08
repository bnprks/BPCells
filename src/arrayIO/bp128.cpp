#include "bp128.h"

namespace BPCells {

BP128UIntReaderBase::BP128UIntReaderBase(uint32_t count) : count(count) {}

uint32_t BP128UIntReaderBase::size() const { return count; }

void BP128UIntReaderBase::seek(uint32_t new_pos) {
    if (pos % 128 != 0 && new_pos % 128 != 0 && new_pos / 128 == pos / 128) {
        pos = new_pos;
    } else {
        pos = new_pos;
        _seek();
    }
}

uint32_t BP128UIntReaderBase::load(uint32_t *out, uint32_t count) {
    uint32_t i = 0;
    if (pos % 128 != 0) {
        i = std::min(count, 128 - pos % 128);
        std::memmove(out, buf.get() + pos % 128, i * sizeof(uint32_t));
        pos += i;
    }
    for (; i + 128 <= count; i += 128) {
        load128(out + i);
        pos += 128;
    }
    if (i == 0) {
        // Handle leftovers. If we get here we know count < 128 and
        // we started aligned to a boundary of 128
        load128(buf.get());
        std::memmove(out, buf.get(), count * sizeof(uint32_t));
        pos += count;
        i += count;
    }
    return i;
}

uint32_t BP128UIntWriterBase::write(uint32_t *in, uint32_t count) {
    uint32_t i = 0;
    if (buf_pos != 0 || count < 128) {
        i = std::min(count, 128 - buf_pos);
        std::memmove(buf.get() + buf_pos, in, i * sizeof(uint32_t));
        buf_pos = (buf_pos + i) % 128;
        if (buf_pos == 0) pack128(buf.get());
    }
    for (; i + 128 <= count; i += 128) {
        pack128(in + i);
    }
    return i;
}

void BP128UIntWriterBase::finalize() {
    if (buf_pos != 0) {
        for (; buf_pos < 128; buf_pos++) {
            buf[buf_pos] = buf[buf_pos - 1];
        }
        pack128(buf.get());
    }
    finalizeWriters();
}

//################### BP128 (Vanilla) #########################################

BP128UIntReader::BP128UIntReader(UIntReader &&data, UIntReader &&idx, uint32_t count)
    : BP128UIntReaderBase(count)
    , data(std::move(data))
    , idx(std::move(idx))
    , prev_idx(this->idx.read_one()) {}

void BP128UIntReader::_seek() {
    idx.seek(pos / 128);
    prev_idx = idx.read_one();
    data.seek(prev_idx);
    if (pos % 128 != 0) load128(buf.get());
}

void BP128UIntReader::load128(uint32_t *out) {
    uint32_t next_idx = idx.read_one();
    uint32_t bits = (next_idx - prev_idx) / 4;
    data.ensureCapacity(bits * 4);
    simdunpack(data.data(), out, bits);
    data.advance(bits * 4);
    prev_idx = next_idx;
}

BP128UIntWriter::BP128UIntWriter(UIntWriter &&_data, UIntWriter &&_idx)
    : BP128UIntWriterBase()
    , data(std::move(_data))
    , idx(std::move(_idx)) {
    idx.write_one(0);
}

void BP128UIntWriter::pack128(uint32_t *in) {
    uint32_t bits = simdmaxbits(in);
    data.ensureCapacity(bits * 4);
    simdpack(in, data.data(), bits);
    data.advance(bits * 4);
    cur_idx += bits * 4;
    idx.write_one(cur_idx);
}

void BP128UIntWriter::finalizeWriters() {
    data.finalize();
    idx.finalize();
}

//######################## BP128 (D1) #########################################

BP128_D1_UIntReader::BP128_D1_UIntReader(
    UIntReader &&data, UIntReader &&idx, UIntReader &&starts, uint32_t count
)
    : BP128UIntReaderBase(count)
    , data(std::move(data))
    , idx(std::move(idx))
    , starts(std::move(starts))
    , prev_idx(this->idx.read_one()) {}

void BP128_D1_UIntReader::_seek() {
    idx.seek(pos / 128);
    starts.seek(pos / 128);
    prev_idx = idx.read_one();
    data.seek(prev_idx);
    if (pos % 128 != 0) load128(buf.get());
}

void BP128_D1_UIntReader::load128(uint32_t *out) {
    uint32_t next_idx = idx.read_one();
    uint32_t start = starts.read_one();
    uint32_t bits = (next_idx - prev_idx) / 4;
    data.ensureCapacity(bits * 4);
    simdunpackd1(start, data.data(), out, bits);
    data.advance(bits * 4);
    prev_idx = next_idx;
}

BP128_D1_UIntWriter::BP128_D1_UIntWriter(
    UIntWriter &&_data, UIntWriter &&_idx, UIntWriter &&_starts
)
    : BP128UIntWriterBase()
    , data(std::move(_data))
    , idx(std::move(_idx))
    , starts(std::move(_starts)) {
    idx.write_one(0);
}

void BP128_D1_UIntWriter::pack128(uint32_t *in) {
    uint32_t start = in[0];
    uint32_t bits = simdmaxbitsd1(start, in);
    data.ensureCapacity(bits * 4);
    simdpackd1(start, in, data.data(), bits);
    data.advance(bits * 4);
    cur_idx += bits * 4;
    idx.write_one(cur_idx);
    starts.write_one(start);
}

void BP128_D1_UIntWriter::finalizeWriters() {
    data.finalize();
    idx.finalize();
    starts.finalize();
}

//######################## BP128 (D1Z) #########################################

BP128_D1Z_UIntReader::BP128_D1Z_UIntReader(
    UIntReader &&data, UIntReader &&idx, UIntReader &&starts, uint32_t count
)
    : BP128UIntReaderBase(count)
    , data(std::move(data))
    , idx(std::move(idx))
    , starts(std::move(starts))
    , prev_idx(this->idx.read_one()) {}

void BP128_D1Z_UIntReader::_seek() {
    idx.seek(pos / 128);
    starts.seek(pos / 128);
    prev_idx = idx.read_one();
    data.seek(prev_idx);
    if (pos % 128 != 0) load128(buf.get());
}

void BP128_D1Z_UIntReader::load128(uint32_t *out) {
    uint32_t next_idx = idx.read_one();
    uint32_t start = starts.read_one();
    uint32_t bits = (next_idx - prev_idx) / 4;
    data.ensureCapacity(bits * 4);
    simdunpackd1z(start, data.data(), out, bits);
    data.advance(bits * 4);
    prev_idx = next_idx;
}

BP128_D1Z_UIntWriter::BP128_D1Z_UIntWriter(
    UIntWriter &&_data, UIntWriter &&_idx, UIntWriter &&_starts
)
    : BP128UIntWriterBase()
    , data(std::move(_data))
    , idx(std::move(_idx))
    , starts(std::move(_starts)) {
    idx.write_one(0);
}

void BP128_D1Z_UIntWriter::pack128(uint32_t *in) {
    uint32_t start = in[0];
    uint32_t bits = simdmaxbitsd1z(start, in);
    data.ensureCapacity(bits * 4);
    simdpackd1z(start, in, data.data(), bits);
    data.advance(bits * 4);
    cur_idx += bits * 4;
    idx.write_one(cur_idx);
    starts.write_one(start);
}

void BP128_D1Z_UIntWriter::finalizeWriters() {
    data.finalize();
    idx.finalize();
    starts.finalize();
}

//######################## BP128 (FOR) #########################################
// This is just FOR encoding using a constant 1 as the frame of reference.
// (i.e. to encode values of an integer sparse matrix)

BP128_FOR_UIntReader::BP128_FOR_UIntReader(UIntReader &&data, UIntReader &&idx, uint32_t count)
    : BP128UIntReaderBase(count)
    , data(std::move(data))
    , idx(std::move(idx))
    , prev_idx(this->idx.read_one()) {}

void BP128_FOR_UIntReader::_seek() {
    idx.seek(pos / 128);
    prev_idx = idx.read_one();
    data.seek(prev_idx);
    if (pos % 128 != 0) load128(buf.get());
}

void BP128_FOR_UIntReader::load128(uint32_t *out) {
    uint32_t next_idx = idx.read_one();
    uint32_t bits = (next_idx - prev_idx) / 4;
    data.ensureCapacity(bits * 4);
    simdunpackFOR(1, data.data(), out, bits);
    data.advance(bits * 4);
    prev_idx = next_idx;
}

BP128_FOR_UIntWriter::BP128_FOR_UIntWriter(UIntWriter &&_data, UIntWriter &&_idx)
    : BP128UIntWriterBase()
    , data(std::move(_data))
    , idx(std::move(_idx)) {
    idx.write_one(0);
}

void BP128_FOR_UIntWriter::pack128(uint32_t *in) {
    uint32_t bits = simdmaxbitsFOR(1, in);
    data.ensureCapacity(bits * 4);
    simdpackFOR(1, in, data.data(), bits);
    data.advance(bits * 4);
    cur_idx += bits * 4;
    idx.write_one(cur_idx);
}

void BP128_FOR_UIntWriter::finalizeWriters() {
    data.finalize();
    idx.finalize();
}

} // end namespace BPCells
