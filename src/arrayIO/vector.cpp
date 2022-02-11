#include "vector.h"
#include <cstring>
namespace BPCells {


ZVecUIntWriter::ZVecUIntWriter(std::vector<uint32_t> &vec): vec(vec) {}
ZVecUIntWriter::ZVecUIntWriter(std::vector<uint32_t> &vec, uint32_t chunk_size) : vec(vec), chunk_size(chunk_size) {}

ZVecUIntWriter::~ZVecUIntWriter() {
    vec.resize(data() - vec.data());
}

void ZVecUIntWriter::_ensureCapacity(size_t new_capacity) {
    new_capacity = std::min((size_t) chunk_size, new_capacity);
    uint32_t num_written = write_buffer.data - vec.data();
    vec.resize(new_capacity + num_written);
    write_buffer.data = vec.data() + num_written;
    write_buffer.size = vec.size() - num_written;
}

void ZVecUIntWriter::next() {
    size_t initial_size = data() + capacity() - vec.data();
    vec.resize(initial_size + chunk_size);
    write_buffer.data = vec.data() + initial_size;
    write_buffer.size = chunk_size;
}

ZVecUIntReader::ZVecUIntReader(const uint32_t *vec, std::size_t capacity) :
    vec(vec) {
    total_size = capacity;
}

void ZVecUIntReader::_ensureCapacity(size_t new_capacity) {
    // An invariant of this reader is data()+capcity() == vec + total_size,
    // So if _ensureCapacity is called, that means there is no space for new_capacity
    throw std::runtime_error("Requested read capacity larger than buffer size");
};

bool ZVecUIntReader::next() {
    if (data() == NULL) {
        read_buffer.data = vec;
        read_buffer.size = total_size;
        return true;
    } else {
        return false;
    }
}

// Seek to a different position in the stream (first integer is position 0)
bool ZVecUIntReader::seek(size_t pos) {
    if (pos > capacity()) return false;
    // invariant: data() + capacity() == vec + total_size
    read_buffer.data = vec + pos;
    read_buffer.size = total_size - pos;
    return true;
};

VecStringWriter::VecStringWriter(std::vector<std::string> &data) : data(data) {}
void VecStringWriter::write(const StringReader &reader) {
    uint32_t i = 0;
    data.resize(0);
    while (true) {
        const char* s = reader.get(i);
        if (s == NULL) break;
        data.push_back(s);
        i++;
    }
}


void VecUIntWriter::write(const uint32_t *buffer, uint32_t count) {
    for (std::size_t i = 0; i < count; i++) {
        vec.push_back(buffer[i]);
    }
}

void VecUIntWriter::finalize() {return;} // no-op since vectors are always written

VecUIntReader::VecUIntReader(const uint32_t *vec, std::size_t capacity) : vec(vec), capacity(capacity) {}

uint32_t VecUIntReader::read(uint32_t *buffer, uint32_t count) {
    uint32_t n = std::min(count, (uint32_t) (capacity - idx));
    std::memmove(buffer, &vec[idx], n * sizeof(uint32_t));
    idx += n;
    return n;
}
uint32_t VecUIntReader::size() {
    return capacity;
}
void VecUIntReader::seek(const std::size_t pos) {
    idx = pos;
}



UnpackedFrags<VecUIntWriter> VecUnpackedFragmentsSaver::chrWriterUnpacked(uint32_t chr_id) {
    if (storage.fragments.size() <= chr_id) {
        storage.fragments.resize(chr_id + 1);
    }
    return UnpackedFrags<VecUIntWriter> {
        VecUIntWriter(storage.fragments[chr_id].start),
        VecUIntWriter(storage.fragments[chr_id].end),
        VecUIntWriter(storage.fragments[chr_id].cell),
        VecUIntWriter(storage.fragments[chr_id].end_max)
    };
}
void VecUnpackedFragmentsSaver::writeCellNames(std::vector<std::string> cell_names) {
    storage.cell_names = cell_names;
}
void VecUnpackedFragmentsSaver::writeChrNames(std::vector<std::string> chr_names) {
    storage.chr_names = chr_names;
}

PackedFrags<VecUIntWriter> VecPackedFragmentsSaver::chrWriterPacked(uint32_t chr_id) {
    if (storage.fragments.size() <= chr_id) {
        storage.fragments.resize(chr_id + 1);
    }
    //T start_data, start_idx, start_starts,
    //    end_data, end_idx, end_max,
    //   cell_data, cell_idx, count;
    return PackedFrags<VecUIntWriter> {
        VecUIntWriter(storage.fragments[chr_id].start_data),
        VecUIntWriter(storage.fragments[chr_id].start_idx),
        VecUIntWriter(storage.fragments[chr_id].start_starts),
        VecUIntWriter(storage.fragments[chr_id].end_data),
        VecUIntWriter(storage.fragments[chr_id].end_idx),
        VecUIntWriter(storage.fragments[chr_id].end_max),
        VecUIntWriter(storage.fragments[chr_id].cell_data),
        VecUIntWriter(storage.fragments[chr_id].cell_idx),
        VecUIntWriter(storage.fragments[chr_id].count)
    };
}
void VecPackedFragmentsSaver::writeCellNames(std::vector<std::string> cell_names) {
    storage.cell_names = cell_names;
}
void VecPackedFragmentsSaver::writeChrNames(std::vector<std::string> chr_names) {
    storage.chr_names = chr_names;
}



} // end namespace BPCells