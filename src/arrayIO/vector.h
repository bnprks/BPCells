#pragma once

#include <stdint.h>
#include <vector>
#include <string>

#include "frags.h"
#include "array_interfaces.h"

namespace BPCells {
    

class ZVecUIntWriter : public UIntWriter {
private:
    std::vector<uint32_t> &vec;
    const uint32_t chunk_size = 8192;
    void _ensureCapacity(size_t capacity) override;
public:
    ~ZVecUIntWriter();
    ZVecUIntWriter() = delete;
    ZVecUIntWriter(std::vector<uint32_t> &vec);
    ZVecUIntWriter(std::vector<uint32_t> &vec, uint32_t chunk_size);
    
    void next() override;
};

class ZVecUIntReader : public UIntReader {
private:
    const uint32_t *vec;
    void _ensureCapacity(size_t capacity) override;
public:
    ZVecUIntReader(const uint32_t *vec, std::size_t capacity);

    // Return a buffer pointing to data available to read.
    // Pointers are only valid until `next` or `seek` gets called.
    // buffer size will always be 0 once there is no data left.
    bool next() override;

    // Seek to a different position in the stream (first integer is position 0)
    bool seek(size_t pos) override;
};

class VecStringWriter : public StringWriter {
private:
    std::vector<std::string> &data;
public:
    VecStringWriter(std::vector<std::string> &data);
    void write(const StringReader &reader) override;
};

class VecUIntWriter {
private:
    std::vector<uint32_t> &vec;
public:
    VecUIntWriter() = delete;
    VecUIntWriter(std::vector<uint32_t> &vec) : vec(vec) {}
    VecUIntWriter(const VecUIntWriter&) = default;
    VecUIntWriter& operator= (const VecUIntWriter&) = delete;

    void write(const uint32_t *buffer, uint32_t count);
    void finalize();
};




class VecUIntReader {
private:
    const uint32_t *vec;
    std::size_t capacity;
    std::size_t idx = 0;
public:
    VecUIntReader() = default;
    VecUIntReader(const uint32_t *vec, std::size_t capacity);

    uint32_t read(uint32_t *buffer, uint32_t count);
    uint32_t size();
    void seek(const std::size_t pos);
};


class VecUnpackedFragmentsSaver {
public:
    using UIntWriter = VecUIntWriter;
    class Storage {
    public:
        std::vector<UnpackedFrags<std::vector<uint32_t>>> fragments;
        std::vector<std::string> cell_names, chr_names;
    };

    VecUnpackedFragmentsSaver() = delete;
    VecUnpackedFragmentsSaver(Storage &storage) : storage(storage) {};
    UnpackedFrags<UIntWriter> chrWriterUnpacked(uint32_t chr_id);

    void writeCellNames(std::vector<std::string> cell_names);
    void writeChrNames(std::vector<std::string> chr_names);
private:
    Storage &storage;
};

class VecPackedFragmentsSaver {
public:
    using UIntWriter = VecUIntWriter;
    class Storage {
    public:
        std::vector<PackedFrags<std::vector<uint32_t>>> fragments;
        std::vector<std::string> cell_names, chr_names;
    };

    VecPackedFragmentsSaver() = delete;
    VecPackedFragmentsSaver(Storage &storage) : storage(storage) {};
    PackedFrags<UIntWriter> chrWriterPacked(uint32_t chr_id);
    void writeCellNames(std::vector<std::string> cell_names);
    void writeChrNames(std::vector<std::string> chr_names);
private:
    Storage &storage;
};

class VecUnpackedFragmentsLoader {
public:
    using UIntReader = VecUIntReader;
    class Storage {
    public:
        std::vector<UnpackedFrags<UIntReader>> fragments;
        std::vector<std::string> cell_names, chr_names;
    };

    VecUnpackedFragmentsLoader() = delete;
    VecUnpackedFragmentsLoader(Storage storage) : storage(storage) {};

    UnpackedFrags<UIntReader> chrReaderUnpacked(uint32_t chr_id) {
        return UnpackedFrags<UIntReader> {
            storage.fragments[chr_id].start,
            storage.fragments[chr_id].end,
            storage.fragments[chr_id].cell,
            storage.fragments[chr_id].end_max
        };
    }

    std::vector<std::string> readCellNames() {return storage.cell_names;}
    std::vector<std::string> readChrNames() {return storage.chr_names;}
private:
    const Storage storage;
};


class VecPackedFragmentsLoader {
public:
    using UIntReader = VecUIntReader;
    class Storage {
    public:
        std::vector<PackedFrags<UIntReader>> fragments;
        std::vector<std::string> cell_names, chr_names;
    };

    VecPackedFragmentsLoader() = delete;
    VecPackedFragmentsLoader(Storage storage) : storage(storage) {};
    PackedFrags<UIntReader> chrReaderPacked(uint32_t chr_id) {
        return PackedFrags<UIntReader> {
            storage.fragments[chr_id].start_data,
            storage.fragments[chr_id].start_idx,
            storage.fragments[chr_id].start_starts,
            storage.fragments[chr_id].end_data,
            storage.fragments[chr_id].end_idx,
            storage.fragments[chr_id].end_max,
            storage.fragments[chr_id].cell_data,
            storage.fragments[chr_id].cell_idx,
            storage.fragments[chr_id].count
        };
    }

    std::vector<std::string> readCellNames() {return storage.cell_names;}
    std::vector<std::string> readChrNames() {return storage.chr_names;}
private:
    const Storage storage;
};


} // end namespace BPCells