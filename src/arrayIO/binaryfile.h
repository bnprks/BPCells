#pragma once

#include <cstring>
#include <fstream>
#include <vector>
#include <filesystem>
#include <arpa/inet.h>

#include "array_interfaces.h"
#include "frags.h"
#include "../fragmentIterators/UnpackedFragments2.h"
#include "../matrixIterators/PackedMatrix.h"
#include "../matrixIterators/UnpackedMatrix.h"

namespace BPCells {

std::string loadVersionMatrixDir(std::string dir);
UnpackedMatrix openUnpackedMatrixDir(std::string dir, uint32_t buffer_size);
UnpackedMatrixWriter createUnpackedMatrixDir(std::string dir, uint32_t buffer_size);
PackedMatrix openPackedMatrixDir(std::string dir, uint32_t buffer_size);
PackedMatrixWriter createPackedMatrixDir(std::string dir, uint32_t buffer_size);

class ZFileUIntWriter : public UIntWriter {
private:
    std::ofstream file;
    std::vector<uint32_t> buf;
    void flush();
    void _ensureCapacity(size_t capacity) override;
public:
    ~ZFileUIntWriter();
    ZFileUIntWriter(const char* path, uint32_t buffer_size = 8192);
    ZFileUIntWriter(const ZFileUIntWriter&) = delete;
    ZFileUIntWriter& operator= (const ZFileUIntWriter&) = delete;
    ZFileUIntWriter& operator= (ZFileUIntWriter&&) = default;

    void next() override;
};


class ZFileUIntReader : public UIntReader {
protected:
    std::ifstream file;
    std::vector<uint32_t> buf;
    void _ensureCapacity(size_t capacity) override;
    size_t readPos(uint32_t* out, uint32_t capacity);
public:
    ZFileUIntReader() = default;
    ZFileUIntReader(const char* path, uint32_t buffer_size = 8192);
    ZFileUIntReader(const ZFileUIntReader&) = delete;
    ZFileUIntReader& operator= (const ZFileUIntReader&) = delete;
    ZFileUIntReader& operator= (ZFileUIntReader&&) = default;

    bool next() override;
    bool seek(size_t pos) override;
};

std::vector<std::string> readLines(std::filesystem::path path);

class FileUIntWriter {
private:
    std::ofstream file;
    uint32_t order_buf[128];
public:
    FileUIntWriter(const char* path);
    FileUIntWriter(const FileUIntWriter&) = delete;
    FileUIntWriter& operator= (const FileUIntWriter&) = delete;
    FileUIntWriter& operator= (FileUIntWriter&&) = default;

    void write(const uint32_t *buffer, uint32_t count);
    void finalize();
};

class FileUIntReader {
private:
    std::ifstream file;
public:
    FileUIntReader() = default;
    FileUIntReader(const char* path);
    FileUIntReader(const FileUIntReader&) = delete;
    FileUIntReader& operator= (const FileUIntReader&) = delete;
    FileUIntReader& operator= (FileUIntReader&&) = default;

    uint32_t read(uint32_t *buffer, uint32_t count);
    uint32_t size();
    void seek(const size_t pos);
};

// Handles both packed and unpacked fragment saving. For consistency,
// users should call only one of chrWriterUnpacked or chrWriterPacked according
// to the desired ouptut format
class FileFragmentsSaver {
private:
    std::filesystem::path dir;
public:
    using UIntWriter = FileUIntWriter;
    FileFragmentsSaver() = delete;
    FileFragmentsSaver(std::string dir);
    UnpackedFrags<UIntWriter> chrWriterUnpacked(uint32_t chr_id);
    PackedFrags<UIntWriter> chrWriterPacked(uint32_t chr_id);
    void writeCellNames(std::vector<std::string> cell_names);
    void writeChrNames(std::vector<std::string> chr_names);
};

// Handles both packed and unpacked fragment loading. For consistency,
// users should call only one of chrReaderUnpacked or chrReaderPacked according
// to the input format
class FileFragmentsLoader {
private:
    std::filesystem::path dir;
public:
    using UIntReader = FileUIntReader;
    FileFragmentsLoader() = default;
    FileFragmentsLoader(std::string dir, bool is_packed);

    UnpackedFrags<UIntReader> chrReaderUnpacked(uint32_t chr_id);
    PackedFrags<UIntReader> chrReaderPacked(uint32_t chr_id);

    std::vector<std::string> readCellNames() {return readLines(dir / "cell_names.txt");}
    std::vector<std::string> readChrNames() {return readLines(dir / "chr_names.txt");}
};

// UnpackedFragments2(std::vector<UnpackedFrags<Reader> > frags,
//         const std::vector<std::string> cell_names, const std::vector<std::string> chr_names)

class FilePackedFragsBuilder {
private:
    std::filesystem::path dir;
public:
    FilePackedFragsBuilder() = delete;
    FilePackedFragsBuilder(std::string dir);
    PackedFrags<FileUIntWriter> chrWriter(uint32_t chr_id);
    void writeCellNames(std::vector<std::string> cell_names);
    void writeChrNames(std::vector<std::string> chr_names);
};




} // end namespace BPCells