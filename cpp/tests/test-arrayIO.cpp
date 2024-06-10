// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include <arrayIO/array_interfaces.h>
#include <arrayIO/binaryfile.h>
#include <arrayIO/bp128.h>
#include <arrayIO/hdf5.h>
#include <arrayIO/vector.h>
#include <gtest/gtest.h>

#include <utils/filesystem_compat.h>

using namespace BPCells;

// Basic Idea for the test:
// Write numbers 1 - 10K in an array
// Read once through in order, checking values
// Read once through while seeking forward by 1, doubling the seek distance

// This isn't particularly stressful on the bitpacking readers/writers, but
// it should provide a basic sanity test

// Write integers 1 to vals in a variety of write patterns
void writeValues(UIntWriter &w, uint32_t vals = 10000) {
    uint32_t i = 1;

    // Writing data in full blocks
    while (i < vals / 2) {
        for (int j = 0; j < w.capacity(); j++) {
            w.data()[j] = i++;
        }
        w.advance(w.capacity());
        w.ensureCapacity();
    }

    // Writing data in specific capacity amounts
    while (i + 200 < vals) {
        w.ensureCapacity(200);
        for (int j = 0; j < 200; j++) {
            w.data()[j] = i++;
        }
        w.advance(200);
    }

    // Writing data one at a time
    while (i <= vals) {
        w.write_one(i++);
    }
    w.finalize();
}

// Read integers from a pre-written array containing values 1-vals with
// a variety of write and seek patterns
void readValues(UIntReader &r, uint32_t vals = 10000) {
    uint32_t i = 0;
    while (i < 200) {
        ASSERT_EQ(i + 1, r.read_one());
        i++;
    }

    // Read data in 200-length blocks
    while (i < vals / 2) {
        r.ensureCapacity(200);
        for (int j = 0; j < 200; j++) {
            ASSERT_EQ(r.data()[j], i + 1);
            i++;
        }
        r.advance(200);
    }

    // Read data in full blocks
    while (i < vals * 2) {
        if (r.capacity() == 0 && !r.requestCapacity()) break;
        ASSERT_EQ(*r.data(), i + 1);
        i += 1;
        r.advance(1);
    }

    ASSERT_EQ(i, vals);

    // Seek data forwards
    for (int j = 0; j < vals * 4; j = j * 2 + 1) {
        r.seek(j);
        bool res = r.requestCapacity();
        ASSERT_EQ(res, j < vals);
        if (!res) break;
        ASSERT_EQ(r.data()[0], j + 1);
    }
    r.seek(vals - 1);
    ASSERT_TRUE(r.requestCapacity());
    ASSERT_EQ(r.data()[0], vals);
    ASSERT_EQ(r.capacity(), 1);

    // Seek data jumping forwards and backwards
    for (uint64_t offset = 0; offset < vals/2; offset += 127) {
        uint32_t val;
        for (auto pos : {offset, vals - offset - 1, offset + 5, vals - offset - 1 - 5}) {
            r.seek(pos);
            bool res = r.requestCapacity();
            ASSERT_EQ(res, pos < vals);
            if (!res) break;
            ASSERT_EQ(r.data()[0], pos + 1);
        }
    }

    // Seek beyond end of data
    r.seek(vals);
    ASSERT_FALSE(r.requestCapacity());
}

TEST(ArrayIO, Vector) {
    SCOPED_TRACE("Vector ArrayIO");
    std::vector<uint32_t> v(0);
    {
        UIntWriter w(std::make_unique<VecUIntWriter>(v), 1024);
        writeValues(w);
    }
    ASSERT_EQ(v.size(), 10000);
    for (int i = 0; i < v.size(); i++) {
        ASSERT_EQ(v[i], i + 1);
    }
    UIntReader r(std::make_unique<VecUIntReader>(v.data(), v.size()), 2040, 1024);
    readValues(r);
}

TEST(ArrayIO, Binaryfile) {
    SCOPED_TRACE("Binaryfile ArrayIO");
    std_fs::path p = std_fs::temp_directory_path() / "BPCells_arrayIO_test/bin_array";
    std_fs::create_directories(p.parent_path());
    if (std_fs::exists(p)) std_fs::remove(p);
    {
        UIntWriter w(std::make_unique<FileUIntWriter>(p.string().c_str()), 1024);
        writeValues(w);
    }
    UIntReader r(std::make_unique<FileUIntReader>(p.string().c_str()), 2040, 1024);
    readValues(r);
}

TEST(ArrayIO, HDF5) {
    SCOPED_TRACE("HDF5 ArrayIO");
    std_fs::path p = std_fs::temp_directory_path() / "BPCells_arrayIO_test/bin_array.hdf5";
    std_fs::create_directories(p.parent_path());
    if (std_fs::exists(p)) std_fs::remove(p);
    H5WriterBuilder write_builder(p.string(), "", 2040, 1024);
    H5ReaderBuilder read_builder(p.string(), "", 2040, 1024);
    {
        UIntWriter w = write_builder.createUIntWriter("my_group");
        writeValues(w);
    }
    UIntReader r = read_builder.openUIntReader("my_group");
    readValues(r);
}

TEST(ArrayIO, BP128) {
    SCOPED_TRACE("BP128 ArrayIO");
    std::vector<uint32_t> data(0);
    std::vector<uint32_t> idx(0);
    std::vector<uint64_t> idx_offsets(0);
    uint32_t vals = 10000;
    {
        UIntWriter w(
            std::make_unique<BP128UIntWriter>(
                UIntWriter(std::make_unique<VecUIntWriter>(data), 1018),
                UIntWriter(std::make_unique<VecUIntWriter>(idx), 1020),
                ULongWriter(std::make_unique<VecNumWriter<uint64_t>>(idx_offsets), 1018)
            ),
            1019
        );
        writeValues(w, vals);
    }
    ASSERT_GT(data.size(), 0);
    ASSERT_GT(idx.size(), 0);
    UIntReader r(
        std::make_unique<BP128UIntReader>(
            UIntReader(std::make_unique<VecUIntReader>(data.data(), data.size()), 2039, 1023),
            UIntReader(std::make_unique<VecUIntReader>(idx.data(), idx.size()), 2040, 1024),
            ULongReader(
                std::make_unique<VecNumReader<uint64_t>>(idx_offsets.data(), idx_offsets.size()),
                1018,
                14
            ),
            vals
        ),
        2041,
        1025
    );
    readValues(r, vals);
}

TEST(ArrayIO, BP128_D1) {
    SCOPED_TRACE("BP128_D1 ArrayIO");
    std::vector<uint32_t> data(0);
    std::vector<uint32_t> idx(0);
    std::vector<uint32_t> starts(0);
    std::vector<uint64_t> idx_offsets(0);
    uint32_t vals = 10000;
    {
        UIntWriter w(
            std::make_unique<BP128_D1_UIntWriter>(
                UIntWriter(std::make_unique<VecUIntWriter>(data), 1020),
                UIntWriter(std::make_unique<VecUIntWriter>(idx), 1019),
                ULongWriter(std::make_unique<VecNumWriter<uint64_t>>(idx_offsets), 1019),
                UIntWriter(std::make_unique<VecUIntWriter>(starts), 1021)
            ),
            1018
        );
        writeValues(w, vals);
    }

    UIntReader r(
        std::make_unique<BP128_D1_UIntReader>(
            UIntReader(std::make_unique<VecUIntReader>(data.data(), data.size()), 2039, 1023),
            UIntReader(std::make_unique<VecUIntReader>(idx.data(), idx.size()), 2040, 1024),
            ULongReader(
                std::make_unique<VecNumReader<uint64_t>>(idx_offsets.data(), idx_offsets.size()),
                1019,
                14
            ),
            UIntReader(std::make_unique<VecUIntReader>(starts.data(), starts.size()), 2039, 1023),
            vals
        ),
        2041,
        1025
    );
    readValues(r, vals);
}

TEST(ArrayIO, BP128_D1Z) {
    SCOPED_TRACE("BP128_D1Z ArrayIO");
    std::vector<uint32_t> data(0);
    std::vector<uint32_t> idx(0);
    std::vector<uint64_t> idx_offsets(0);
    std::vector<uint32_t> starts(0);
    uint32_t vals = 10000;
    {
        UIntWriter w(
            std::make_unique<BP128_D1Z_UIntWriter>(
                UIntWriter(std::make_unique<VecUIntWriter>(data), 1020),
                UIntWriter(std::make_unique<VecUIntWriter>(idx), 1019),
                ULongWriter(std::make_unique<VecNumWriter<uint64_t>>(idx_offsets), 1019),
                UIntWriter(std::make_unique<VecUIntWriter>(starts), 1021)
            ),
            1018
        );
        writeValues(w, vals);
    }

    UIntReader r(
        std::make_unique<BP128_D1Z_UIntReader>(
            UIntReader(std::make_unique<VecUIntReader>(data.data(), data.size()), 2039, 1023),
            UIntReader(std::make_unique<VecUIntReader>(idx.data(), idx.size()), 2040, 1024),
            ULongReader(
                std::make_unique<VecNumReader<uint64_t>>(idx_offsets.data(), idx_offsets.size()),
                1019,
                14
            ),
            UIntReader(std::make_unique<VecUIntReader>(starts.data(), starts.size()), 2039, 1023),
            vals
        ),
        2041,
        1025
    );
    readValues(r, vals);
}

TEST(ArrayIO, BP128_FOR) {
    SCOPED_TRACE("BP128_FOR ArrayIO");
    std::vector<uint32_t> data(0);
    std::vector<uint32_t> idx(0);
    std::vector<uint64_t> idx_offsets(0);
    uint32_t vals = 10000;
    {
        UIntWriter w(
            std::make_unique<BP128_FOR_UIntWriter>(
                UIntWriter(std::make_unique<VecUIntWriter>(data), 1020),
                UIntWriter(std::make_unique<VecUIntWriter>(idx), 1019),
                ULongWriter(std::make_unique<VecNumWriter<uint64_t>>(idx_offsets), 1019)
            ),
            1018
        );
        writeValues(w, vals);
    }

    UIntReader r(
        std::make_unique<BP128_FOR_UIntReader>(
            UIntReader(std::make_unique<VecUIntReader>(data.data(), data.size()), 2039, 1023),
            UIntReader(std::make_unique<VecUIntReader>(idx.data(), idx.size()), 2040, 1024),
            ULongReader(
                std::make_unique<VecNumReader<uint64_t>>(idx_offsets.data(), idx_offsets.size()),
                1019,
                14
            ),
            vals
        ),
        2041,
        1025
    );
    readValues(r, vals);
}

TEST(ArrayIO, BP128_idx_offset_incompressible) {
    // Steps: 
    // 1. Set the idx_offset wraparound point abnormally low (4096)
    // 2. Make a compressed dataset, where the values at index i are UINT32_MAX - i % 65537
    // 3. Read forward through whole data
    // 4. Seek backwards and forwards and read data
    std::vector<uint32_t> data, idx;
    std::vector<uint64_t> idx_offsets;

    // 1. Set the idx_offset wraparound point abnormally low (4096 - 1)
    uint64_t OFFSET = 4096-1;
    BP128UIntReader::setOffsetIncrement(OFFSET);
    BP128UIntWriter::setOffsetIncrement(OFFSET);

    // 2. Make a compressed dataset, where the values at index i are UINT32_MAX - i % 65537
    BP128UIntWriter w(
        UIntWriter(std::make_unique<VecUIntWriter>(data), 1023),
        UIntWriter(std::make_unique<VecUIntWriter>(idx), 1022),
        ULongWriter(std::make_unique<VecNumWriter<uint64_t>>(idx_offsets), 1025)
    );
    uint64_t count = 100000;
    uint64_t chunk_size = 127;
    uint32_t buf[1024];
    for (uint64_t i = 0; i < count; ) {
        for (int j = 0; j < chunk_size; j++) {
            buf[j] = UINT32_MAX - (i + j) % 65537;
        }
        i += w.write(buf, std::min(chunk_size, count - i));
    }
    w.finalize();

    for (auto x : idx) ASSERT_LT(x, OFFSET);
    for (uint64_t i = 0; i < idx.size(); i++) {
        ASSERT_EQ(idx[i], (i*128) % OFFSET);
    }

    BP128UIntReader r(
        UIntReader(std::make_unique<VecUIntReader>(data.data(), data.size()), 2039, 1023),
        UIntReader(std::make_unique<VecUIntReader>(idx.data(), idx.size()), 2040, 1024),
        ULongReader(
            std::make_unique<VecNumReader<uint64_t>>(idx_offsets.data(), idx_offsets.size()),
            1019,
            14
        ),
        count
    );

    // 3. Read forward through whole data
    uint64_t read = 0;
    while (read < count) {
        uint64_t loaded = r.load(buf, 1023);
        for (uint64_t i = 0; i < loaded; i++) {
            ASSERT_EQ(buf[i], UINT32_MAX - (read + i) % 65537);
        }
        read += count;
    }
    ASSERT_EQ(read, count);

    // 4. Seek backwards and forwards and read data
    // Seek symmetric offsets from
    for (uint64_t offset = 0; offset < count/2; offset += 127) {
        uint32_t val;
        for (auto pos : {offset, count - offset - 1, offset + 5, count - offset - 1 - 5}) {
            r.seek(pos);
            r.load(&val, 1);
            ASSERT_EQ(val, UINT32_MAX - pos % 65537);
        }
    }

    BP128UIntReader::setOffsetIncrement(UINT32_MAX);
    BP128UIntWriter::setOffsetIncrement(UINT32_MAX);
}


TEST(ArrayIO, BP128_idx_offset_compressible) {
    SCOPED_TRACE("BP128 ArrayIO");

    uint64_t OFFSET = 4096-1;
    BP128UIntReader::setOffsetIncrement(OFFSET);
    BP128UIntWriter::setOffsetIncrement(OFFSET);

    std::vector<uint32_t> data(0);
    std::vector<uint32_t> idx(0);
    std::vector<uint64_t> idx_offsets(0);
    uint32_t vals = 100000;
    {
        UIntWriter w(
            std::make_unique<BP128UIntWriter>(
                UIntWriter(std::make_unique<VecUIntWriter>(data), 1018),
                UIntWriter(std::make_unique<VecUIntWriter>(idx), 1020),
                ULongWriter(std::make_unique<VecNumWriter<uint64_t>>(idx_offsets), 1018)
            ),
            1019
        );
        writeValues(w, vals);
    }
    ASSERT_GT(data.size(), 0);
    ASSERT_GT(idx.size(), 0);
    ASSERT_GT(idx_offsets.size(), 2);
    UIntReader r(
        std::make_unique<BP128UIntReader>(
            UIntReader(std::make_unique<VecUIntReader>(data.data(), data.size()), 2039, 1023),
            UIntReader(std::make_unique<VecUIntReader>(idx.data(), idx.size()), 2040, 1024),
            ULongReader(
                std::make_unique<VecNumReader<uint64_t>>(idx_offsets.data(), idx_offsets.size()),
                1018,
                14
            ),
            vals
        ),
        2041,
        1025
    );
    readValues(r, vals);
    BP128UIntReader::setOffsetIncrement(UINT32_MAX);
    BP128UIntWriter::setOffsetIncrement(UINT32_MAX);
}