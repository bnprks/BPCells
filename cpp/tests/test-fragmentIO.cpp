// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include <utils/filesystem_compat.h>
#include <sstream>

#include <gtest/gtest.h>

#include <arrayIO/vector.h>
#include <fragmentIterators/BedFragments.h>
#include <fragmentIterators/ChrSelect.h>
#include <fragmentIterators/FragmentIterator.h>
#include <fragmentIterators/StoredFragments.h>

#include "utils-fragments.h"

using namespace BPCells;
using namespace ::testing;

char **my_argv;
int my_argc;

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    my_argc = argc;
    my_argv = argv;

    return RUN_ALL_TESTS();
}

TEST(FragmentIO, BedRoundtrip) {
    uint32_t max_cell = 50;
    auto frags_vec = Testing::generateFrags(2000, 3, 400, max_cell - 1, 100, 1336);
    std::unique_ptr<VecReaderWriterBuilder> v = writeFragmentTuple(frags_vec);
    StoredFragments frags = StoredFragments::openUnpacked(*v);

    std_fs::path p = std_fs::temp_directory_path() / "BPCells_fragmentIO_test/fragments.tsv.gz";
    std_fs::create_directories(p.parent_path());
    if (std_fs::exists(p)) std_fs::remove(p);

    BedFragmentsWriter w(p.string().c_str());
    w.write(frags);

    BedFragments bed1(p.string().c_str());

    std_fs::path p2 = std_fs::temp_directory_path() / "BPCells_fragmentIO_test/fragments2.tsv.gz";
    if (std_fs::exists(p2)) std_fs::remove(p2);
    BedFragmentsWriter w2(p2.string().c_str());
    w2.write(bed1);

    std::stringstream command;
    command << "bash -c \"diff -q <(gunzip -c '" << p2.string() << "' | cut -f 1-4) "
            << "<(gunzip -c '" << p.string() << "')\"";
    EXPECT_EQ(0, std::system(command.str().c_str()));
}

TEST(FragmentIO, UnpackedVec) {
    // This test is a little redundant, since writeFragmentTuple already assumes
    // the unpacked round-trip works
    uint32_t max_cell = 50;
    auto frags_vec = Testing::generateFrags(2000, 3, 400, max_cell - 1, 100, 1336);
    std::unique_ptr<VecReaderWriterBuilder> v = writeFragmentTuple(frags_vec);
    StoredFragments in = StoredFragments::openUnpacked(*v);

    // Write properly to memory
    VecReaderWriterBuilder vb1(1024);
    auto w1 = StoredFragmentsWriter::createUnpacked(vb1);
    w1.write(in);

    auto loader = StoredFragments::openUnpacked(vb1);
    in.restart();
    ASSERT_TRUE(Testing::fragments_identical(loader, in));
}

TEST(FragmentIO, PackedVec) {
    uint32_t max_cell = 50;
    auto frags_vec = Testing::generateFrags(2000, 3, 400, max_cell - 1, 100, 1336);
    std::unique_ptr<VecReaderWriterBuilder> v = writeFragmentTuple(frags_vec);
    StoredFragments in = StoredFragments::openUnpacked(*v);

    // Read from disk to memory
    VecReaderWriterBuilder vb1(1024);
    auto w1 = StoredFragmentsWriter::createPacked(vb1);
    w1.write(in);

    SCOPED_TRACE("PackedVec");
    auto loader = StoredFragmentsPacked::openPacked(vb1);
    in.restart();
    ASSERT_TRUE(Testing::fragments_identical(loader, in));
}

TEST(FragmentIO, ReducedCapacityWrite) {
    // Test writing StoredFragments when the reader loads more at a time than the
    // chunk size of the writer
    uint32_t max_cell = 50;
    auto frags_vec = Testing::generateFrags(2000, 3, 400, max_cell - 1, 100, 1336);
    // Write with chunk size of 1024
    std::unique_ptr<VecReaderWriterBuilder> v = writeFragmentTuple(frags_vec, 0, false, 1024);
    StoredFragments in = StoredFragments::openUnpacked(*v);

    // Write properly to memory, but with chunk size 50
    VecReaderWriterBuilder vb1(50);
    auto w1 = StoredFragmentsWriter::createUnpacked(vb1);
    w1.write(in);

    auto loader = StoredFragments::openUnpacked(vb1);
    in.restart();
    ASSERT_TRUE(Testing::fragments_identical(loader, in));

    // Can't write packed with size < 128, but we can try 129
    VecReaderWriterBuilder vb2(129);
    auto w2 = StoredFragmentsWriter::createPacked(vb2);
    in.restart();
    w2.write(in);

    auto loader2 = StoredFragmentsPacked::openPacked(vb2);
    in.restart();
    ASSERT_TRUE(Testing::fragments_identical(loader2, in));
}

TEST(FragmentIO, ChrSelectRoundTrip) {
    // Slight hack to get path of an input test data file via CMake config
    ASSERT_EQ(my_argc, 2);
    BedFragments in(my_argv[1]);
    ChrNameSelect l1(std::make_unique<BedFragments>(my_argv[1]), {"chr1", "chr2", "chr3", "chrX", "chrY"});

    VecReaderWriterBuilder v;
    StoredFragmentsWriter::createPacked(v).write(in);
    StoredFragmentsPacked l2 = StoredFragmentsPacked::openPacked(v);

    VecReaderWriterBuilder v2;
    StoredFragmentsWriter::createPacked(v2).write(l2);
    StoredFragmentsPacked l3 = StoredFragmentsPacked::openPacked(v2);

    in.restart();
    ASSERT_TRUE(Testing::fragments_identical(in, l3));
}