#include <filesystem>
#include <sstream>

#include <gtest/gtest.h>

#include <fragmentIterators/FragmentIterator.h>
#include <fragmentIterators/BedFragments.h>
#include <fragmentIterators/StoredFragments.h>
#include <fragmentIterators/ChrSelect.h>
#include <arrayIO/vector.h>

#include "utils-fragments.h"

namespace fs = std::filesystem;
using namespace BPCells;
using namespace ::testing;



TEST(FragmentIO, BedRoundtrip) {
    uint32_t max_cell = 50;
    auto frags_vec = Testing::generateFrags(2000, 3, 400, max_cell-1, 100, 1336);
    std::unique_ptr<VecReaderWriterBuilder> v = writeFragmentTuple(frags_vec);
    StoredFragments frags = StoredFragments::openUnpacked(*v);

    fs::path p = fs::temp_directory_path() / "BPCells_fragmentIO_test/fragments.tsv.gz";
    fs::create_directories(p.parent_path());
    if (fs::exists(p)) fs::remove(p);

    BedFragmentsWriter w(p.c_str());
    w.write(frags);

    BedFragments bed1(p.c_str());

    fs::path p2 = fs::temp_directory_path() / "BPCells_fragmentIO_test/fragments2.tsv.gz";
    if (fs::exists(p2)) fs::remove(p2);
    BedFragmentsWriter w2(p2.c_str());
    w2.write(bed1);

    std::stringstream command;
    command << "bash -c \"diff -q <(gunzip -c '" << p2.string() << "' | cut -f 1-4) "
        << "<(gunzip -c '" << p.string() << "')\"";
    EXPECT_EQ(0, std::system(command.str().c_str()));
}

TEST(FragmentIO, UnpackedVec) {
    // This test is a little redundant, since writeFragmentTuple already assumes the unpacked round-trip works
    uint32_t max_cell = 50;
    auto frags_vec = Testing::generateFrags(2000, 3, 400, max_cell-1, 100, 1336);
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
    auto frags_vec = Testing::generateFrags(2000, 3, 400, max_cell-1, 100, 1336);
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
    auto frags_vec = Testing::generateFrags(2000, 3, 400, max_cell-1, 100, 1336);
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
    // Failing thing
    BedFragments in("../../data/mini_fragments.tsv.gz");
    ChrNameSelect l1(in, {"chr1", "chr2", "chr3", "chrX", "chrY"});

    VecReaderWriterBuilder v;
    StoredFragmentsWriter::createPacked(v).write(in);
    StoredFragmentsPacked l2 = StoredFragmentsPacked::openPacked(v);

    VecReaderWriterBuilder v2;
    StoredFragmentsWriter::createPacked(v2).write(l2);
    StoredFragmentsPacked l3 = StoredFragmentsPacked::openPacked(v2);

    in.restart();
    ASSERT_TRUE(Testing::fragments_identical(in, l3));
}