#include <filesystem>
#include <sstream>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <fragmentIterators/FragmentIterator.h>
#include <fragmentIterators/BedFragments.h>
#include <fragmentIterators/StoredFragments.h>
#include <arrayIO/vector.h>

namespace fs = std::filesystem;
using namespace BPCells;
using namespace ::testing;


const char* data_path;

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    
    if (argc <= 1) {
        data_path = "/Users/ben/Dropbox/greenleaf/playground/fragment_io/BPCells/tests/data/mini_fragments.tsv.gz";
    } else {
        data_path = argv[1];
    }

    return RUN_ALL_TESTS();
}

//
//TEST(FragmentIO, BedRoundtrip) {
//    BedFragments bed(data_path);
//    FragmentIterator in(bed);
//
//    fs::path p = fs::temp_directory_path() / "BPCells_fragmentIO_test/fragments.tsv.gz";
//    fs::create_directories(p.parent_path());
//    if (fs::exists(p)) fs::remove(p);
//
//    BedFragmentsWriter w(p.c_str());
//    w.write(in);
//
//    std::stringstream command;
//    command << "bash -c \"diff -q <(gunzip -c '" << std::string_view(data_path) << "' | cut -f 1-4) "
//        << "<(gunzip -c '" << p.string() << "')\"";
//    EXPECT_EQ(0, std::system(command.str().c_str()));
//}

void equal_vec(std::vector<uint32_t> v1, std::vector<uint32_t> v2) {
    ASSERT_EQ(v1.size(), v2.size());
    for (int i = 0; i < v1.size(); i++) {
        ASSERT_EQ(v1[i], v2[i]);
    }
}

void assert_equal_fragments(FragmentLoader &l1, FragmentLoader &l2) {
    l1.restart(); l2.restart();
    FragmentIterator i1(l1);
    FragmentIterator i2(l2);

    while(true) {
        bool res1 = i1.nextChr();
        bool res2 = i2.nextChr();
        ASSERT_EQ(res1, res2);
        if (!res1) break;
        ASSERT_EQ(i1.currentChr(), i2.currentChr());
        while(true) {
            bool res1 = i1.nextFrag();
            bool res2 = i2.nextFrag();
            ASSERT_EQ(res1, res2);
            if (!res1) break;
            ASSERT_EQ(i1.cell(), i2.cell());
            ASSERT_EQ(i1.start(), i2.start());
            ASSERT_EQ(i1.end(), i2.end());
        }
    }
    for (uint32_t i = 0; ;i++) {
        ASSERT_STREQ(i1.cellNames(i), i2.cellNames(i));
        if (i1.cellNames(i) == NULL) break;
    }
    for (uint32_t i = 0; ;i++) {
        ASSERT_STREQ(i1.chrNames(i), i2.chrNames(i));
        if (i1.chrNames(i) == NULL) break;
    }
}

TEST(FragmentIO, UnpackedVec) {
    BedFragments bed(data_path);
    FragmentIterator in(bed);

    // Read from disk to memory
    VecReaderWriterBuilder vb1(1024);
    auto w1 = StoredFragmentsWriter::createUnpacked(vb1);
    w1.write(in);

    SCOPED_TRACE("UnpackedVec");
    auto loader = StoredFragments::openUnpacked(vb1);
    assert_equal_fragments(
        loader,
        bed
    );
}

TEST(FragmentIO, PackedVec) {
    BedFragments bed(data_path);
    FragmentIterator in(bed);

    // Read from disk to memory
    VecReaderWriterBuilder vb1(1024);
    auto w1 = StoredFragmentsWriter::createPacked(vb1);
    w1.write(in);

    SCOPED_TRACE("PackedVec");
    auto loader = StoredFragmentsPacked::openPacked(vb1);
    assert_equal_fragments(
        loader,
        bed
    );
}