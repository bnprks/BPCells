#include <gtest/gtest.h>

#include "utils-fragments.h"

#include <fragmentIterators/FragmentIterator.h>
#include <fragmentIterators/StoredFragments.h>
#include <fragmentIterators/RegionSelect.h>
#include <arrayIO/vector.h>

using namespace BPCells;



TEST(FragmentSeeking, SeekStoredFrags) {
    using namespace Testing;
    std::vector<Frag> frags_vec;
    // Write chr1 test data, zig-zag up+down
    
    for (uint32_t start = 0; start < 2000; start += 500) {
        for (uint32_t end = start+1; end < start+250; end++) {
            frags_vec.push_back({0, start, end, 0});
        }
        for (uint32_t end = start+250; end > start; end--) {
            frags_vec.push_back({0, start, end, 0});
        }
    }
    // Write chr2 test data, 128 frags not aligned on a 128-block boundary in stored array
    // Smaller ends than previous chr
    for (uint32_t i = 0; i < 128; i++) {
        frags_vec.push_back({1, i+10, i+1000, 0});
    }

    // Write chr2 test data, 128 frags not aligned on a 128-block boundary in stored array
    // Larger ends than previous chr
    for (uint32_t i = 0; i < 128; i++) {
        frags_vec.push_back({2, i+10, i+2000, 0});
    }

    std::unique_ptr<VecReaderWriterBuilder> v = writeFragmentTuple(frags_vec);
    StoredFragments frags = StoredFragments::openUnpacked(*v);

    // Try seeking a couple locations in chr1 to make sure basic functionality works
    frags.seek(0, 1249);
    EXPECT_TRUE(frags.load());
    EXPECT_EQ(frags.startData()[0], 1000);
    EXPECT_GT(frags.endData()[0], 1249-128);
    frags.seek(0, 250);
    EXPECT_TRUE(frags.load());
    EXPECT_GT(frags.capacity(), 127);
    EXPECT_EQ(frags.startData()[127], 500);
    EXPECT_GT(frags.endData()[127], 500);

    // Seek to high section of chr3
    frags.seek(2, 2126);
    EXPECT_TRUE(frags.load());
    EXPECT_GT(frags.startData()[0], 56);
    EXPECT_GT(frags.endData()[0], 2050);
    frags.seek(2, 0);
    EXPECT_TRUE(frags.load());
    EXPECT_EQ(frags.startData()[0], 10);
    // For chr2, because of the end_max leftover from chr1, seeking should go to start of chr
    frags.seek(1, 1126);
    EXPECT_TRUE(frags.load());
    EXPECT_EQ(frags.startData()[0], 10);
}

TEST(FragmentSeeking, RegionSelect) {
    using namespace Testing;

    // Strategy -- cell 0 has no overlaps, cell 1 has all overlaps. 
    // Edge cases to consider:
    // - Chromosome that is entirely covered by an overlap (chr1)
    // - Chromosome that has no regions in it (chr0)
    // - Overlapping regions (chr2)
    // - Having to match on chromosome names rather than IDs
    
    std::vector<std::string> region_chr_levels{"chr2", "chr1"};
    std::vector<uint32_t> region_chr{1, 0, 0, 0};
    std::vector<uint32_t> region_start{0, 10, 20, 50};
    std::vector<uint32_t> region_end{UINT32_MAX, 30, 40, 60};

    std::vector<Frag> c0, c1;
    // - Chromosome that has no regions in it (chr0)
    for(uint32_t i = 0; i < 20; i++) c0.push_back({0, i, i+10, 0});
    // - Chromosome that is entirely covered by an overlap (chr1)
    for(uint32_t i = 21; i < 40; i++) c1.push_back({1, i, i+10, 1});

    c0.push_back({2, 0, 10, 0});
    c0.push_back({2, 41, 50, 0});
    c0.push_back({2, 61, 62, 0});
    c1.push_back({2, 0, 41, 1});
    c1.push_back({2, 41, 61, 1});
    for (uint32_t i = 10; i <= 39; i++) {
        c1.push_back({2, 0, i+1, 1});
        c1.push_back({2, i, 42, 1});
    }
    for (uint32_t i = 50; i <= 59; i++) {
        c1.push_back({2, 42, i+1, 1});
        c1.push_back({2, i, 70, 1});
    }

    std::vector<Frag> both;
    both.insert(both.begin(), c0.begin(), c0.end());
    both.insert(both.begin(), c1.begin(), c1.end());

    std::unique_ptr<VecReaderWriterBuilder> v_both = writeFragmentTuple(both);
    StoredFragments frags_both = StoredFragments::openUnpacked(*v_both);

    std::unique_ptr<VecReaderWriterBuilder> v_c0 = writeFragmentTuple(c0, 2U);
    StoredFragments frags_c0 = StoredFragments::openUnpacked(*v_c0);

    std::unique_ptr<VecReaderWriterBuilder> v_c1 = writeFragmentTuple(c1);
    StoredFragments frags_c1 = StoredFragments::openUnpacked(*v_c1);
    
    std::unique_ptr<StringReader> chr_levels1 = std::make_unique<VecStringReader>(region_chr_levels);
    std::unique_ptr<StringReader> chr_levels2 = std::make_unique<VecStringReader>(region_chr_levels);
    RegionSelect inclusive(frags_both, region_chr, region_start, region_end, std::move(chr_levels1), false);
    RegionSelect exclusive(frags_both, region_chr, region_start, region_end, std::move(chr_levels2), true);


    EXPECT_TRUE(fragments_identical(frags_c1, inclusive));
    frags_both.restart();
    EXPECT_TRUE(fragments_identical(frags_c0, exclusive));
}