// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include <gtest/gtest.h>

#include "utils-fragments.h"

#include <array>
#include <arrayIO/vector.h>
#include <fragmentIterators/CellSelect.h>
#include <fragmentIterators/FragmentIterator.h>
#include <fragmentIterators/MergeFragments.h>
#include <fragmentIterators/RegionSelect.h>
#include <fragmentIterators/Rename.h>
#include <fragmentIterators/StoredFragments.h>
#include <fragmentUtils/InsertionIterator.h>

using namespace BPCells;

TEST(FragmentUtils, SeekStoredFrags) {
    using namespace Testing;
    std::vector<Frag> frags_vec;
    // Write chr1 test data, zig-zag up+down

    for (uint32_t start = 0; start < 2000; start += 500) {
        for (uint32_t end = start + 1; end < start + 250; end++) {
            frags_vec.push_back({0, start, end, 0});
        }
        for (uint32_t end = start + 250; end > start; end--) {
            frags_vec.push_back({0, start, end, 0});
        }
    }
    // Write chr2 test data, 128 frags not aligned on a 128-block boundary in
    // stored array Smaller ends than previous chr
    for (uint32_t i = 0; i < 128; i++) {
        frags_vec.push_back({1, i + 10, i + 1000, 0});
    }

    // Write chr2 test data, 128 frags not aligned on a 128-block boundary in
    // stored array Larger ends than previous chr
    for (uint32_t i = 0; i < 128; i++) {
        frags_vec.push_back({2, i + 10, i + 2000, 0});
    }

    std::unique_ptr<VecReaderWriterBuilder> v = writeFragmentTuple(frags_vec);
    StoredFragments frags = StoredFragments::openUnpacked(*v);

    // Try seeking a couple locations in chr1 to make sure basic functionality
    // works
    frags.seek(0, 1249);
    EXPECT_TRUE(frags.load());
    EXPECT_EQ(frags.startData()[0], 1000);
    EXPECT_GT(frags.endData()[0], 1249 - 128);
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
    // For chr2, because of the end_max leftover from chr1, seeking should go to
    // start of chr
    frags.seek(1, 1126);
    EXPECT_TRUE(frags.load());
    EXPECT_EQ(frags.startData()[0], 10);
}

TEST(FragmentUtils, RegionSelect) {
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
    for (uint32_t i = 0; i < 20; i++)
        c0.push_back({0, i, i + 10, 0});
    // - Chromosome that is entirely covered by an overlap (chr1)
    for (uint32_t i = 21; i < 40; i++)
        c1.push_back({1, i, i + 10, 1});

    c0.push_back({2, 0, 10, 0});
    c0.push_back({2, 41, 50, 0});
    c0.push_back({2, 61, 62, 0});
    c1.push_back({2, 0, 41, 1});
    c1.push_back({2, 41, 61, 1});
    for (uint32_t i = 10; i <= 39; i++) {
        c1.push_back({2, 0, i + 1, 1});
        c1.push_back({2, i, 42, 1});
    }
    for (uint32_t i = 50; i <= 59; i++) {
        c1.push_back({2, 42, i + 1, 1});
        c1.push_back({2, i, 70, 1});
    }

    std::vector<Frag> both;
    both.insert(both.begin(), c0.begin(), c0.end());
    both.insert(both.begin(), c1.begin(), c1.end());

    std::unique_ptr<VecReaderWriterBuilder> v_both = writeFragmentTuple(both);

    std::unique_ptr<VecReaderWriterBuilder> v_c0 = writeFragmentTuple(c0, 2U);
    StoredFragments frags_c0 = StoredFragments::openUnpacked(*v_c0);

    std::unique_ptr<VecReaderWriterBuilder> v_c1 = writeFragmentTuple(c1);
    StoredFragments frags_c1 = StoredFragments::openUnpacked(*v_c1);

    std::unique_ptr<StringReader> chr_levels1 =
        std::make_unique<VecStringReader>(region_chr_levels);
    std::unique_ptr<StringReader> chr_levels2 =
        std::make_unique<VecStringReader>(region_chr_levels);
    RegionSelect inclusive(
        std::make_unique<StoredFragments>(StoredFragments::openUnpacked(*v_both)),
        region_chr,
        region_start,
        region_end,
        std::move(chr_levels1),
        false
    );
    RegionSelect exclusive(
        std::make_unique<StoredFragments>(StoredFragments::openUnpacked(*v_both)),
        region_chr,
        region_start,
        region_end,
        std::move(chr_levels2),
        true
    );

    EXPECT_TRUE(fragments_identical(frags_c1, inclusive));
    EXPECT_TRUE(fragments_identical(frags_c0, exclusive));
}

TEST(FragmentUtils, MergeFragments) {
    uint32_t max_cell = 50;
    auto v1 = Testing::generateFrags(1000, 3, 200, max_cell - 1, 25, 1336);
    auto v2 = Testing::generateFrags(1000, 3, 200, max_cell - 1, 25, 1334);
    auto v3 = Testing::generateFrags(1000, 3, 200, max_cell - 1, 25, 1227);

    // Make sure that there are no matching start coordinates between v1, v2, v3
    // since the ordering is undefined 
    for (auto &f : v1) {
        uint32_t width = f.end - f.start;
        f.start = f.start * 3;
        f.end = f.start + width;
    }
    for (auto &f : v2) {
        uint32_t width = f.end - f.start;
        f.start = f.start * 3 + 1;
        f.end = f.start + width;
    }
    for (auto &f : v3) {
        uint32_t width = f.end - f.start;
        f.start = f.start * 3 + 2;
        f.end = f.start + width;
    }

    std::sort(v1.begin(), v1.end(), [](const Testing::Frag &a, const Testing::Frag &b) {
        if (a.chr != b.chr) return a.chr < b.chr;
        return a.start < b.start;
    });
    std::sort(v2.begin(), v2.end(), [](const Testing::Frag &a, const Testing::Frag &b) {
        if (a.chr != b.chr) return a.chr < b.chr;
        return a.start < b.start;
    });
    std::sort(v3.begin(), v3.end(), [](const Testing::Frag &a, const Testing::Frag &b) {
        if (a.chr != b.chr) return a.chr < b.chr;
        return a.start < b.start;
    });

    std::vector<Testing::Frag> v;
    v.insert(v.end(), v1.begin(), v1.end());
    v.insert(v.end(), v2.begin(), v2.end());
    v.insert(v.end(), v3.begin(), v3.end());
    uint32_t idx = 0;
    for (int i = 0; i < v.size(); i++) {
        v[i].cell += max_cell * (i / 1000);
    }
    std::stable_sort(v.begin(), v.end(), [](const Testing::Frag &a, const Testing::Frag &b) {
        if (a.chr != b.chr) return a.chr < b.chr;
        return a.start < b.start;
    });

    std::unique_ptr<VecReaderWriterBuilder> v_expect = writeFragmentTuple(v, max_cell, true);
    StoredFragments expected = StoredFragments::openUnpacked(*v_expect);

    std::unique_ptr<VecReaderWriterBuilder> v1_data = writeFragmentTuple(v1, max_cell, true);

    std::unique_ptr<VecReaderWriterBuilder> v2_data = writeFragmentTuple(v2, max_cell, true);
    std::vector<std::string> &names_v2 = v2_data->getStringVecs().at("cell_names");
    for (uint32_t i = 0; i < max_cell; i++) {
        names_v2[i] = std::string("c") + std::to_string(i + max_cell);
    }

    std::unique_ptr<VecReaderWriterBuilder> v3_data = writeFragmentTuple(v3, max_cell, true);
    std::vector<std::string> &names_v3 = v3_data->getStringVecs().at("cell_names");
    for (uint32_t i = 0; i < max_cell; i++) {
        names_v3[i] = std::string("c") + std::to_string(i + 2 * max_cell);
    }

    std::vector<std::unique_ptr<FragmentLoader>> merge_vec;
    merge_vec.push_back(std::make_unique<StoredFragments>(StoredFragments::openUnpacked(*v1_data)));
    merge_vec.push_back(std::make_unique<StoredFragments>(StoredFragments::openUnpacked(*v2_data)));
    merge_vec.push_back(std::make_unique<StoredFragments>(StoredFragments::openUnpacked(*v3_data)));

    MergeFragments merge(std::move(merge_vec), v_expect->getStringVecs().at("chr_names"));

    EXPECT_TRUE(Testing::fragments_identical(expected, merge));
}

TEST(FragmentUtils, InsertionIterator) {
    uint32_t max_cell = 50;
    auto v = Testing::generateFrags(2000, 3, 400, max_cell - 1, 100, 1336);
    std::sort(v.begin(), v.end(), [](const Testing::Frag &a, const Testing::Frag &b) {
        if (a.chr != b.chr) return a.chr < b.chr;
        return a.start < b.start;
    });
    std::unique_ptr<VecReaderWriterBuilder> d = writeFragmentTuple(v, max_cell, true, 50);
    StoredFragments frags = StoredFragments::openUnpacked(*d);

    InsertionIterator it(frags);

    std::vector<std::array<uint32_t, 3>> insert;

    for (const auto &f : v)
        insert.push_back({f.chr, f.start, f.cell});
    for (const auto &f : v)
        insert.push_back({f.chr, f.end - 1, f.cell});

    std::stable_sort(
        insert.begin(),
        insert.end(),
        [](const std::array<uint32_t, 3> &a, const std::array<uint32_t, 3> &b) {
            if (a[0] != b[0]) return a[0] < b[0];
            return a[1] < b[1];
        }
    );

    uint32_t current_chr = 0;
    ASSERT_TRUE(it.nextChr());
    for (int i = 0; i < insert.size(); i++) {
        if (insert[i][0] != current_chr) {
            ASSERT_FALSE(it.nextInsertion());
            ASSERT_FALSE(it.nextInsertion());
            ASSERT_TRUE(it.nextChr());
            current_chr = insert[i][0];
        }
        ASSERT_TRUE(it.nextInsertion());
        ASSERT_EQ(it.chr(), current_chr);
        ASSERT_EQ(it.coord(), insert[i][1]);
        ASSERT_EQ(it.cell(), insert[i][2]);
    }
    ASSERT_FALSE(it.nextInsertion());
    ASSERT_FALSE(it.nextChr());
    ASSERT_FALSE(it.nextChr());
}

TEST(FragmentUtils, CellSelect) {
    uint32_t max_cell = 50;
    auto v = Testing::generateFrags(200, 3, 400, max_cell - 1, 100, 1336);
    std::sort(v.begin(), v.end(), [](const Testing::Frag &a, const Testing::Frag &b) {
        if (a.chr != b.chr) return a.chr < b.chr;
        return a.start < b.start;
    });

    std::vector<Testing::Frag> v2;
    for (auto f : v) {
        if (f.cell > 15 || f.cell < 10) continue;
        f.cell = 15 - f.cell;
        v2.push_back(f);
    }

    std::unique_ptr<VecReaderWriterBuilder> d1 = writeFragmentTuple(v, max_cell, true);
    std::unique_ptr<VecReaderWriterBuilder> d2 = writeFragmentTuple(v2, max_cell, true);
    std::vector<std::string> &names_d2 = d2->getStringVecs().at("cell_names");
    for (uint32_t i = 0; i <= 5; i++) {
        names_d2[i] = std::string("c") + std::to_string(15 - i);
    }
    names_d2.resize(6);

    StoredFragments out = StoredFragments::openUnpacked(*d2);
    CellNameSelect select1(
        std::make_unique<StoredFragments>(StoredFragments::openUnpacked(*d1)),
        {"c15", "c14", "c13", "c12", "c11", "c10"}
    );

    ASSERT_TRUE(Testing::fragments_identical(select1, out));
    CellIndexSelect select2(
        std::make_unique<StoredFragments>(StoredFragments::openUnpacked(*d1)),
        {15, 14, 13, 12, 11, 10}
    );
    out.restart();
    ASSERT_TRUE(Testing::fragments_identical(select2, out));
}

TEST(FragmentUtils, CellPrefix) {
    uint32_t max_cell = 50;
    auto v = Testing::generateFrags(200, 3, 400, max_cell - 1, 100, 1336);

    std::unique_ptr<VecReaderWriterBuilder> d1 = writeFragmentTuple(v, max_cell, true);
    std::vector<std::string> &names_d1 = d1->getStringVecs().at("cell_names");
    for (uint32_t i = 0; i <= 5; i++) {
        names_d1[i] = std::string("c") + std::to_string(15 - i);
    }
    names_d1.resize(6);

    StoredFragments in1 = StoredFragments::openUnpacked(*d1);
    PrefixCells in2(
        std::make_unique<StoredFragments>(StoredFragments::openUnpacked(*d1)), "prefix#"
    );
    PrefixCells in3(std::make_unique<StoredFragments>(StoredFragments::openUnpacked(*d1)), "");

    uint32_t i = 0;
    while (in1.cellNames(i) != NULL) {
        ASSERT_EQ(std::string(in1.cellNames(i)), std::string("c") + std::to_string(15 - i));
        ASSERT_EQ("prefix#" + std::string(in1.cellNames(i)), std::string(in2.cellNames(i)));
        ASSERT_EQ(std::string(in1.cellNames(i)), std::string(in3.cellNames(i)));
        i += 1;
    }
}