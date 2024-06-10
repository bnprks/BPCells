// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include <gtest/gtest.h>

#include <arrayIO/vector.h>
#include <fragmentIterators/BedFragments.h>
#include <fragmentIterators/FragmentIterator.h>
#include <fragmentIterators/StoredFragments.h>
#include <matrixIterators/CSparseMatrix.h>
#include <matrixIterators/PeakMatrix.h>
#include <matrixIterators/StoredMatrix.h>
#include <matrixIterators/TileMatrix.h>

#include <Eigen/SparseCore>

#include "utils-fragments.h"

using namespace BPCells;
using namespace Eigen;

Map<SparseMatrix<double>> get_map(const SparseMatrix<double> &mat) {
    return Map<SparseMatrix<double>>(
        mat.rows(),
        mat.cols(),
        mat.nonZeros(),
        (int *)mat.outerIndexPtr(),
        (int *)mat.innerIndexPtr(),
        (double *)mat.valuePtr()
    );
}

bool matrix_identical_cpp(MatrixLoader<uint32_t> &mat1, MatrixLoader<uint32_t> &mat2);

TEST(PeakMatrix, PeakMatrix) {
    VecReaderWriterBuilder v;

    UIntWriter w_cell = v.createUIntWriter("cell");
    UIntWriter w_start = v.createUIntWriter("start");
    UIntWriter w_end = v.createUIntWriter("end");
    UIntWriter w_end_max = v.createUIntWriter("end_max");
    UIntWriter w_chr_ptr = v.createUIntWriter("chr_ptr");
    std::unique_ptr<StringWriter> w_chr_names = v.createStringWriter("chr_names");
    std::unique_ptr<StringWriter> w_cell_names = v.createStringWriter("cell_names");
    v.writeVersion("unpacked-fragments-v1");
    // Write chr1 test data
    uint32_t count = 0;
    // i = cell, j = start coord. The iteration order is to keep things
    // start-sorted
    for (int j = 0; j < 5; j++) {
        for (int i = 0; i <= j; i++) {
            for (int k = 0; k < i + 1; k++) {
                w_cell.write_one(i);
                w_start.write_one(j);
                w_end.write_one(1002 + i);
                count += 1;
            }
        }
    }
    w_chr_ptr.write_one(0);
    w_chr_ptr.write_one(count);
    // Write chr2 test data
    w_cell.write_one(0);
    w_start.write_one(9);
    w_end.write_one(21);
    w_cell.write_one(1);
    w_start.write_one(9);
    w_end.write_one(20);
    w_cell.write_one(2);
    w_start.write_one(10);
    w_end.write_one(21);
    w_cell.write_one(3);
    w_start.write_one(10);
    w_end.write_one(20);
    w_chr_ptr.write_one(count);
    w_chr_ptr.write_one(count + 4);
    w_end_max.write_one(1001 + 4);

    w_chr_names->write(VecStringReader(std::vector<std::string>{"chr1", "chr2"}));
    w_cell_names->write(VecStringReader(std::vector<std::string>{"c0", "c1", "c2", "c3", "c4"}));

    w_cell.finalize();
    w_start.finalize();
    w_end.finalize();
    w_end_max.finalize();
    w_chr_ptr.finalize();

    std::vector<uint32_t> chr = {0, 0, 0, 0, 0, 1};
    std::vector<uint32_t> start = {2, 1002, 1004, 1000, 2000, 10};
    std::vector<uint32_t> end = {4, 1005, 1006, 1008, 2010, 20};
    PeakMatrix m(
        std::make_unique<StoredFragments>(StoredFragments::openUnpacked(v)),
        chr,
        start,
        end,
        std::make_unique<VecStringReader>(std::vector<std::string>{"chr1", "chr2"}),
        0 // insertion overlaps
    );

    // MatrixIterator it(m);

    // while (it.nextCol()) {
    //     while (it.nextValue()) {
    //         printf("(%d, %d, %d)\n", it.row(), it.col(), it.val());
    //     }
    // }
    std::vector<Eigen::Triplet<double>> expected_triplets = {
        {0, 0, 2},
        {1, 0, 4},
        {2, 0, 6},
        {3, 0, 4}, // Peak 1
        {1, 1, 8},
        {2, 1, 9},
        {3, 1, 8}, // Peak 2
        {3, 2, 8},
        {4, 2, 5}, // Peak 3
        {0, 3, 5},
        {1, 3, 8},
        {2, 3, 9},
        {3, 3, 8},
        {4, 3, 5}, // Peak 4
        {1, 5, 1},
        {2, 5, 1},
        {3, 5, 2} // Peak 6
    };
    Eigen::SparseMatrix<double> eigen_mat(5, 6);
    eigen_mat.setFromTriplets(expected_triplets.begin(), expected_triplets.end());

    MatrixConverterLoader<double, uint32_t> expected_int(
        std::make_unique<CSparseMatrix<double>>(get_map(eigen_mat))
    );

    EXPECT_TRUE(matrix_identical_cpp(expected_int, m));
}

TEST(PeakMatrix, PeakMatrixSeek) {
    uint32_t chrs = 5;
    uint32_t max_coord = 1000;
    auto v = Testing::writeFragmentTuple(Testing::generateFrags(50000, chrs, max_coord, 50, 125));

    std::vector<uint32_t> v_chr, v_start, v_end;
    std::vector<std::string> chr_levels;
    for (uint32_t chr = 0; chr <= chrs; chr++) {
        chr_levels.push_back(std::string("chr") + std::to_string(chr));
        for (uint32_t start = 100; start < max_coord; start += 100) {
            v_chr.push_back(chr);
            v_start.push_back(start);
            v_end.push_back(start + 25);
        }
    }

    auto peak_mat = std::make_unique<PeakMatrix>(
        std::make_unique<StoredFragments>(StoredFragments::openUnpacked(*v)),
        v_chr,
        v_start,
        v_end,
        std::make_unique<VecStringReader>(chr_levels),
        0 // Count insertions
    );
    auto peak_mat_double =
        std::make_unique<MatrixConverterLoader<uint32_t, double>>(std::move(peak_mat));
    CSparseMatrixWriter<double> mat_writer;
    mat_writer.write(*peak_mat_double);

    Eigen::SparseMatrix<double> eigen_mat = mat_writer.getMat();

    MatrixIterator<double> peak_mat_it(std::move(peak_mat_double));
    peak_mat_it.restart();
    // Check that the first column upon restart matches the sparse matrix first
    // column
    ASSERT_TRUE(peak_mat_it.nextCol());
    bool has_nonzero = false;
    for (SparseMatrix<double>::InnerIterator it(eigen_mat, 0); it; ++it) {
        has_nonzero = true;
        ASSERT_TRUE(peak_mat_it.nextValue());
        ASSERT_EQ(peak_mat_it.row(), it.row());
        ASSERT_EQ(peak_mat_it.val(), it.value());
    }
    EXPECT_FALSE(peak_mat_it.nextValue());
    EXPECT_TRUE(has_nonzero);

    std::mt19937 gen(1337);
    std::uniform_int_distribution col(0, (int)v_chr.size() - 1);
    // Check that 50 random col seeks all work
    for (int i = 0; i < 50; i++) {
        uint32_t c = col(gen);
        bool has_nonzero = false;
        peak_mat_it.seekCol(c);
        for (SparseMatrix<double>::InnerIterator it(eigen_mat, c); it; ++it) {
            has_nonzero = true;
            ASSERT_TRUE(peak_mat_it.nextValue());
            ASSERT_EQ(peak_mat_it.row(), it.row());
            ASSERT_EQ(peak_mat_it.val(), it.value());
        }
        EXPECT_FALSE(peak_mat_it.nextValue());
        ASSERT_TRUE(has_nonzero);
    }
}

TEST(PeakMatrix, TileMatrix) {

    // Cases to test:
    // - Reads that span more than one tile region
    // - Correctly truncate the last tile for regions that aren't an even multiple
    // of tile width
    // - Handle reads with and without overlaps for a tile

    std::vector<uint32_t> chr = {0, 0, 0, 1};
    std::vector<uint32_t> start = {10, 30, 50, 70};
    std::vector<uint32_t> end = {20, 40, 60, 80};
    std::vector<uint32_t> width = {5, 3, 5, 12};

    VecReaderWriterBuilder v;
    UIntWriter w_cell = v.createUIntWriter("cell");
    UIntWriter w_start = v.createUIntWriter("start");
    UIntWriter w_end = v.createUIntWriter("end");
    UIntWriter w_end_max = v.createUIntWriter("end_max");
    UIntWriter w_chr_ptr = v.createUIntWriter("chr_ptr");
    std::unique_ptr<StringWriter> w_chr_names = v.createStringWriter("chr_names");
    std::unique_ptr<StringWriter> w_cell_names = v.createStringWriter("cell_names");
    v.writeVersion("unpacked-fragments-v1");

    uint32_t count = 0;

    // No overlaps on cell 0
    w_cell.write_one(0);
    w_start.write_one(9);
    w_end.write_one(21);
    count += 1;
    w_cell.write_one(0);
    w_start.write_one(9);
    w_end.write_one(10);
    count += 1;
    // Overlap spanning regions on cell 1
    w_cell.write_one(1);
    w_start.write_one(12);
    w_end.write_one(78);
    count += 1;
    // Tile middle region by end coord on cell 2
    for (int i = 0; i < 12; i++) {
        for (int j = 0; j <= i; j++) {
            w_cell.write_one(2);
            w_start.write_one(11 + i);
            w_end.write_one(30 + i);
            count += 1;
        }
    }
    // More no overlaps on cell 0
    w_cell.write_one(0);
    w_start.write_one(20);
    w_end.write_one(21);
    count += 1;
    // Tile middle region by start coord on cell 3
    for (int i = 0; i < 12; i++) {
        for (int j = 0; j <= i + 1; j++) {
            w_cell.write_one(3);
            w_start.write_one(29 + i);
            w_end.write_one(50 + i);
            count += 1;
        }
    }
    w_chr_ptr.write_one(0);
    w_chr_ptr.write_one(count);
    // Write chr2 test data
    w_cell.write_one(0);
    w_start.write_one(69);
    w_end.write_one(81);
    w_cell.write_one(1);
    w_start.write_one(69);
    w_end.write_one(80);
    w_cell.write_one(2);
    w_start.write_one(70);
    w_end.write_one(81);
    w_cell.write_one(3);
    w_start.write_one(70);
    w_end.write_one(80);
    w_chr_ptr.write_one(count);
    w_chr_ptr.write_one(count + 4);
    for (int i = 0; i < count; i += 128) {
        w_end_max.write_one(81);
    }

    w_chr_names->write(VecStringReader(std::vector<std::string>{"chr1", "chr2"}));
    w_cell_names->write(VecStringReader(std::vector<std::string>{"c0", "c1", "c2", "c3", "c4"}));

    w_cell.finalize();
    w_start.finalize();
    w_end.finalize();
    w_end_max.finalize();
    w_chr_ptr.finalize();
    StoredFragments frags = StoredFragments::openUnpacked(v);

    TileMatrix m(
        std::make_unique<StoredFragments>(StoredFragments::openUnpacked(v)),
        chr,
        start,
        end,
        width,
        std::make_unique<VecStringReader>(std::vector<std::string>{"chr1", "chr2"})
    );

    // MatrixIterator it(m);
    // print_matrix(it);
    // while (it.nextCol()) {
    //     while (it.nextValue()) {
    //         printf("(%d, %d, %d)\n", it.row(), it.col(), it.val());
    //     }
    // }

    std::vector<Eigen::Triplet<double>> expected_triplets = {
        {1, 0, 1},
        {2, 0, 10}, // Tile 1.1
        {2, 1, 35}, // Tile 1.2
        {2, 2, 9},
        {3, 2, 12}, // Tile 2.1
        {2, 3, 18},
        {3, 3, 21}, // Tile 2.2
        {2, 4, 27},
        {3, 4, 30}, // Tile 2.3
        {2, 5, 11},
        {3, 5, 12}, // Tile 2.4
        {3, 6, 25}, // Tile 3.1
        {3, 7, 50}, // Tile 3.2
        {1, 8, 1},
        {2, 8, 1},
        {3, 8, 2} // Tile 4.1
    };
    Eigen::SparseMatrix<double> eigen_mat(5, 9);
    eigen_mat.setFromTriplets(expected_triplets.begin(), expected_triplets.end());

    MatrixConverterLoader<double, uint32_t> expected_int(
        std::make_unique<CSparseMatrix<double>>(get_map(eigen_mat))
    );

    EXPECT_TRUE(matrix_identical_cpp(expected_int, m));
}

TEST(PeakMatrix, TileMatrixSeek) {
    uint32_t chrs = 5;
    uint32_t max_coord = 1000;
    auto v = Testing::writeFragmentTuple(Testing::generateFrags(50000, chrs, max_coord, 50, 125));

    std::vector<uint32_t> v_chr, v_start, v_end, v_width;
    std::vector<std::string> chr_levels;
    // Use some slightly variable tile choices, 1 tile per chromosome
    for (uint32_t chr = 0; chr <= chrs; chr++) {
        v_chr.push_back(chr);
        v_start.push_back(chr);
        v_end.push_back(max_coord - chr);
        v_width.push_back(25 * (chr + 1));
        chr_levels.push_back(std::string("chr") + std::to_string(chr));
    }

    auto tile_mat = std::make_unique<TileMatrix>(
        std::make_unique<StoredFragments>(StoredFragments::openUnpacked(*v)),
        v_chr,
        v_start,
        v_end,
        v_width,
        std::make_unique<VecStringReader>(chr_levels)
    );
    std::unique_ptr<MatrixLoader<double>> tile_mat_double =
        std::make_unique<MatrixConverterLoader<uint32_t, double>>(std::move(tile_mat));
    CSparseMatrixWriter<double> mat_writer;

    mat_writer.write(*tile_mat_double);
    Eigen::SparseMatrix<double> eigen_mat = mat_writer.getMat();

    MatrixIterator<double> tile_mat_it(std::move(tile_mat_double));
    tile_mat_it.restart();
    // Check that the first column upon restart matches the sparse matrix first
    // column
    ASSERT_TRUE(tile_mat_it.nextCol());
    bool has_nonzero = false;
    for (SparseMatrix<double>::InnerIterator it(eigen_mat, 0); it; ++it) {
        has_nonzero = true;
        ASSERT_TRUE(tile_mat_it.nextValue());
        ASSERT_EQ(tile_mat_it.row(), it.row());
        ASSERT_EQ(tile_mat_it.val(), it.value());
    }
    EXPECT_FALSE(tile_mat_it.nextValue());
    EXPECT_TRUE(has_nonzero);

    std::mt19937 gen(1337);
    std::uniform_int_distribution col(0, (int)tile_mat_it.cols() - 1);

    // Check that 50 random col seeks all work
    for (int i = 0; i < 50; i++) {
        uint32_t c = col(gen);
        bool has_nonzero = false;
        tile_mat_it.seekCol(c);
        for (SparseMatrix<double>::InnerIterator it(eigen_mat, c); it; ++it) {
            has_nonzero = true;
            ASSERT_TRUE(tile_mat_it.nextValue());
            ASSERT_EQ(tile_mat_it.row(), it.row());
            ASSERT_EQ(tile_mat_it.val(), it.value());
        }
        EXPECT_FALSE(tile_mat_it.nextValue());
        ASSERT_TRUE(has_nonzero);
    }
}

bool matrix_identical_cpp(MatrixLoader<uint32_t> &mat1, MatrixLoader<uint32_t> &mat2) {
    mat1.restart();
    mat2.restart();
    MatrixIterator<uint32_t> i1((std::unique_ptr<MatrixLoader<uint32_t>>(&mat1)));
    MatrixIterator<uint32_t> i2((std::unique_ptr<MatrixLoader<uint32_t>>(&mat2)));
    i1.preserve_input_loader();
    i2.preserve_input_loader();

    while (true) {
        bool res1 = i1.nextCol();
        bool res2 = i2.nextCol();
        if (res1 != res2) {
            std::cerr << "Different number of columns." << std::endl;
            return false;
        }
        if (!res1) break;
        if (i1.currentCol() != i2.currentCol()) {
            std::cerr << "Different columnloaded" << std::endl;
            return false;
        }
        while (true) {
            bool res1 = i1.nextValue();
            bool res2 = i2.nextValue();
            if (res1 != res2) {
                std::cerr << "Different number of entries in column." << std::endl;
                return false;
            }
            if (!res1) break;
            if (i1.row() != i2.row() || i1.col() != i2.col() || i1.val() != i2.val()) {
                printf(
                    "Mismatched entries: (%d,%d=%d) vs. (%d,%d=%d)\n",
                    i1.row(),
                    i1.col(),
                    i1.val(),
                    i2.row(),
                    i2.col(),
                    i2.val()
                );
                return false;
            }
        }
    }
    // TODO: Check row/col names
    return true;
}