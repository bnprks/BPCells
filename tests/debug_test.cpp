//
// Created by Ben Parks on 9/16/21.
//

#include <iostream>

#include "../src/matrixIterators/TileMatrix.h"
#include "../src/matrixIterators/TileMatrix2.h"
#include "../src/fragmentIterators/BedFragments.h"
#include "../src/fragmentIterators/CellSelect.h"
#include "../src/matrixIterators/CSparseMatrix.h"

int main() {
    BPCells::BedFragments bed(
            "/Users/ben/Downloads/PackedInsertionsData/01_raw_data/atac_pbmc_10k_2.0/fragments_chr1_chr10.tsv.gz",
            "#"
    );

    std::vector<uint32_t> cell_indices;
    for (int i = 0; i < 1000000; i++) cell_indices.push_back(i);

    BPCells::CellIndexSelect selected_cells(bed, cell_indices);

    std::vector<uint32_t> chr = {0, 1};
    std::vector<uint32_t> start = {0, 1};
    std::vector<uint32_t> end = {248956422, 133797422};
    std::vector<uint32_t> tile_width = {500, 500};
    std::vector<std::string> chr_names = {"chr1", "chr10"};

    BPCells::TileMatrix mat(selected_cells, chr, start, end, tile_width, chr_names);
    BPCells::MatrixConverterLoader<uint32_t, double> double_mat(mat);
    BPCells::CSparseMatrixWriter writer;
    writer.write(double_mat);
    printf("%ld\n", writer.getMat().nonZeros());
}