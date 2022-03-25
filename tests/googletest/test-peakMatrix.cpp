
#include <gtest/gtest.h>

#include <arrayIO/vector.h>
#include <fragmentIterators/FragmentIterator.h>
#include <fragmentIterators/StoredFragments.h>
#include <fragmentIterators/BedFragments.h>
#include <matrixIterators/PeakMatrix.h>
#include <matrixIterators/StoredMatrix.h>
#include <matrixIterators/CSparseMatrix.h>

#include <Eigen/SparseCore>


using namespace BPCells;
using namespace Eigen;

Map<SparseMatrix<double>> get_map(const SparseMatrix<double> &mat) {
   return Map<SparseMatrix<double>>(
        mat.rows(), 
        mat.cols(), 
        mat.nonZeros(), 
        (int *) mat.outerIndexPtr(), 
        (int *) mat.innerIndexPtr(), 
        (double *) mat.valuePtr()
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
    std::unique_ptr<StringWriter> w_chr_names =  v.createStringWriter("chr_names");
    std::unique_ptr<StringWriter> w_cell_names =  v.createStringWriter("cell_names");
    v.writeVersion("unpacked-fragments-v1");
    // Write chr1 test data
    uint32_t count = 0;
    // i = cell, j = start coord. The iteration order is to keep things start-sorted
    for (int j = 0; j < 5; j++) {
        for (int i = 0; i <= j; i++) {
            for (int k = 0; k < i + 1; k++) {
                w_cell.write_one(i);
                w_start.write_one(j);
                w_end.write_one(1002+i);    
                count += 1;    
            }
        }
    }
    w_chr_ptr.write_one(0); w_chr_ptr.write_one(count);
    // Write chr2 test data
    w_cell.write_one(0); w_start.write_one(9); w_end.write_one(21);
    w_cell.write_one(1); w_start.write_one(9); w_end.write_one(20);
    w_cell.write_one(2); w_start.write_one(10); w_end.write_one(21);
    w_cell.write_one(3); w_start.write_one(10); w_end.write_one(20);
    w_chr_ptr.write_one(count); w_chr_ptr.write_one(count+4);
    w_end_max.write_one(1001+4);

   
    w_chr_names->write(VecStringReader(std::vector<std::string>{"chr1", "chr2"}));
    w_cell_names->write(VecStringReader(std::vector<std::string>{"c0", "c1", "c2", "c3", "c4"}));
    

    w_cell.finalize(); w_start.finalize(); w_end.finalize(); w_end_max.finalize(); w_chr_ptr.finalize();
    StoredFragments frags = StoredFragments::openUnpacked(v);

    std::vector<uint32_t> chr = {0,0,0,1};
    std::vector<uint32_t> start = {2, 1002, 1004, 10};
    std::vector<uint32_t> end = {4, 1005, 1006, 20};
    PeakMatrix m(
        frags, chr, start, end, 
        std::make_unique<VecStringReader>(std::vector<std::string>{"chr1", "chr2"})
    );

    // MatrixIterator it(m);

    // while (it.nextCol()) {
    //     while (it.nextValue()) {
    //         printf("(%d, %d, %d)\n", it.row(), it.col(), it.val());
    //     }
    // }
    std::vector<Eigen::Triplet<double>> expected_triplets = {
        {0,0,2}, {1,0,4}, {2,0,6}, {3,0,4}, //Peak 1
        {1,1,8}, {2,1,9}, {3,1,8}, // Peak 2
        {3,2,8}, {4,2,5}, // Peak 3
        {1,3,1}, {2,3,1}, {3,3,2} // Peak 4
    };
    Eigen::SparseMatrix<double> eigen_mat(5,4);
    eigen_mat.setFromTriplets(expected_triplets.begin(), expected_triplets.end());
    
    CSparseMatrix expected_csparse(get_map(eigen_mat));
    MatrixConverterLoader<double, uint32_t> expected_int(expected_csparse);

    EXPECT_TRUE(matrix_identical_cpp(expected_int, m));
}



bool matrix_identical_cpp(MatrixLoader<uint32_t> &mat1, MatrixLoader<uint32_t> &mat2) {
    mat1.restart(); mat2.restart();
    MatrixIterator<uint32_t> i1(mat1);
    MatrixIterator<uint32_t> i2(mat2);

    while(true) {
        bool res1 = i1.nextCol();
        bool res2 = i2.nextCol();
        if(res1 != res2) {
            std::cerr << "Different number of columns." << std::endl;
            return false;
        }
        if (!res1) break;
        if(i1.currentCol() != i2.currentCol()) {
            std::cerr << "Different columnloaded" << std::endl;
            return false;
        }
        while(true) {
            bool res1 = i1.nextValue();
            bool res2 = i2.nextValue();
            if (res1 != res2) {
                std::cerr << "Different number of fragments in chromosome." << std::endl;
                return false;
            }
            if (!res1) break;
            if (i1.row() != i2.row() || i1.col() != i2.col() || i1.val() != i2.val()) {
                printf("Mismatched fragments: (%d,%d=%d) vs. (%d,%d=%d)\n",
                    i1.row(), i1.col(), i1.val(),
                    i2.row(), i2.col(), i2.val()
                );
                return false;
            }
        }
    }
    // TODO: Check row/col names
    return true;
}