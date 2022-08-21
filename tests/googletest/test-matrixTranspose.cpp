#include <filesystem>
#include <random>

#include <gtest/gtest.h>

#include <matrixIterators/StoredMatrixTransposeWriter.h>
#include <matrixIterators/StoredMatrixWriter.h>
#include <matrixIterators/CSparseMatrix.h>
#include <matrixIterators/ImportMatrixHDF5.h>
#include <arrayIO/vector.h>

#include <Eigen/Core>

namespace fs = std::filesystem;
using namespace BPCells;
using namespace ::testing;
using namespace Eigen;

SparseMatrix<double> generate_mat(uint32_t n_row, uint32_t n_col, uint32_t seed=125124) {
    std::mt19937 gen(seed); //Standard mersenne_twister_engine
    std::uniform_int_distribution<> distrib(1, 20);
    std::uniform_int_distribution<> nonzero(0, 4); // 1/5 chance of being non-zero
 
    std::vector<Triplet<double>> triplets;

    for (int i = 0; i < n_row; i++) {
        for (int j = 0; j < n_col; j++) {
            if (0 == nonzero(gen)) {
                triplets.push_back({i, j, distrib(gen)});
            }
        }
    }
    
    SparseMatrix<double> mat(n_row,n_col);
    mat.setFromTriplets(triplets.begin(), triplets.end());
    return mat;
}

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


void test_transpose(const Eigen::SparseMatrix<double> orig_mat) {
    CSparseMatrix mat_d (get_map(orig_mat));
    MatrixConverterLoader<double, uint32_t> mat_i(mat_d);

    VecReaderWriterBuilder vb1(1024);

    fs::remove_all(fs::temp_directory_path() / "tmp_storage_uint");
    fs::remove_all(fs::temp_directory_path() / "tmp_storage_double");
    // Use small load sizes to help boost the number of rounds used for merging
    StoredMatrixTransposeWriter<uint32_t> w_uint((fs::temp_directory_path() / "tmp_storage_uint").c_str(), 512, 16384);
    StoredMatrixTransposeWriter<double> w_double((fs::temp_directory_path() / "tmp_storage_double").c_str(), 512, 16384);

    w_double.write(mat_d);
    mat_d.restart();
    w_uint.write(mat_i);

    StoredMatrix<uint32_t> trans_i = w_uint.read();
    auto loader1 = MatrixConverterLoader<uint32_t, double>(trans_i);

    StoredMatrix<double> loader2 = w_double.read();
    
    CSparseMatrixWriter mem1, mem2;
    mem1.write(loader1);
    mem2.write(loader2);
    EXPECT_TRUE(mem1.getMat().isApprox(orig_mat.transpose()));
    EXPECT_TRUE(mem2.getMat().isApprox(orig_mat.transpose()));
}

TEST(MatrixTranspose, SmallIntMatrix) {
    test_transpose(generate_mat(100, 2000));
    test_transpose(generate_mat(3, 100000));
    test_transpose(generate_mat(100000, 3));
    test_transpose(generate_mat(2000, 100));
}

// TEST(MatrixTranspose, MatFile) {
//     auto mat = open10xFeatureMatrix("/Users/ben/Downloads/20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5", 16384);
//     std::filesystem::remove_all(std::filesystem::path("test-dir-deleteme"));
//     StoredMatrixTransposeWriter<uint32_t> mat_t("test-dir-deleteme", 4194304, 1073741824);
//     mat_t.write(mat);
//     auto mat_t_read = mat_t.read();
//     VecReaderWriterBuilder data;
//     StoredMatrixWriter<uint32_t>::createPacked(data).write(mat_t_read);
// }