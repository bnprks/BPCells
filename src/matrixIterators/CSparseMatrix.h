#include "MatrixIterator.h"

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

namespace BPCells {

// Get Eigen and iterate over an Eigen sparse matrix
class CSparseMatrix : public MatrixLoader<double> {
    typedef Eigen::Map<Eigen::SparseMatrix<double>> EigenMat;
private:
    const EigenMat mat;
    EigenMat::InnerIterator colIter;
    uint32_t col;
public:
    CSparseMatrix(const EigenMat mat) : mat(mat) {
        restart();
    }
    
    // Reset the iterator to start from the beginning
    void restart() override {
        col = UINT32_MAX;
    };

    // Return the count of rows and columns
    uint32_t rows() const override { return mat.rows(); };
    uint32_t cols() const override { return mat.cols(); };

    // child class implementation of _nextCol. Should return
    // the index of the new column, or UINT32_MAX if there are no more columns
    bool nextCol() override {
        if (col+1 >= mat.cols()) return false;
        col++;
        colIter = EigenMat::InnerIterator(mat, col);
        return true;
    }

    uint32_t currentCol() const override {return col;}

    // Return number of matrix entries loaded. Should repeatedly return 0 at the end of the matrix
    // Return -1 for error
    int32_t load(uint32_t count, SparseVector<double> buffer) override {
        for(int i = 0; i < count; i++) {
            if(!colIter) {
                return i;
            }
            buffer.idx[i] = colIter.row();
            buffer.val[i] = colIter.value();
            ++colIter;
        }
        return count;
    };
};


class CSparseMatrixWriter : public MatrixWriter<double> {
private:
    Eigen::SparseMatrix<double> eigen_mat;
public:
    bool write(MatrixIterator<double> &mat, void (*checkInterrupt)(void) = NULL) override {
        uint32_t count = 0;
        std::vector<Eigen::Triplet<double>> triplets;
        
        while(mat.nextCol()) {
            while (mat.nextValue()) {
                triplets.push_back(Eigen::Triplet<double> (
                    mat.row(), mat.col(), mat.val()
                ));
                if(count++ % 10000 == 0 && checkInterrupt != NULL) checkInterrupt();
            }
        }
        
        eigen_mat = Eigen::SparseMatrix<double>(mat.rows(), mat.cols());
        eigen_mat.setFromTriplets(triplets.begin(), triplets.end());
        eigen_mat.makeCompressed();
        return true;
    };

    const Eigen::SparseMatrix<double> getMat() {
        return eigen_mat;
    }
};

} // end namespace BPCells