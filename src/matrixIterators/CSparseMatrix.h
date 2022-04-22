#include "MatrixIterator.h"

#ifndef RCPP_EIGEN
#include <Eigen/SparseCore>
#else
#include <RcppEigen.h>
#endif
// [[Rcpp::depends(RcppEigen)]]

namespace BPCells {

// Get Eigen and iterate over an Eigen sparse matrix
class CSparseMatrix : public MatrixLoader<double> {
    typedef Eigen::Map<Eigen::SparseMatrix<double>> EigenMat;
private:
    const EigenMat mat;
    std::vector<uint32_t> row_buf;
    std::vector<double>  val_buf;
    uint32_t idx;
    uint32_t load_size;
    uint32_t col;
public:
    CSparseMatrix(const EigenMat mat, uint32_t load_size = 1024) : mat(mat), load_size(load_size) {
        restart();
    }
    
    // Return the count of rows and columns
    uint32_t rows() const override { return mat.rows(); };
    uint32_t cols() const override { return mat.cols(); };

    // Return name for a given row or column.
    // If a matrix doesn't have assigned names this will return NULL
    const char* rowNames(uint32_t row) const override {
        return NULL;
    }
    const char* colNames(uint32_t col) const override {
        return NULL;
    }

    // Reset the iterator to start from the beginning
    void restart() override {
        col = UINT32_MAX;
    };

    void seekCol(uint32_t new_col) override {
        col = new_col;
        idx = mat.outerIndexPtr()[col];
        row_buf.resize(0);
        val_buf.resize(0);
    }

    bool nextCol() override {
        if (col+1 >= mat.cols()) return false;
        col++;
        idx = mat.outerIndexPtr()[col];
        row_buf.resize(0);
        val_buf.resize(0);
        return true;
    }

    uint32_t currentCol() const override {return col;}

    bool load() override {
        idx += capacity();
        if (idx == mat.outerIndexPtr()[col+1])  {
            row_buf.resize(0);
            val_buf.resize(0);
            return false;
        }
        uint32_t cap = std::min(load_size, mat.outerIndexPtr()[col+1] - idx);
        row_buf.resize(cap);
        val_buf.resize(cap);
        std::memmove(row_buf.data(), mat.innerIndexPtr() + idx, sizeof(uint32_t)*cap);
        std::memmove(val_buf.data(), mat.valuePtr() + idx, sizeof(double)*cap);
        return true;
    };

    uint32_t capacity() const override {return row_buf.size();}

    uint32_t* rowData() override {return row_buf.data();}
    double* valData() override {return val_buf.data();}

};


class CSparseMatrixWriter : public MatrixWriter<double> {
private:
    Eigen::SparseMatrix<double> eigen_mat;
public:
    void write(MatrixLoader<double> &loader, void (*checkInterrupt)(void) = NULL) override {
        MatrixIterator<double> mat(loader);
        uint32_t count = 0;
        std::vector<Eigen::Triplet<double>> triplets;
        
        while(mat.nextCol()) {
            while (mat.nextValue()) {
                triplets.push_back(Eigen::Triplet<double> (
                    mat.row(), mat.col(), mat.val()
                ));
                if(count++ % 8192 == 0 && checkInterrupt != NULL) checkInterrupt();
            }
        }
        
        eigen_mat = Eigen::SparseMatrix<double>(mat.rows(), mat.cols());
        eigen_mat.setFromTriplets(triplets.begin(), triplets.end());
        eigen_mat.makeCompressed();
    };

    const Eigen::SparseMatrix<double> getMat() {
        return eigen_mat;
    }
};

class CSparseTransposeMatrixWriter : public MatrixWriter<double> {
private:
    Eigen::SparseMatrix<double> eigen_mat;
public:
    void write(MatrixLoader<double> &loader, void (*checkInterrupt)(void) = NULL) override {
        MatrixIterator<double> mat(loader);
        uint32_t count = 0;
        std::vector<Eigen::Triplet<double>> triplets;
        
        while(mat.nextCol()) {
            while (mat.nextValue()) {
                triplets.push_back(Eigen::Triplet<double> (
                    mat.col(), mat.row(), mat.val()
                ));
                if(count++ % 8192 == 0 && checkInterrupt != NULL) checkInterrupt();
            }
        }
        
        eigen_mat = Eigen::SparseMatrix<double>(mat.cols(), mat.rows());
        eigen_mat.setFromTriplets(triplets.begin(), triplets.end());
        eigen_mat.makeCompressed();
    };

    const Eigen::SparseMatrix<double> getMat() {
        return eigen_mat;
    }
};

} // end namespace BPCells