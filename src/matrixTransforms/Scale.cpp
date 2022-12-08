#include "Scale.h"

namespace BPCells {

// Scale rows and/or columns of a matrix
// out[i,j] = in[i,j] * row_params[0,i] * col_params[0,j]
// If row_params or col_params have 0 rows, then skip scaling along that dimension

Scale::Scale(MatrixLoader<double> &mat, TransformFit fit) : MatrixTransform(mat, fit) {}

bool Scale::load() {
    if (!loader.load()) return false;

    double *val_data = valData();
    const uint32_t *row_data = rowData();
    const uint32_t col = currentCol();
    const uint32_t cap = capacity();

    if (fit.col_params.size() > 0) {
        for (uint32_t i = 0; i < cap; i++) {
            val_data[i] *= fit.col_params(col);
        }
    }
    if (fit.row_params.size() > 0) {
        for (uint32_t i = 0; i < cap; i++) {
            val_data[i] *= fit.row_params(row_data[i]);
        }
    }

    return true;
}

Eigen::MatrixXd
Scale::denseMultiplyRight(const Eigen::Map<Eigen::MatrixXd> B, void (*checkInterrupt)(void)) {
    Eigen::MatrixXd res;

    // Scale input by col scale
    if (fit.col_params.size() > 0) {
        Eigen::MatrixXd B2(B.array().colwise() * fit.col_params.row(0).transpose());
        res = loader.denseMultiplyRight(
            Eigen::Map<Eigen::MatrixXd>(B2.data(), B2.rows(), B2.cols()), checkInterrupt
        );
    } else {
        res = loader.denseMultiplyRight(B, checkInterrupt);
    }

    // Scale output by row scale
    if (fit.row_params.size() > 0) {
        return res.array().colwise() * fit.row_params.row(0).transpose();
    }
    return res;
}

Eigen::MatrixXd
Scale::denseMultiplyLeft(const Eigen::Map<Eigen::MatrixXd> B, void (*checkInterrupt)(void)) {
    Eigen::MatrixXd res;

    // Scale input by row scale
    if (fit.row_params.size() > 0) {
        Eigen::MatrixXd B2(B.array().rowwise() * fit.row_params.row(0));
        res = loader.denseMultiplyLeft(
            Eigen::Map<Eigen::MatrixXd>(B2.data(), B2.rows(), B2.cols()), checkInterrupt
        );
    } else {
        res = loader.denseMultiplyLeft(B, checkInterrupt);
    }

    // Scale output by col scale
    if (fit.col_params.size() > 0) {
        return res.array().rowwise() * fit.col_params.row(0);
    }
    return res;
}
// Calculate matrix-vector product A*v where A (this) is sparse and B is a dense matrix.
Eigen::VectorXd
Scale::vecMultiplyRight(const Eigen::Map<Eigen::VectorXd> v, void (*checkInterrupt)(void)) {
    Eigen::VectorXd res;

    // Scale input by col scale
    if (fit.col_params.size() > 0) {
        Eigen::VectorXd v2(v.array() * fit.col_params.row(0).transpose());
        res = loader.vecMultiplyRight(
            Eigen::Map<Eigen::VectorXd>(v2.data(), v2.size()), checkInterrupt
        );
    } else {
        res = loader.vecMultiplyRight(v, checkInterrupt);
    }

    // Scale output by row scale
    if (fit.row_params.size() > 0) {
        return res.array() * fit.row_params.row(0).transpose();
    }
    return res;
}

Eigen::VectorXd
Scale::vecMultiplyLeft(const Eigen::Map<Eigen::VectorXd> v, void (*checkInterrupt)(void)) {
    Eigen::VectorXd res;

    // Scale input by row scale
    if (fit.row_params.size() > 0) {
        Eigen::VectorXd v2(v.array() * fit.row_params.row(0).transpose());
        res = loader.vecMultiplyLeft(
            Eigen::Map<Eigen::VectorXd>(v2.data(), v2.size()), checkInterrupt
        );
    } else {
        res = loader.vecMultiplyLeft(v, checkInterrupt);
    }

    // Scale output by col scale
    if (fit.col_params.size() > 0) {
        return res.array() * fit.col_params.row(0).transpose();
    }
    return res;
}

} // end namespace BPCells
