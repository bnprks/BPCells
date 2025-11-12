#pragma once
#include "Pseudobulk.h"
#include <unordered_map>

namespace BPCells {

struct PseudobulkStatsSparse {
  Eigen::SparseMatrix<double> non_zeros;
  Eigen::SparseMatrix<double> sum;
  Eigen::SparseMatrix<double> mean;
  Eigen::SparseMatrix<double> var;
};

struct PseudobulkStatsTriplet {
  std::vector<Eigen::Triplet<double>> non_zeros;
  std::vector<Eigen::Triplet<double>> sum;
  std::vector<Eigen::Triplet<double>> mean;
  std::vector<Eigen::Triplet<double>> var;
};

struct PseudobulkStatsTemp {
  std::vector<double> non_zeros;
  std::vector<double> sum;
  std::vector<double> mean;
  std::vector<double> var;
};

template <typename T>
PseudobulkStatsSparse pseudobulk_matrix_sparse(std::unique_ptr<MatrixLoader<T>> &&mat,
                                               const std::vector<uint32_t>& cell_groups,
                                               PseudobulkStatsMethod method,
                                               bool transpose,
                                               std::atomic<bool> *user_interrupt);

} // namespace BPCells
