
#include "PseudobulkSparse.h"
namespace BPCells {

template <typename T>
PseudobulkStatsSparse pseudobulk_matrix_sparse(std::unique_ptr<MatrixLoader<T>> &&mat,
                                               const std::vector<uint32_t>& cell_groups,
                                               PseudobulkStatsMethod method,
                                               bool transpose,
                                               std::atomic<bool> *user_interrupt) {
  MatrixIterator<T> it(std::move(mat));
  if (transpose && (it.rows() != cell_groups.size())) {
    throw std::invalid_argument("pseudobulk_matrix_sparse(): Cell groups must match the number of columns in the matrix");
  } else if (!transpose && (it.cols() != cell_groups.size())) {
    throw std::invalid_argument("pseudobulk_matrix_sparse(): Cell groups must match the number of columns in the matrix");
  }
  // Count each group
  std::vector<double> group_count_vec;
  uint32_t last_cell = 0;
  for (auto& cell_group : cell_groups) {
    if (cell_group < last_cell) {
      throw std::invalid_argument("pseudobulk_matrix_sparse(): Cell groups must be ordered");
    }
    last_cell = cell_group;
    if (cell_group >= group_count_vec.size()) group_count_vec.resize(cell_group + 1);
    group_count_vec[cell_group]++;
  }
  Eigen::Map<Eigen::Array<double, 1, Eigen::Dynamic>> group_count(group_count_vec.data(), group_count_vec.size());
  // If transposed, also keep matrix transposed with groups as rows for cache-friendliness
  uint32_t num_rows = transpose ? it.cols() : it.rows();
  uint32_t num_cols = group_count_vec.size();
  uint32_t tmp_size;
  if (transpose) {
    tmp_size = num_cols;
  } else {
    tmp_size = num_rows;
  }
  struct PseudobulkStatsSparse res;
  struct PseudobulkStatsTriplet trip;
  struct PseudobulkStatsTemp tmp;
  if ((method & PseudobulkStatsMethod::NonZeros) == PseudobulkStatsMethod::NonZeros) {
    res.non_zeros = Eigen::SparseMatrix<double>(num_rows, num_cols);
    trip.non_zeros.reserve(num_rows * num_cols);
    tmp.non_zeros = std::vector<double>(tmp_size, 0.0);
  }
  if ((method & PseudobulkStatsMethod::Sum) == PseudobulkStatsMethod::Sum) {
    res.sum = Eigen::SparseMatrix<double>(num_rows, num_cols);
    trip.sum.reserve(num_rows * num_cols);
    tmp.sum = std::vector<double>(tmp_size, 0.0);
  }
  if ((method & PseudobulkStatsMethod::Mean) == PseudobulkStatsMethod::Mean) {
    res.mean = Eigen::SparseMatrix<double>(num_rows, num_cols);
    trip.mean.reserve(num_rows * num_cols);
    tmp.non_zeros = std::vector<double>(tmp_size, 0.0);
    tmp.mean = std::vector<double>(tmp_size, 0.0);
  }
  if ((method & PseudobulkStatsMethod::Variance)== PseudobulkStatsMethod::Variance) {
    res.var = Eigen::SparseMatrix<double>(num_rows, num_cols);
    trip.var.reserve(num_rows * num_cols);
    tmp.non_zeros = std::vector<double>(tmp_size, 0.0);
    tmp.mean = std::vector<double>(tmp_size, 0.0);
    tmp.var = std::vector<double>(tmp_size, 0.0);
  }
  uint32_t ncells = 0;
  uint32_t g_idx = 0;
  // Use a one pass strategy that can do non-zero, mean, and variance at the same time
  while (it.nextCol()) {
    const uint32_t col = it.currentCol();
    if (user_interrupt != NULL && *user_interrupt) return res;
    while (it.load()) {
      const uint32_t *row_data = it.rowData();
      const T *val_data = it.valData();
      const uint32_t count = it.capacity();
      // Do transpose logic to calculate row_idx, col_idx out of loop, to avoid repeating conditional checks
      if (transpose) {
        // Doing variance calculation also does mean and non-zero count, so specific order is required to make sure we don't double count
        if ((method & PseudobulkStatsMethod::Variance) == PseudobulkStatsMethod::Variance) {
          for (uint32_t i = 0; i < count; i++) update_mean_and_variance(val_data[i], tmp.non_zeros[cell_groups[row_data[i]]], tmp.mean[cell_groups[row_data[i]]], tmp.var[cell_groups[row_data[i]]]);
        } else if ((method & PseudobulkStatsMethod::Mean) == PseudobulkStatsMethod::Mean) {
          for (uint32_t i = 0; i < count; i++) update_mean(val_data[i], tmp.non_zeros[cell_groups[row_data[i]]], tmp.mean[cell_groups[row_data[i]]]);
        } else if ((method & PseudobulkStatsMethod::NonZeros) == PseudobulkStatsMethod::NonZeros) {
          for (uint32_t i = 0; i < count; i++) tmp.non_zeros[cell_groups[row_data[i]]] += 1;
        }
        if ((method & PseudobulkStatsMethod::Sum) == PseudobulkStatsMethod::Sum) {
          for (uint32_t i = 0; i < count; i++) tmp.sum[cell_groups[row_data[i]]] += val_data[i];
        }
      } else {
        if ((method & PseudobulkStatsMethod::Variance) == PseudobulkStatsMethod::Variance) {
          for (uint32_t i = 0; i < count; i++) update_mean_and_variance(val_data[i], tmp.non_zeros[row_data[i]], tmp.mean[row_data[i]], tmp.var[row_data[i]]);
        } else if ((method & PseudobulkStatsMethod::Mean) == PseudobulkStatsMethod::Mean) {
          for (uint32_t i = 0; i < count; i++) update_mean(val_data[i], tmp.non_zeros[row_data[i]], tmp.mean[row_data[i]]);
        } else if ((method & PseudobulkStatsMethod::NonZeros) == PseudobulkStatsMethod::NonZeros) {
          for (uint32_t i = 0; i < count; i++) tmp.non_zeros[row_data[i]] += 1;
        }
        if ((method & PseudobulkStatsMethod::Sum) == PseudobulkStatsMethod::Sum) {
          for (uint32_t i = 0; i < count; i++) tmp.sum[row_data[i]] += val_data[i];
        }
      }
    }
    if (transpose) {
      for (uint32_t i = 0; i < tmp.var.size(); i++) {
        if (tmp.mean[i] > 0) trip.var.emplace_back(col, i, (tmp.var[i] + (tmp.mean[i] * tmp.mean[i] * tmp.non_zeros[i] * ((group_count_vec[i] - tmp.non_zeros[i]) / group_count_vec[i]))) / (group_count_vec[i] - 1));
      }
      for (uint32_t i = 0; i < tmp.mean.size(); i++) {
        if (tmp.mean[i] > 0) trip.mean.emplace_back(col, i, tmp.mean[i] * tmp.non_zeros[i] / group_count_vec[i]);
      }
      for (uint32_t i = 0; i < tmp.non_zeros.size(); i++) {
        if (tmp.non_zeros[i] > 0) trip.non_zeros.emplace_back(col, i, tmp.non_zeros[i]);
      }
      for (uint32_t i = 0; i < tmp.sum.size(); i++) {
        if (tmp.sum[i] > 0) trip.sum.emplace_back(col, i, tmp.sum[i]);
      }
      std::fill(tmp.var.begin(), tmp.var.end(), 0.0);
      std::fill(tmp.mean.begin(), tmp.mean.end(), 0.0);
      std::fill(tmp.sum.begin(), tmp.sum.end(), 0.0);
      std::fill(tmp.non_zeros.begin(), tmp.non_zeros.end(), 0.0);
      continue;
    }
    // transpose == false
    ncells += 1;
    if (ncells < group_count_vec[g_idx]) continue;
    for (uint32_t i = 0; i < tmp.var.size(); i++) {
      if (tmp.mean[i] > 0) trip.var.emplace_back(i, g_idx, (tmp.var[i] + (tmp.mean[i] * tmp.mean[i] * tmp.non_zeros[i] * (group_count_vec[g_idx] - tmp.non_zeros[i]) / group_count_vec[g_idx])) / (group_count_vec[g_idx] - 1));
    }
    for (uint32_t i = 0; i < tmp.mean.size(); i++) {
      if (tmp.mean[i] > 0) trip.mean.emplace_back(i, g_idx, tmp.mean[i] * tmp.non_zeros[i] / group_count_vec[g_idx]);
    }
    for (uint32_t i = 0; i < tmp.non_zeros.size(); i++) {
      if (tmp.non_zeros[i] > 0) trip.non_zeros.emplace_back(i, g_idx, tmp.non_zeros[i]);
    }
    for (uint32_t i = 0; i < tmp.sum.size(); i++) {
      if (tmp.sum[i] > 0) trip.sum.emplace_back(i, g_idx, tmp.sum[i]);
    }
    std::fill(tmp.var.begin(), tmp.var.end(), 0.0);
    std::fill(tmp.mean.begin(), tmp.mean.end(), 0.0);
    std::fill(tmp.sum.begin(), tmp.sum.end(), 0.0);
    std::fill(tmp.non_zeros.begin(), tmp.non_zeros.end(), 0.0);
    g_idx += 1;
    ncells = 0;
  }
  if ((method & PseudobulkStatsMethod::NonZeros) == PseudobulkStatsMethod::NonZeros) res.non_zeros.setFromTriplets(trip.non_zeros.begin(), trip.non_zeros.end());
  if ((method & PseudobulkStatsMethod::Sum) == PseudobulkStatsMethod::Sum) res.sum.setFromTriplets(trip.sum.begin(), trip.sum.end());
  if ((method & PseudobulkStatsMethod::Mean) == PseudobulkStatsMethod::Mean) res.mean.setFromTriplets(trip.mean.begin(), trip.mean.end());
  if ((method & PseudobulkStatsMethod::Variance) == PseudobulkStatsMethod::Variance) res.var.setFromTriplets(trip.var.begin(), trip.var.end());
  return res;
};
template PseudobulkStatsSparse pseudobulk_matrix_sparse<uint32_t>(std::unique_ptr<MatrixLoader<uint32_t>>&& mat,
                                                                  const std::vector<uint32_t>& cell_groups,
                                                                  PseudobulkStatsMethod method,
                                                                  bool transpose,
                                                                  std::atomic<bool> *user_interrupt);
template PseudobulkStatsSparse pseudobulk_matrix_sparse<float>(std::unique_ptr<MatrixLoader<float>>&& mat,
                                                               const std::vector<uint32_t>& cell_groups,
                                                               PseudobulkStatsMethod method,
                                                               bool transpose,
                                                               std::atomic<bool> *user_interrupt);
template PseudobulkStatsSparse pseudobulk_matrix_sparse<double>(std::unique_ptr<MatrixLoader<double>>&& mat,
                                                                const std::vector<uint32_t>& cell_groups,
                                                                PseudobulkStatsMethod method,
                                                                bool transpose,
                                                                std::atomic<bool> *user_interrupt);
} // namespace BPCells
