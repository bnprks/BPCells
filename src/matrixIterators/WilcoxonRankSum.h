#pragma once

#ifndef RCPP_EIGEN
#include <Eigen/SparseCore>
#else
#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <RcppEigen.h>
#endif

#include "ColwiseRank.h"
namespace BPCells {

// Perform Wilcoxon rank sum test (aka Mann-Whitney U test) on each column of a matrix.
// Given the groupings specified by `groups`, perform a 1-vs-rest test, and return a matrix
// of p-values with dimensions (# groups) x (# columns)
template<typename T>
Eigen::MatrixXd wilcoxon_rank_sum(std::unique_ptr<MatrixLoader<T>> &&mat, const std::vector<uint32_t> &groups, std::atomic<bool> *user_interrupt);

}
