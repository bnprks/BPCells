# Copyright 2023 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

#' Test for marker features
#'
#' Given a features x cells matrix, perform one-vs-all differential
#' tests to find markers.
#' 
#' Tips for using the values from this function:  
#' - Use `dplyr::mutate()` to add columns for e.g. adjusted p-value and log fold change.
#' - Use `dplyr::filter()` to get only differential genes above some given threshold
#' - To get adjusted p-values, use R `p.adjust()`, recommended method is "BH"
#' - To get log2 fold change: if your input matrix was already log-transformed,
#'   calculate `(foreground_mean - background_mean)/log(2)`. If your input
#'   matrix was not log-transformed, calculate `log2(forground_mean/background_mean)`
#'
#' @param mat IterableMatrix object of dimensions features x cells
#' @param groups Character/factor vector of cell groups/clusters. Length #cells
#' @param method Test method to use. Current options are:  
#'   - `wilcoxon`: Wilconxon rank-sum test a.k.a Mann-Whitney U test
#' @return tibble with the following columns:  
#'  - **foreground**: Group ID used for the foreground
#'  - **background**: Group ID used for the background (or NA if comparing to rest of cells)
#'  - **feature**: ID of the feature
#'  - **p_val_raw**: Unadjusted p-value for differential test
#'  - **foreground_mean**: Average value in the foreground group
#'  - **background_mean**: Average value in the background group
#' @export
marker_features <- function(mat, groups, method="wilcoxon") {
    assert_is(mat, "IterableMatrix")
    assert_true(length(groups) == ncol(mat))
    method <- match.arg(method)
    groups <- as.factor(groups)

    if (storage_order(mat) != "row") {
        rlang::inform(c(
            "Warning: marker features calculation requires row-major storage",
            "Consider using transpose_storage_order() if running marker_features repeatedly"
        ), .frequency = "regularly", .frequency_id = "marker_features_transpose")
        outdir <- tempfile("transpose")
        rlang::inform(sprintf("Writing transposed storage order to %s", outdir))
        mat <- transpose_storage_order(mat, outdir=outdir)
    }

    test_fn <- get(sprintf("wilcoxon_rank_sum_pval_%s_cpp", matrix_type(mat)))
    p_vals <- test_fn(iterate_matrix(mat), as.integer(groups) - 1)

    group_membership <- cluster_membership_matrix(groups, levels(groups))
    group_membership <- t(as(t(group_membership), "IterableMatrix"))

    group_sums <- as.matrix(as(t(mat %*% group_membership), "dgCMatrix")) # dim groups x features
    foreground_means <- multiply_rows(group_sums, 1 / as.numeric(table(groups)))
    background_means <- add_cols(-group_sums, colSums(group_sums)) %>%
      multiply_rows(1 / (length(groups) - as.numeric(table(groups))))
    
    feature_names <- if (is.null(rownames(mat))) seq_len(nrow(mat)) else rownames(mat)

    tibble::tibble(
        foreground = rep.int(levels(groups), nrow(mat)),
        background = NA_character_,
        feature = rep(feature_names, each=length(levels(groups))),
        p_val_raw = as.numeric(p_vals),
        foreground_mean = as.numeric(foreground_means),
        background_mean = as.numeric(background_means)
    )
}

#' Aggregate counts matrices by cell group or feature.
#' 
#' Given a features x cells matrix, group cells by `cell_groups` and aggregate counts by `method` for each
#' feature.
#' @param method (string) Method to aggregate counts.  Current options are:
#'  - `non-zeros`
#'  - `sum`
#'  - `mean`
#'  - `var`
#' @param clip_values (logical) If TRUE, clip high values to the 99th percentile of each cell.
#' @param threads (integer) Number of threads to use.
#' @return (Named list) A list of matrices with each matrix representing a pseudobulk matrix with a different aggregation method.
#' Each matrix is of shape (features x groups), and names are one of `non_zeros`, `sum`, `mean`, `var`.
#' @details The stats are ordered by complexity: nonzero, sum, mean, then variance. All
#' less complex stats are calculated in the process of calculating a more complicated stat.
#' So to calculate mean and variance simultaneously, just ask for variance,
#' which will compute mean and nonzero counts as a side-effect.
#' @inheritParams marker_features
#' @export
pseudobulk_matrix <- function(mat, cell_groups, method = c("sum", "mean", "var", "non-zeros"), clip_values = FALSE, threads = 1L) {
  assert_is(mat, "IterableMatrix")
  assert_is(cell_groups, c("factor", "character", "numeric"))
  cell_groups <- as.factor(cell_groups)
  method <- match.arg(method)
  assert_is(clip_values, "logical")
  assert_is(threads, "integer")

  if (clip_values) {
    quantile_values <- matrix_quantile_per_cell(mat, 0.99)
    mat <- min_by_col(mat, quantile_values)
  }
  iter <- iterate_matrix(parallel_split(mat, threads, threads*4))
  cell_groups <- as.factor(cell_groups)
  res <- pseudobulk_matrix_cpp(iter, cell_groups = as.integer(cell_groups), method = method, transpose = mat@transpose)
  # give colnames and rownames for each matrix in res, which is a named list
  for (res_slot in names(res)) {
    if (length(res[[res_slot]]) == 0) {
      res[[res_slot]] <- NULL
    } else {
      colnames(res[[res_slot]]) <- levels(cell_groups)
      rownames(res[[res_slot]]) <- rownames(mat)
    }
  }
  return(res)
}

#' Find the nth quantile value of each cell in a matrix.
#' @param quantile (numeric) quantile value to be found from each cell, between [0, 1].
#' @return (Numeric) cell quantile values with number of entries equal to number of colums in the matrix.
#' @inheritParams marker_features
#' @export
matrix_quantile_per_cell <- function(mat, quantile = 0.99) {
  assert_is(mat, "IterableMatrix")
  assert_is(quantile, "numeric")
  assert_true(quantile >= 0 && quantile < 1)
  
  # NOTE: Have to keep the entire matrix in memory if doing `quantile_per_row()`
  # for this operation which obviously doesn't mean it doesn't work if transposed.
  if (mat@transpose) stop("matrix_quantile_per_cell() does not support transposed matrices. \nPlease use transpose_storage_order().")
  iter <- iterate_matrix(convert_matrix_type(mat, "double"))
  res <- matrix_quantile_per_col_cpp(iter, quantile)
  names(res) <- colnames(mat)
  return(res)
}