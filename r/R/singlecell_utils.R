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
#' feature.  If the matrix is instead transposed, aggregated counts by `method` for each cell group.
#' @param method (string) Method to aggregate counts.  Current options are:
#'  - `sum` 
#'  - `mean`
#' @param clip_values (logical) If TRUE, clip high values to the 99th percentile of the data.
#' @return (data.frame) A data frame with row corresponding to the cell groups or features, and columns
#' coresponding to the aggregated counts, with each column named by the cell group.
#' @inheritParams marker_features
pseudobulk_counts <- function(mat, cell_groups, method = c("sum", "mean"), clip_values = FALSE) {
  assert_is(mat, "IterableMatrix")
  assert_is(cell_groups, c("factor", "character", "numeric"))
  cell_groups <- as.factor(cell_groups)
  method <- match.arg(method)
  assert_is(clip_values, "logical")
  iter <- iterate_matrix(convert_matrix_type(mat, "double"))
  if (clip_values) {
    if (mat@transpose) {
      percentile_values <- find_matrix_percentile_per_row_cpp(iter, 0.99)
      mat <- min_by_row(mat, percentile_values)
    } else {
      percentile_values <- find_matrix_percentile_per_col_cpp(iter, 0.99)
      mat <- min_by_col(mat, percentile_values)
    }
  }
  res <- pseudobulk_counts_cpp(iter, as.integer(cell_groups), method, clip_values, mat@transpose)
  return(res)
}


#' Aggregate counts matrices by cell group or feature.
#' 
#' Given a features x cells matrix, group cells by `cell_groups` and aggregate counts by `method` for each 
#' feature.
#' @param method (string) Method to aggregate counts.  Current options are:
#'  - `sum` 
#'  - `mean`
#' @param clip_values (logical) If TRUE, clip high values to the 99th percentile of the data.
#' @param threads (integer) Number of threads to use.
#' @return (data.frame) A data frame with row corresponding to features, and columns
#' coresponding to the aggregated counts, with each column named by the cell group.
#' @inheritParams marker_features
pseudobulk_counts_matrix_multiply <- function(mat, cell_groups, method = c("sum", "mean"), clip_values = FALSE, threads = 1L) {
  assert_is(mat, "IterableMatrix")
  assert_is(cell_groups, c("factor", "character", "numeric"))
  cell_groups <- as.factor(cell_groups)
  method <- match.arg(method)
  assert_is(clip_values, "logical")
  assert_is(threads, "integer")
  iter <- iterate_matrix(parallel_split(mat, threads, threads*4))
  if (clip_values) {
    if (mat@transpose) {
      percentile_values <- find_matrix_percentile_per_row_cpp(iter, 0.99)
      mat <- min_by_row(mat, percentile_values)
    } else {
      percentile_values <- find_matrix_percentile_per_col_cpp(iter, 0.99)
      mat <- min_by_col(mat, percentile_values)
    }
  }
  groups_membership_matrix <- cluster_membership_matrix(cell_groups, levels(cell_groups))
  if (method == "mean") {
    groups_membership_matrix <- multiply_cols(groups_membership_matrix, 1 / as.numeric(table(cell_groups)))
  }
  res <- as.matrix(as((mat %*% groups_membership_matrix), "dgCMatrix"))
  res <- as.data.frame(res)
  return(res)
}

#' Find the nth percentile value of each cell in a matrix.
matrix_percentile_per_cell <- function(mat, percentile = 0.99, threads = 1L) {
  assert_is(mat, "IterableMatrix")
  assert_is(percentile, "numeric")
  assert_true(percentile >= 0 && percentile < 1)
  iter <- iterate_matrix(convert_matrix_type(mat, "double"))
  if (mat@transpose) {
    res <- matrix_percentile_per_row_cpp(iter, percentile)
  } else {
    res <- matrix_percentile_per_col_cpp(iter, percentile)
    
  }
  names(res) <- colnames(mat)
  return(res)
}