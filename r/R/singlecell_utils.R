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
#' Given a `(features x cells)` matrix, group cells by `cell_groups` and aggregate counts by `method` for each
#' feature.
#' @param cell_groups (Character/factor vector)  of cell groups/clusters. Length must be equal to the number of cells.
#' @param method (Character vector) Method(s) to aggregate counts. If one method is provided, the output will be a matrix. If multiple methods are provided, the output will be a named list of matrices.
#' 
#' Current options are:
#'  - `nonzeros`
#'  - `sum`
#'  - `mean`
#'  - `var`
#' @param clip_values (integer) If set, clip matrix values to the `<clip_values>`th quantile of each cell, between 0 and 1 inclusive.
#' @param threads (integer) Number of threads to use.
#' @return 
#'  - If `method` is length `1`, returns a matrix of shape `(features x groups)`.
#'  - If `method` is greater than length `1`, returns A list of matrices with each matrix representing a pseudobulk matrix with a different aggregation method.
#' Each matrix is of shape `(features x groups)`, and names are one of `nonzeros`, `sum`, `mean`, `var`.
#' @details The stats are ordered by complexity: nonzero, sum, mean, then variance. All
#' less complex stats are calculated in the process of calculating a more complicated stat.
#' That means it will not take additional time to calculate `mean` if `var` is also requested.
#' @inheritParams marker_features
#' @export
pseudobulk_matrix <- function(mat, cell_groups, method = "sum", clip_values = 1, threads = 1L) {
  assert_is(mat, "IterableMatrix")
  assert_is(cell_groups, c("factor", "character", "numeric"))
  cell_groups <- as.factor(cell_groups)
  assert_is(method, "character")
  methods <- c("var", "mean", "sum", "nonzeros")
  method <- match.arg(method, choices = methods, several.ok = TRUE)
  assert_is(clip_values, "numeric")
  assert_is(threads, "integer")
  # if multiple methods are provided, only need to pass in the top method as it will also calculate the less complex stats
  top_method <- intersect(methods, method)[1]
  top_method_number <- match(top_method, rev(methods)) - 1
  if (clip_values != 1) {
    quantile_values <- colQuantiles(mat, clip_values)
    # if all values are 0, skip
    mat <- min_by_col(mat, quantile_values)
  }
  iter <- iterate_matrix(parallel_split(mat, threads, threads*4))
  cell_groups <- as.factor(cell_groups)
  res <- pseudobulk_matrix_cpp(iter, cell_groups = as.integer(cell_groups), method = top_method_number, transpose = mat@transpose)
  # if res is a single matrix, return with colnames and rownames
  if (length(method) == 1) {
    colnames(res[[method]]) <- levels(cell_groups)
    rownames(res[[method]]) <- rownames(mat)
    return(res[[method]])
  }
  # give colnames and rownames for each matrix in res, which is a named list
  for (res_slot in names(res)) {
    if ((length(res[[res_slot]]) == 0) || !(res_slot %in% method)) {
      res[[res_slot]] <- NULL
    } else {
      colnames(res[[res_slot]]) <- levels(cell_groups)
      rownames(res[[res_slot]]) <- rownames(mat)
    }
  }
  return(res)
}


#' Find the nth quantile value of each column in a matrix.
#' @param x IterableMatrix object or a matrix-like object.
#' @param probs (Numeric) Quantile value(s) to be computed, between 0 and 1.
#' @return Numeric vector of quantile values for each column.
#' @describeIn IterableMatrix-methods Calculate colQuantiles (replacement for `matrixStats::colQuantiles`)
#' @export
colQuantiles <- function(x, probs = 0.99, ...) UseMethod("colQuantiles")
#' @export
colQuantiles.default <- function(x, probs = 0.99, ..., na.rm = FALSE, type = 7) {
  if (requireNamespace("MatrixGenerics", quietly = TRUE)) {
    MatrixGenerics::colQuantiles(x, probs = probs, na.rm = na.rm, type = type, ...)
  } else if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::colQuantiles(x, probs = probs, na.rm = na.rm, type = type, ...)
  } else {
    stop("Cannot run colQuantiles on a non-BPCells object unless MatrixGenerics or matrixStats is installed.")
  }
}
#' @export
colQuantiles.IterableMatrix <- function(x, probs = 0.99, ...) {
  assert_is(x, "IterableMatrix")
  assert_is(probs, "numeric")
  assert_true(all(probs >= 0 & probs <= 1))
  if (x@transpose) {
    stop("colQuantiles() does not support transposed matrices.\nPlease use transpose_storage_order().")
  }
  iter <- iterate_matrix(convert_matrix_type(x, "double"))
  res <- matrix_quantile_per_col_cpp(iter, probs)
  names(res) <- colnames(x)
  return(res)
}
rlang::on_load({
  if (requireNamespace("MatrixGenerics", quietly = TRUE)) {
    setMethod("colQuantiles", "IterableMatrix", colQuantiles.IterableMatrix)
  }
})