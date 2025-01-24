# Copyright 2023 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.


#################
# Feature selection
#################

#' Feature selection functions
#' 
#' Apply a feature selection method to a `(features x cells)` matrix.
#' @rdname feature_selection
#' @param mat (IterableMatrix) dimensions features x cells
#' @param num_feats (float) Number of features to return.  If the number given is between 0 and 1, treat as a proportion of 
#' the number of rows, rounded down.  Otherwise, treat as an absolute number.
#' If the number is higher than the number of features in the matrix, 
#' all features will be returned.
#' @param normalize (function) Normalize matrix using a given function. If `NULL`, no normalization is performed.
#' @param threads (integer) Number of threads to use.
#' @returns
#' Return a dataframe with the following columns, sorted descending by score:
#' - `names`: Feature name.
#' - `score`: Scoring of the feature, depending on the method used.  
#' - `highly_variable`: Logical vector of whether the feature is highly variable.
#'
#' Each different feature selection method will have a different scoring method:
#' - `select_features_by_variance`: Score representing variance of each feature.
#' @details 
#' `select_features_by_variance` Calculates the variance of each feature using the following process:
#'  1. Perform an optional term frequency + log normalization, for each feature.
#'  2. Find `num_feats` features with the highest variance.
#' @export
select_features_variance <- function(
  mat, num_feats = 0.05, 
  normalize = NULL,
  threads = 1L
) {
  assert_greater_than_zero(num_feats)
  assert_is_wholenumber(num_feats)
  assert_len(num_feats, 1)
  assert_is(num_feats, "numeric")
  if (rlang::is_missing(mat)) {
    return(create_partial(
      missing_args = list(
        num_feats = missing(num_feats),
        normalize = missing(normalize), 
        threads = missing(threads)
      )
    ))
  }
  assert_is(mat, "IterableMatrix")
  num_feats <- min(max(num_feats, 0), nrow(mat))
  if (num_feats < 1 && num_feats > 0) num_feats <- floor(nrow(mat) * num_feats)
  if (!is.null(normalize)) mat <- partial_apply(normalize, threads = threads)(mat)
  features_df <- tibble::tibble(
    names = rownames(mat),
    score = matrix_stats(mat, row_stats = "variance", threads = threads)$row_stats["variance",]
  ) %>% 
    dplyr::arrange(desc(score)) %>% 
    dplyr::mutate(highly_variable = dplyr::row_number() <= num_feats)
  return(features_df)
}


#' @rdname feature_selection
#' @returns
#'  - `select_features_by_dispersion`: Score representing the dispersion of each feature.
#' @details 
#' `select_features_by_dispersion` calculates the dispersion of each feature using the following process:
#'  1. Perform an optional term frequency + log normalization, for each feature.
#'  2. Find the dispersion (variance/mean) of each feature.
#'  3. Find `num_feats` features with the highest dispersion.
select_features_dispersion <- function(
  mat, num_feats = 0.05, 
  normalize = NULL,
  threads = 1L
) {
  assert_greater_than_zero(num_feats)
  assert_is_wholenumber(num_feats)
  assert_len(num_feats, 1)
  assert_is(num_feats, "numeric")
  if (rlang::is_missing(mat)) {
    return(create_partial(
      missing_args = list(
        num_feats = missing(num_feats),
        normalize = missing(normalize), 
        threads = missing(threads)
      )
    ))
  }
  num_feats <- min(max(num_feats, 0), nrow(mat))
  if (num_feats < 1 && num_feats > 0) num_feats <- floor(nrow(mat) * num_feats)
  assert_is(mat, "IterableMatrix")
  if (!is.null(normalize)) mat <- partial_apply(normalize, threads = threads)(mat)
  mat_stats <- matrix_stats(mat, row_stats = "variance", threads = threads)
  features_df <- tibble::tibble(
    names = rownames(mat),
    score = mat_stats$row_stats["variance", ] / mat_stats$row_stats["mean", ]
  ) %>% 
    dplyr::arrange(desc(score)) %>% 
    dplyr::mutate(highly_variable = dplyr::row_number() <= num_feats)
  return(features_df)
}


#' @rdname feature_selection
#' @returns
#' - `select_features_by_mean`: Score representing the mean accessibility of each feature.
#' @details 
#' `select_features_by_mean` calculates the mean accessibility of each feature using the following process:
#' 1. Get the sum of each binarized feature.
#' 2. Find `num_feats` features with the highest accessibility.
#' @export
select_features_mean <- function(mat, num_feats = 0.05, normalize = NULL, threads = 1L) {
  assert_is_wholenumber(num_feats)
  assert_greater_than_zero(num_feats)
  assert_is(num_feats, "numeric")
  if (rlang::is_missing(mat)) {
    return(create_partial(
      missing_args = list(
        num_feats = missing(num_feats),
        normalize = missing(normalize), 
        threads = missing(threads)
      )
    ))
  }
  assert_is(mat, "IterableMatrix")
  num_feats <- min(max(num_feats, 0), nrow(mat))
  if (num_feats < 1 && num_feats > 0) num_feats <- floor(nrow(mat) * num_feats)
  if (!is.null(normalize)) mat <- partial_apply(normalize, threads = threads)(mat)
  # get the sum of each feature, binarized
  # get the top features
  
  features_df <- tibble::tibble(
    names = rownames(mat),
    score = matrix_stats(mat, row_stats = "nonzero", threads = threads)$row_stats["nonzero", ]
  ) %>%
    dplyr::arrange(desc(score)) %>%
    dplyr::mutate(highly_variable = dplyr::row_number() <= num_feats)
  return(features_df)
}


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
#' @param cell_groups (Character/factor) Vector of group/cluster assignments for each cell. Length must be `ncol(mat)`.
#' @param method (Character vector) Method(s) to aggregate counts. If one method is provided, the output will be a matrix. If multiple methods are provided, the output will be a named list of matrices.
#'
#' Current options are: `nonzeros`, `sum`, `mean`, `variance`.
#' @param threads (integer) Number of threads to use.
#' @return
#'  - If `method` is length `1`, returns a matrix of shape `(features x groups)`.
#'  - If `method` is greater than length `1`, returns a list of matrices with each matrix representing a pseudobulk matrix with a different aggregation method.
#' Each matrix is of shape `(features x groups)`, and names are one of `nonzeros`, `sum`, `mean`, `variance`.
#' @details Some simpler stats are calculated in the process of calculating more complex
#' statistics. So when calculating `variance`, `nonzeros` and `mean` can be included with no
#' extra calculation time, and when calculating `mean`, adding `nonzeros` will take no extra time.
#' @inheritParams marker_features
#' @export
pseudobulk_matrix <- function(mat, cell_groups, method = "sum", threads = 0L) {
  assert_is(mat, "IterableMatrix")
  assert_is(cell_groups, c("factor", "character", "numeric"))
  assert_true(length(cell_groups) == ncol(mat))
  cell_groups <- as.factor(cell_groups)
  assert_is(method, "character")
  methods <- c("variance", "mean", "sum", "nonzeros")
  for (m in method) {
    if (!(m %in% methods)) {
      rlang::abort(sprintf("method must be one of: %s", paste(methods, collapse = ", ")))
    }
  }
  assert_is(threads, "integer")

  it <- mat %>%
    convert_matrix_type("double") %>%
    parallel_split(threads, threads*4) %>%
    iterate_matrix()
  
  res <- pseudobulk_matrix_cpp(it, cell_groups = as.integer(cell_groups) - 1, method = method, transpose = mat@transpose)
  # if res is a single matrix, return with colnames and rownames
  if (length(method) == 1) {
    colnames(res[[method]]) <- levels(cell_groups)
    rownames(res[[method]]) <- rownames(mat)
    return(res[[method]])
  }
  # give colnames and rownames for each matrix in res, which is a named list
  for (res_slot in names(res)) {
    # Filter out methods that weren't requested
    if (!(res_slot %in% method)) {
      res[[res_slot]] <- NULL
    } else {
      colnames(res[[res_slot]]) <- levels(cell_groups)
      rownames(res[[res_slot]]) <- rownames(mat)
    }
  }
  return(res)
}