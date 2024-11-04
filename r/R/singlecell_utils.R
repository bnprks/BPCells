# Copyright 2024 BPCells contributors
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


#' Perform latent semantic indexing (LSI) on a matrix.
#' @param mat (IterableMatrix) dimensions features x cells
#' @param n_dimensions (integer) Number of dimensions to keep during PCA.
#' @param scale_factor (integer) Scale factor for the tf-idf log transform.
#' @param save_in_memory (logical) If TRUE, save the log(tf-idf) matrix in memory.  
#' If FALSE, save to a temporary location in disk.  Saving in memory will result in faster downstream operations,
#' but will require in higher memory usage.  Comparison of memory usage and speed is in the details section.
#' @param threads (integer) Number of threads to use.
#' @return dgCMatrix of shape (n_dimensions, ncol(mat)).
#' @details Compute LSI through first doing a log(tf-idf) transform, z-score normalization, then PCA.  Tf-idf implementation is from Stuart & Butler et al. 2019.
#' 
#' ** Saving in memory vs disk: **
#' Following the log(tf-idf) transform, the matrix is stored into a temporary location, as the next step will break the sparsity pattern of the matrix.
#' This is done to prevent re-calculation of queued operations during PCA optimization.
#' 
#' Running on a 2600 cell dataset with 50000 peaks and 4 threads, as an example:
#' - Saving in memory: 233 MB memory usage, 22.7 seconds runtime
#' - Saving in disk: 17.1 MB memory usage, 25.1 seconds runtime
#' 
#' @export
lsi <- function(mat, n_dimensions = 50L, scale_factor = 1e4, save_in_memory = FALSE, threads = 1L) {
  assert_is(mat, "IterableMatrix")
  assert_is_wholenumber(n_dimensions)
  assert_len(n_dimensions, 1)
  assert_greater_than_zero(n_dimensions)
  assert_true(n_dimensions < min(ncol(mat), nrow(mat)))
  assert_is_wholenumber(threads)

  # log(tf-idf) transform
  npeaks <- colSums(mat) # Finding that sums are non-multithreaded and there's no interface to pass it in, but there is implementation in `ConcatenateMatrix.h`
  tf <- mat %>% multiply_cols(1 / npeaks)
  idf_ <- ncol(mat) / rowSums(mat)
  mat_tfidf <- tf %>% multiply_rows(idf_)
  mat_log_tfidf <- log1p(scale_factor * mat_tfidf)
  # Save to prevent re-calculation of queued operations
  if (save_in_memory) {
    mat_log_tfidf <- write_matrix_memory(mat_log_tfidf, compress = FALSE)
  } else {
    mat_log_tfidf <- write_matrix_dir(mat_log_tfidf, tempfile("mat_log_tfidf"), compress = FALSE)
  } 
  # Z-score normalization
  cell_peak_stats <- matrix_stats(mat_log_tfidf, col_stats="variance", threads = threads)$col_stats
  cell_means <- cell_peak_stats["mean",]
  cell_vars <- cell_peak_stats["variance",]
  mat_lsi_norm <- mat_log_tfidf %>%
    add_cols(-cell_means) %>%
    multiply_cols(1 / cell_vars)
  # Run pca
  svd_attr_ <- svds(mat_lsi_norm, k = n_dimensions, threads = threads)
  pca_res <- t(svd_attr_$u) %*% mat_lsi_norm
  return(pca_res)
}

#' Get the most variable features within a matrix
#' @param num_feats (integer) Number of features to return.  If the number is higher than the number of features in the matrix, 
#' ll features will be returned.
#' @param n_bins (integer) Number of bins for binning mean gene expression.  Normalizing dispersion is done with respect to each bin, 
#' and if the number of features
#' within a bin is less than 2, the dispersion is set to 1.
#' @returns IterableMatrix subset of the most variable features.
#' @inheritParams lsi
#' @details The formula for calculating the most variable features is from the Seurat package (Satjia et al. 2015).
#' 
#' Calculate using the following process:
#'  1. Calculate the dispersion of each feature (variance / mean)
#'  2. Log normalize dispersion and mean
#'  3. Bin the features by their means, and normalize dispersion within each bin
#' @export
highly_variable_features <- function(mat, num_feats, n_bins, threads = 1L) {
  assert_is(mat, "IterableMatrix")
  assert_greater_than_zero(num_feats)
  assert_is_wholenumber(num_feats)
  assert_len(num_feats, 1)
  assert_is_wholenumber(n_bins)
  assert_len(n_bins, 1)
  assert_greater_than_zero(n_bins)
  if (nrow(mat) <= num_feats) {
    log_progress(sprintf("Number of features (%s) is less than num_feats (%s), returning all features", nrow(mat), num_feats))
    return(mat)
  }
  
  feature_means <- matrix_stats(mat, row_stats = c("mean"))$row_stats["mean", ]
  feature_vars <- matrix_stats(mat, row_stats = c("variance"))$row_stats["variance", ]
  feature_means[feature_means == 0] <- 1e-12
  feature_dispersion <- feature_vars / feature_means
  feature_dispersion[feature_dispersion == 0] <- NA
  feature_dispersion <- log(feature_dispersion)
  feature_means <- log1p(feature_means)
  mean_bins <- cut(feature_means, n_bins, labels = FALSE)
  
  bin_mean <- tapply(feature_dispersion, mean_bins, function(x) mean(x, na.rm = TRUE))
  bin_sd <- tapply(feature_dispersion, mean_bins, function(x) sd(x, na.rm = TRUE))
  # Set feats that are in bins with only one feat to have a norm dispersion of 1
  one_gene_bin <- is.na(bin_sd)
  bin_sd[one_gene_bin] <- bin_mean[one_gene_bin]
  bin_mean[one_gene_bin] <- 0
  # map mean_bins indices to bin_stats
  # Do a character search as bins without features mess up numeric indexing
  feature_dispersion_norm <- (feature_dispersion - bin_mean[as.character(mean_bins)]) / bin_sd[as.character(mean_bins)]
  names(feature_dispersion_norm) <- names(feature_dispersion)
  feature_dispersion_norm <- sort(feature_dispersion_norm) # sorting automatically removes NA values
  if (length(feature_dispersion_norm) < num_feats) log_progress(sprintf("Number of features (%s) is less than num_feats (%s), returning all non-zero features", length(feature_dispersion_norm), num_feats))
  variable_features_ <- feature_dispersion_norm[max(1, (length(feature_dispersion_norm) - num_feats + 1)):length(feature_dispersion_norm)]
  return(mat[names(variable_features_), ])
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
pseudobulk_matrix <- function(mat, cell_groups, method = "sum", threads = 1L) {
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
  # if multiple methods are provided, only need to pass in the top method as it will also calculate the less complex stats
  iter <- iterate_matrix(parallel_split(mat, threads, threads*4))
  res <- pseudobulk_matrix_cpp(iter, cell_groups = as.integer(cell_groups) - 1, method = method, transpose = mat@transpose)
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