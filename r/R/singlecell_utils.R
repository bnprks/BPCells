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
#' @param z_score_norm (logical) If TRUE, z-score normalize the matrix before PCA.
#' @param n_dimensions (integer) Number of dimensions to keep during PCA.
#' @param scale_factor (integer) Scale factor for the tf-idf log transform.
#' @param threads (integer) Number of threads to use.
#' @param save_lsi (logical) If TRUE, save the SVD attributes for the matrix, as well as the idf normalization vector.
#' @return 
#' - If save_lsi is FALSE, return a dgCMatrix of shape (n_dimensions, ncol(mat)).
#' - If save_lsi is TRUE, return a list with the following elements:
#'  - **pca_res**: dgCMatrix of shape (n_dimensions, ncol(mat))
#' - **svd_attr**: List of SVD attributes
#' - **idf**: Inverse document frequency vector
#' @details Compute LSI through first doing a log(tf-idf) transform, z-score normalization, then PCA.  Tf-idf implementation is from Stuart & Butler et al. 2019.
#' 
#' Running on a 2600 cell dataset with 50000 peaks and 4 threads, as an example:
#' - 17.1 MB memory usage, 25.1 seconds runtime
#' @export
lsi <- function(
  mat,
  z_score_norm = TRUE, n_dimensions = 50L, scale_factor = 1e4,
  save_lsi = FALSE,
  threads = 1L
) {
  assert_is(mat, "IterableMatrix")
  assert_is_wholenumber(n_dimensions)
  assert_len(n_dimensions, 1)
  assert_greater_than_zero(n_dimensions)
  assert_true(n_dimensions < min(ncol(mat), nrow(mat)))
  assert_is_wholenumber(threads)

  # log(tf-idf) transform
  mat_stats <- matrix_stats(mat, row_stats = c("mean"), col_stats = c("mean"))
  
  npeaks <- colSums(mat) # Finding that sums are non-multithreaded and there's no interface to pass it in, but there is implementation in `ConcatenateMatrix.h`
  tf <- mat %>% multiply_cols(1 / npeaks)
  idf_ <- ncol(mat) / rowSums(mat)
  mat_tfidf <- tf %>% multiply_rows(idf_)
  mat_log_tfidf <- log1p(scale_factor * mat_tfidf)
  # Save to prevent re-calculation of queued operations
  mat_log_tfidf <- write_matrix_dir(mat_log_tfidf, tempfile("mat_log_tfidf"), compress = FALSE)
  # Z-score normalization
  if (z_score_norm) {
    cell_peak_stats <- matrix_stats(mat_log_tfidf, col_stats = "variance", threads = threads)$col_stats
    cell_means <- cell_peak_stats["mean",]
    cell_vars <- cell_peak_stats["variance",]
    mat_log_tfidf <- mat_log_tfidf %>%
      add_cols(-cell_means) %>%
      multiply_cols(1 / cell_vars)
  }
  # Run pca
  svd_attr_ <- svds(mat_log_tfidf, k = n_dimensions, threads = threads)
  pca_res <- t(svd_attr_$u) %*% mat_log_tfidf
  if(save_lsi) {
    return(list(
      pca_res = pca_res,
      svd_attr = svd_attr_,
      idf = idf_
    ))
  }
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
  mat_stats <- matrix_stats(mat, row_stats = c("variance"), threads = threads)
  feature_means <- mat_stats$row_stats["mean", ]
  feature_vars <- mat_stats$row_stats["variance", ]
  feature_means[feature_means == 0] <- 1e-12
  feature_dispersion <- feature_vars / feature_means
  feature_dispersion[feature_dispersion == 0] <- NA
  feature_dispersion <- log(feature_dispersion)
  feature_means <- log1p(feature_means)
  features_df <- data.frame(
    name = names(feature_means),
    vars = feature_vars, 
    means = feature_means,
    dispersion = feature_dispersion
  ) 
  features_df <- features_df %>% 
    dplyr::mutate(bin = cut(means, n_bins, labels=FALSE)) %>% 
    dplyr::group_by(bin) %>% 
    dplyr::mutate( 
      bin_mean = mean(dispersion, na.rm = TRUE), 
      bin_sd = sd(dispersion, na.rm = TRUE),
      bin_sd_is_na = is.na(bin_sd), 
      bin_sd = ifelse(bin_sd_is_na, bin_mean, bin_sd), # Set feats that are in bins with only one feat to have a norm dispersion of 1
      bin_mean = ifelse(bin_sd_is_na, 0, bin_mean),
      feature_dispersion_norm = (dispersion - bin_mean) / bin_sd
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(name, feature_dispersion_norm) %>%
    dplyr::arrange(desc(feature_dispersion_norm)) %>% 
    dplyr::slice(1:min(num_feats, nrow(.)))
  return(mat[features_df$name,])
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