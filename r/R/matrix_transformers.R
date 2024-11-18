# Copyright 2024 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.


#' Perform latent semantic indexing (LSI) on a matrix.
#' TODO: Add more details from when upstream LSI PR is reviewed.
#' @name LSITransformer
#' @export
setClass("LSITransformer",
  contains = "Transformer",
  slots = list(
    idf_ = "numeric",
    svd_attr_ = "list",
    z_score_norm = "logical",
    n_dimensions = "numeric",
    scale_factor = "numeric",
    threads = "integer"
  )
)

#' Create a new LSITransformer object
#' @export
LSITransformer <- function(z_score_norm = TRUE, n_dimensions = 50L, scale_factor = 1e4L, threads = 1L) {
  return(new(
    "LSITransformer", z_score_norm = z_score_norm, n_dimensions = n_dimensions, 
    scale_factor = scale_factor, threads = threads, step_name = "LSITransformer"))
}

#' @noRd
#' @export
setMethod("fit", signature(object = "LSITransformer", x = "IterableMatrix"), function(object, x, ...) {
  ret <- lsi(
    x, z_score_norm = object@z_score_norm, n_dimensions = object@n_dimensions, 
    scale_factor = object@scale_factor, threads = object@threads,
    save_lsi = TRUE
  )
  object@idf_ <- ret$idf
  object@svd_attr_ <- ret$svd_attr
  object@fitted <- TRUE
  return(object)
})


#' @noRd
#' @export
setMethod("project", signature(object = "LSITransformer", x = "IterableMatrix"), function(object, x, ...) {
  # rudimentary implementation -- Works but is duplicate code.
  assert_true(object@fitted)
  # Wait until LSI PR has been reviewed
  npeaks <- colSums(x) # Finding that sums are non-multithreaded and there's no interface to pass it in, but there is implementation in `ConcatenateMatrix.h`
  tf <- x %>% multiply_cols(1 / npeaks)
  mat_tfidf <- tf %>% multiply_rows(object@idf_)
  mat_log_tfidf <- log1p(object@scale_factor * mat_tfidf)
  mat_log_tfidf <- write_matrix_dir(mat_log_tfidf, tempfile("mat_log_tfidf"), compress = FALSE)
  if (object@z_score_norm) {
    cell_peak_stats <- matrix_stats(mat_log_tfidf, col_stats = "variance", threads = object@threads)$col_stats
    cell_means <- cell_peak_stats["mean",]
    cell_vars <- cell_peak_stats["variance",]
    mat_log_tfidf <- mat_log_tfidf %>%
      add_cols(-cell_means) %>%
      multiply_cols(1 / cell_vars)
  }
  pca_res <- t(object@svd_attr_$u) %*% mat_log_tfidf
  return(pca_res)
})


setMethod("short_description", "LSITransformer", function(x) {
  return(sprintf("LSITransformer(z_score_norm=%s, n_dimensions=%d, scale_factor=%d, threads=%d)",
    x@z_score_norm, x@n_dimensions, x@scale_factor, x@threads))
})


#' Perform feature selection on a matrix using dispersion.
#' TODO: Add more details from when upstream lsi PR is reviewed.
#' @name VarFeatSelectorTransformer
#' @export
setClass("VarFeatSelectorTransformer",
  contains = "Transformer",
  slots = list(
    features_ = "character",
    num_feats = "numeric",
    n_bins = "numeric",
    threads = "integer"
  )
)

#' Get the most variable features within a matrix.
#' @param num_feats (integer) Number of features to return.  If the number is higher than the number of features in the matrix,
#' all features will be returned.
#' @param n_bins (integer) Number of bins for binning mean gene expression.  Normalizing dispersion is done with respect to each bin,
#' and if the number of features
#' within a bin is less than 2, the dispersion is set to 1.
#' @param threads (integer) Number of threads to use.
#' @details The formula for calculating the most variable features is from the Seurat package (Satjia et al. 2015).
#'
#' Calculate using the following process:
#'  1. Calculate the dispersion of each feature (variance / mean)
#'  2. Log normalize dispersion and mean
#'  3. Bin the features by their means, and normalize dispersion within each bin
#' @export
VarFeatSelectorTransformer <- function(num_feats, n_bins = 20L, threads = 1L) {
  return(new(
    "VarFeatSelectorTransformer",
    num_feats = num_feats,
    n_bins = n_bins,
    threads = threads,
    step_name = "VarFeatSelectorTransformer")
  )
}

#' @noRd
#' @export
setMethod("fit", signature(object = "VarFeatSelectorTransformer", x = "IterableMatrix"), function(object, x, ...) {
  # Not sure what we should do in the case that x does not have rownames
  assert_true(!is.null(rownames(x)))
  ret <- highly_variable_features(x, num_feats = object@num_feats, n_bins = object@n_bins, save_feat_selection = TRUE, threads = object@threads)
  object@features_ <- ret$feature_selection$name
  object@fitted <- TRUE
  return(object)
})

#' @noRd
#' @export
setMethod("project", signature(object = "VarFeatSelectorTransformer", x = "IterableMatrix"), function(object, x, ...) {
  assert_true(object@fitted)
  return(x[object@features_,])
})

#' @noRd
#' @export
setMethod("short_description", "VarFeatSelectorTransformer", function(x) {
  return(sprintf("VarFeatSelectorTransformer(num_feats=%d, n_bins=%d, threads=%d)",
    x@num_feats, x@n_bins, x@threads))
})
