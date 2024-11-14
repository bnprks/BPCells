# Copyright 2024 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

# #' Perform latent semantic indexing (LSI) on a matrix.
# #' @name LSITransformer
# #' @export
# setClass("LSITransformer",
#   contains = "Transformer",
#   slots = list(
#     idf_ = "numeric",
#     svd_attr_ = "list",
#     z_score_norm = "logical",
#     n_dimensions = "integer",
#     scale_factor = "integer",
#     threads = "integer"
#   ),
#   prototype = list(
#     idf_ = numeric(0),
#     svd_attr_ = list(),
#     z_score_norm = FALSE,
#     n_dimensions = 20L,
#     scale_factor = 1e4L,
#     threads = 1L
#   )
# )

# #' Create a new LSITransformer object
# #' @export
# LSITransformer <- function(z_score_norm, n_dimensions, scale_factor, threads) {
#   return(new(
#     "LSITransformer", z_score_norm = z_score_norm, n_dimensions = n_dimensions, 
#     scale_factor = scale_factor, threads = threads, step_name = "LSITransformer"))
# }

# setMethod("fit", signature(object = "LSITransformer", x = "IterableMatrix"), function(object, x, ...) {
#   ret <- lsi(
#     x, z_score_norm = object@z_score_norm, n_dimensions = object@n_dimensions, 
#     scale_factor = object@scale_factor, threads = object@threads,
#     save_lsi = TRUE
#   )
#   object@idf_ <- ret$idf
#   object@svd_attr_ <- ret$svd_attr
#   object@fitted <- TRUE
#   return(object)
# })

# setMethod("transform", signature(object = "LSITransformer", x = "IterableMatrix"), function(object, x, ...) {
#   # rudimentary implementation -- Works but is duplicate code.  
#   assert_true(object@fitted)
#   # Wait until LSI PR has been reviewed
#   npeaks <- colSums(x) # Finding that sums are non-multithreaded and there's no interface to pass it in, but there is implementation in `ConcatenateMatrix.h`
#   tf <- x %>% multiply_cols(1 / npeaks)
#   mat_tfidf <- tf %>% multiply_rows(object@idf_)
#   mat_log_tfidf <- log1p(object@scale_factor * mat_tfidf)
#   mat_log_tfidf <- write_matrix_dir(mat_log_tfidf, tempfile("mat_log_tfidf"), compress = FALSE)
#   if (object@z_score_norm) {
#     cell_peak_stats <- matrix_stats(mat_log_tfidf, col_stats = "variance", threads = object@threads)$col_stats
#     cell_means <- cell_peak_stats["mean",]
#     cell_vars <- cell_peak_stats["variance",]
#     mat_log_tfidf <- mat_log_tfidf %>%
#       add_cols(-cell_means) %>%
#       multiply_cols(1 / cell_vars)
#   }
#   pca_res <- t(object@svd_attr_$u) %*% mat_log_tfidf
#   return(pca_res)
# })

# setMethod("short_description", "LSITransformer", function(x) {
#   return(sprintf("LSITransformer(z_score_norm=%s, n_dimensions=%d, scale_factor=%d, threads=%d)",
#     x@z_score_norm, x@n_dimensions, x@scale_factor, x@threads))
# })

# setClass("VarFeatSelectorTransformer",
#   contains = "Transformer",
#   slots = list(
#     num_feats = "integer",
#     n_bins = "integer"
#   )
# )