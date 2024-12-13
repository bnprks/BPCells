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
#' 
#' Given a `(features x cells)` matrix, perform LSI to perform tf-idf, z-score normalization, and PCA to create a latent space representation of the matrix of shape `(n_dimensions, ncol(mat))`.
#' @param mat (IterableMatrix) dimensions features x cells
#' @param n_dimensions (integer) Number of dimensions to keep during PCA.
#' @param corr_cutoff (numeric) Numeric filter for the correlation of a PC to the sequencing depth.  If the PC has a correlation that is great or equal to
#' the corr_cutoff, it will be excluded from the final PCA matrix.
#' @param scale_factor (integer) Scale factor for the tf-idf log transform.
#' @param save_lsi (logical) If `TRUE`, save the SVD attributes for the matrix, as well as the idf normalization vector.
#' @param threads (integer) Number of threads to use.
#' @return 
#' - If `save_lsi` is `FALSE`, return a dgCMatrix of shape `(n_dimensions, ncol(mat))`.
#' - If `save_lsi` is `TRUE`, return a list with the following elements:
#'    - `pca_res`: dgCMatrix of shape `(n_dimensions, ncol(mat))`
#'    - `unused_pcs`: Integer vector of PCs that were filtered out due to high correlation with sequencing depth
#'    - `svd_attr`: List of SVD attributes
#'    - `idf`: Inverse document frequency vector
#' @details Compute LSI through first doing a log(tf-idf) transform, z-score normalization, then PCA.  Tf-idf implementation is from Stuart & Butler et al. 2019.
#' 
#' Running on a 2600 cell dataset with 50000 peaks and 4 threads, as an example:
#' - 17.1 MB memory usage, 25.1 seconds runtime
#' @export
lsi <- function(
  mat,
  n_dimensions = 50L, corr_cutoff = 1, scale_factor = 1e4,
  save_lsi = FALSE,
  threads = 1L
) {
  assert_is(mat, "IterableMatrix")
  assert_is_wholenumber(n_dimensions)
  assert_len(n_dimensions, 1)
  assert_greater_than_zero(n_dimensions)
  assert_true(n_dimensions < min(ncol(mat), nrow(mat)))
  assert_true((corr_cutoff >= 0) && (corr_cutoff <= 1))
  assert_is_wholenumber(threads)

  # log(tf-idf) transform
  mat_stats <- matrix_stats(mat, row_stats = c("mean"), col_stats = c("mean"))
  read_depth <- mat_stats$col_stats["mean",] * nrow(mat) 
  tf <- mat %>% multiply_cols(1 / read_depth)
  idf_ <- 1 / mat_stats$row_stats["mean",]
  mat_tfidf <- tf %>% multiply_rows(idf_)
  mat_log_tfidf <- log1p(scale_factor * mat_tfidf)
  # Save to prevent re-calculation of queued operations
  mat_log_tfidf <- write_matrix_dir(
    convert_matrix_type(mat_log_tfidf, type = "float"), 
    tempfile("mat_log_tfidf"), compress = TRUE
  )
  # Run pca
  svd_attr_ <- svds(mat_log_tfidf, k = n_dimensions, threads = threads)
  pca_res <- t(svd_attr_$u) %*% mat_log_tfidf
  
  # Filter out PCs that are highly correlated with sequencing depth
  pca_corrs <- abs(cor(read_depth, t(pca_res)))
  pca_feats_to_keep <- which(pca_corrs < corr_cutoff)
  if (length(pca_feats_to_keep) != n_dimensions) {
    # not sure if this is the route we want to take.  Should we just leave the PCs in,
    # and not use them in downstream analysis?
    pca_res <- t(svd_attr_$u[, pca_feats_to_keep]) %*% mat_log_tfidf
  }

  if(save_lsi) {
    return(list(
      pca_res = pca_res,
      svd_attr = svd_attr_,
      unused_pcs <- which(pca_corrs >= corr_cutoff),
      idf = idf_
    ))
  }
  return(pca_res)
}

#' Get the most variable features within a matrix.
#' @param num_feats (integer) Number of features to return.  If the number is higher than the number of features in the matrix, 
#' all features will be returned.
#' @param n_bins (integer) Number of bins for binning mean gene expression.  Normalizing dispersion is done with respect to each bin, 
#' and if the number of features
#' within a bin is less than 2, the dispersion is set to 1.
#' @param save_feat_selection (logical) If `TRUE`, save the dispersions, means, and the features selected.
#' @returns 
#' - If `save_feat_selection` is `FALSE`, return an IterableMatrix subset of the most variable features of shape `(num_variable_features, ncol(mat))`.
#' - If `save_feat_selection` is `TRUE`, return a list with the following elements:
#'    - `mat`: IterableMatrix subset of the most variable features of shape `(num_variable_features, ncol(mat))`.
#'    - `feature_selection`: Dataframe with columns `name`, `mean`, `dispersion`, `bin`, `feature_dispersion_norm`.
#' @inheritParams lsi
#' @details The formula for calculating the most variable features is from the Seurat package (Satjia et al. 2015).
#' 
#' Calculate using the following process:
#'  1. Calculate the dispersion of each feature (variance / mean)
#'  2. Log normalize dispersion and mean
#'  3. Bin the features by their means, and normalize dispersion within each bin
#' @export
highly_variable_features_by_bin_variance <- function(
  mat, num_feats = 25000, n_bins = 20,
  save_feat_selection = FALSE,
  threads = 1L
) {
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
  # Calculate row information for dispersion
  mat_stats <- matrix_stats(mat, row_stats = c("variance"), threads = threads)
  feature_means <- mat_stats$row_stats["mean", ]
  feature_vars <- mat_stats$row_stats["variance", ]
  # Calculate dispersion, and log normalize
  feature_dispersion <- feature_vars / feature_means
  feature_dispersion[feature_dispersion == 0] <- NA
  feature_dispersion <- log(feature_dispersion)
  feature_dispersion[feature_means == 0] <- 0
  feature_means <- log1p(feature_means)
  features_df <- data.frame(
    name = names(feature_means),
    var = feature_vars, 
    mean = feature_means,
    dispersion = feature_dispersion
  ) 

  # Bin by mean, and normalize dispersion with each bin
  features_df <- features_df %>% 
    dplyr::mutate(bin = cut(mean, n_bins, labels=FALSE)) %>% 
    dplyr::group_by(bin) %>% 
    dplyr::mutate( 
      feature_dispersion_norm = (dispersion - mean(dispersion)) / sd(dispersion),
      feature_dispersion_norm = if (dplyr::n() == 1) {1} else {feature_dispersion_norm} # Set feats that are in bins with only one feat to have a norm dispersion of 1  
    ) %>%  
    dplyr::ungroup() %>%
    dplyr::slice_max(order_by = feature_dispersion_norm, n = num_feats, with_ties = FALSE)
  feats_of_interest <- which(rownames(mat) %in% features_df$name) # get rownames to get original sorted order
  if (save_feat_selection) {
    # get rownames that are in features_df$name
    return(list(
      mat = mat[feats_of_interest,],
      feature_selection = features_df
    ))
  }
  return(mat[feats_of_interest,])
}

#' Get the most variable features within a matrix, across clusters
#' @param clusters (factor) Vector of cluster assignments for each cell. Length must be `ncol(mat)`.
#' @param num_feats (integer) Number of features to return.  If the number is higher than the number of features in the matrix, 
#' all features will be returned.
#' @param log_normalize (logical) If `TRUE`, log normalize the matrix pseudobulked matrix prior to variance calculation
#' @param scale_factor (numeric) Scale factor for log normalization.  Unused if `log_normalize` is `FALSE`.
#' @param return_subsetted_matrix (logical) If `TRUE`, return the subsetted matrix of the most variable features.
#' @returns
#' - If `return_subsetted_matrix` is `FALSE`, return a dataframe with the following columns, sorted descending by variance:
#'   - `names`: Feature name
#'   - `score`: Variance of the feature
#'   - `highly_variable`: Logical vector of whether the feature is highly variable
#' - Else, return the original matrix subsetted to the most variable features of shape `(num_variable_features, ncol(mat))`.
#' @inheritParams highly_variable_features
#' @details 
#' Calculate using the following process:
#'  1. Pseudobulk cells by cluster
#'  2. Normalize by term frequency for each pseudobulk
#'  3. Do an optional log normalization of the pseudobulk matrix
#'  4. Find `num_feats` features with the highest variance across clusters
highly_variable_features_by_cluster_variance <- function(
  mat, clusters, num_feats = 25000, 
  log_normalize = TRUE, scale_factor = 1e4, return_subsetted_matrix = FALSE,
  threads = 1L, verbose = TRUE
) {
  assert_is(mat, "IterableMatrix")
  assert_greater_than_zero(num_feats)
  assert_is_wholenumber(num_feats)
  assert_len(num_feats, 1)
  assert_true(num_feats <= nrow(mat))
  assert_is(clusters, c("factor", "character", "numeric"))
  assert_true(length(clusters) == ncol(mat))
  assert_true(is.logical(log_normalize))
  if(log_normalize) assert_greater_than_zero(scale_factor)
  assert_true(is.logical(return_subsetted_matrix))
  assert_true(is.logical(verbose))

  # Pseudobulk the matrix by cluster
  group_mat <- pseudobulk_matrix(mat, clusters, method = "sum", threads = threads)
  group_mat_normalized <- t(t(group_mat) / colSums(group_mat))
  
  if (log_normalize) group_mat_normalized <- log1p(group_mat_normalized * scale_factor)
  features_df <- tibble::tibble(
    names = rownames(mat),
    score = rowVars(group_mat_normalized)
  ) %>% # Give the first num_feats features a highly_variable of TRUE
    dplyr::arrange(desc(score)) %>% 
    dplyr::mutate(highly_variable = dplyr::row_number() <= num_feats)
  if (return_subsetted_matrix) {
    # Do not alter feature order from original
    filtered_df <- features_df %>% dplyr::filter(highly_variable)
    feats_of_interest <- which(rownames(mat) %in% filtered_df$name)
    return(mat[feats_of_interest,])
  }
  return(features_df)
}

#' Get the top features from a matrix, based on the sum of each feature.
#' @param num_feats Number of features to return.  If the number is higher than the number of features in the matrix,
#' all features will be returned.
#' @inheritParams highly_variable_features
#' @return Integer vector of indices of the top features in mat.
#' @export
top_features <- function(mat, num_feats, threads = 1L) {
  assert_is(mat, "IterableMatrix")
  assert_is_wholenumber(num_feats)
  assert_greater_than_zero(num_feats)
  assert_is_wholenumber(threads)
  assert_greater_than_zero(threads)
  assert_true(num_feats <= nrow(mat))

  # get the sum of each feature, binarized
  feature_sums <- matrix_stats(mat, row_stats = "nonzero", threads = threads)$row_stats["nonzero",]

  # get the top features
  top_features <- order(feature_sums, decreasing = TRUE)[1:num_feats]
  return(top_features)
}

#' Run iterative LSI on a matrix.
#' 
#' This function will compute an iterative LSI dimensionality reduction on an ArchRProject.
#' @param mat The name of the data matrix to retrieve from the ArrowFiles associated with the `ArchRProject`. Valid options are
#' "TileMatrix" or "PeakMatrix".
#' @param n_iterations The number of LSI iterations to perform.
#' @param first_feature_selection_method First iteration selection method for features to use for LSI. Either "Top" for the top accessible/average or "Var" for the top variable features. 
#' "Top" should be used for all scATAC-seq data (binary) while "Var" should be used for all scRNA/other-seq data types (non-binary).
#' @param first_feature_selection_args Arguments to pass to the first feature selection method. If using "Top", refer to `top_features()`. If using "Var", refer to `highly_variable_features()`.
#' @param var_feature_selection_args Arguments to pass to the variable feature selection method. Refer to `highly_variable_features()`.
#' @param lsi_args Arguments to pass to the LSI function. Refer to `lsi()`.
#' @param cluster_method Method to use for clustering. Current options are "leiden", "louvain", and "seurat".
#' @param clustering_args Arguments to pass to the clustering function. Refer to `cluster_graph_leiden()`, `cluster_graph_louvain()`, or `cluster_graph_seurat()`.
#' @return A simple list of LSI results and features with the following format:
#' - `pca_res`: dgCMatrix of shape `(n_dimensions, ncol(mat))`
#' 
#' @seealso `lsi()`, `top_features()`, `highly_variable_features()`
#' @inheritParams lsi
iterative_lsi <- function(
    mat, 
    n_iterations = 2,
    first_feature_selection_method = c("Top", "Var"),
    first_feature_selection_args = list(),
    var_feature_selection_args = list(),
    lsi_args = list(),
    cluster_method = c("leiden", "louvain", "seurat"),
    clustering_args = list(),
    saveIterations, threads = 1L
) {
  assert_is(mat, "IterableMatrix")
  assert_true(n_iterations > 0)
  assert_is_wholenumber(n_iterations)
  assert_is(firstSelection, "character")
  assert_is(threads, "integer")

  cluster_func <- get(sprintf("cluster_graph_%s", cluster_method))
  lsi_args <- c(lsi_args, mat = mat, threads = threads)
  clustering_args <- c(clustering_args, data = mat, threads = threads)
  var_feature_selection_args <- c(var_feature_selection_args, mat = mat, threads = threads)
  for (i in seq_len(iterations)) {
    # For the first iteration, we need to select the features used by LSI to create the first dim reduction
    if ((i == 1)) {
      
      features <- ifelse(
        tolower(firstSelection) == "top", 
        do.call(top_features, first_feature_selection_args),
        do.call(highly_variable_features, first_feature_selection_args)
      )
    }
    # Perform LSI
    pca_res <- do.call(lsi, lsi_args)
    if (n_iterations == 1) {
      return(pca_res)
    }
    # Every feature selection, if done iteratively, will be based on variable feature selection based off of clustering.
    # Cluster LSI
    # Placeholder here for now, just to check if everything works
    clusts <- mat %>% knn_hnsw(threads = threads) %>% do.call(cluster_func, clustering_args)
  }
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