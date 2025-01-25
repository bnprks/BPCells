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

#' @rdname feature_selection
#' @param n_bins (integer) Number of bins for binning mean gene expression.  Normalizing dispersion is done with respect to each bin, 
#' and if the number of features
#' within a bin is less than 2, the dispersion is set to 1.
#' @returns
#'  - `select_features_binned_dispersion`: Score representing the bin normalized dispersion of each feature.
#' @details 
#' `select_features_binned_dispersion` calculates the bin normalized dispersion of each feature using the following process, given by the Seurat package (Satjia et al. 2015):
#'  1. Calculate the dispersion of each feature (variance / mean)
#'  2. Log normalize dispersion and mean
#'  3. Bin the features by their means, and normalize dispersion within each bin
#' @export
select_features_binned_dispersion <- function(
  mat, num_feats = 25000, n_bins = 20,
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
    features_df <- tibble::tibble(
      names = rownames(mat),
      score = rep(0, nrow(mat)),
      highly_variable = rep(TRUE, nrow(mat))
    )
    return(mat)
  }
  # Calculate row information for dispersion
  mat_stats <- matrix_stats(mat, row_stats = c("variance"), threads = threads)
  feature_means <- mat_stats$row_stats["mean", ]
  feature_vars <- mat_stats$row_stats["variance", ]
  # Calculate dispersion, and log normalize
  feature_dispersion <- feature_vars / feature_means
  feature_dispersion[feature_dispersion == 0] <- NA
  feature_dispersion <- log1p(feature_dispersion)
  feature_dispersion[feature_means == 0] <- 0
  feature_means <- log1p(feature_means)
  features_df <- tibble::tibble(
    names = names(feature_means),
    var = feature_vars, 
    mean = feature_means,
    dispersion = feature_dispersion
  )
  # Bin by mean, and normalize dispersion with each bin
  features_df <- features_df %>% 
    dplyr::mutate(bin = cut(mean, n_bins, labels=FALSE)) %>% 
    dplyr::group_by(bin) %>% 
    dplyr::mutate( 
      score = (dispersion - mean(dispersion)) / sd(dispersion),
      score = if (dplyr::n() == 1) {1} else {score} # Set feats that are in bins with only one feat to have a norm dispersion of 1  
    ) %>% 
    dplyr::ungroup() %>%
    dplyr::arrange(desc(score)) %>%
    dplyr::mutate(highly_variable = dplyr::row_number() <= num_feats) %>% 
    dplyr::select(c("names", "score", "highly_variable"))
  return(features_df)
}


#################
# DimReduction Class Definition
#################

#' Barebones definition of a DimReduction class.
#' 
#' Represents a latent space output of a matrix after a transformation function, with any required information to reproject other inputs using this object.
#' Child classes should implement a `project` method to allow for the projection of other matrices using
#' the fitted transformation object.
#' @field cell_embeddings (IterableMatrix, dgCMatrix, matrix) The projected data
#' @field fitted_params (list) A list of parameters used for the transformation of a matrix.
#' @export
DimReduction <- function(x, fitted_params = list(), ...) {
  assert_is(x, c("IterableMatrix", "dgCMatrix", "matrix"))
  assert_is(fitted_params, "list")
  structure(
    list(
      cell_embeddings = x,
      fitted_params = fitted_params,
      ...
    ),
    class = "DimReduction"
  )
}

#' Perform a dimensionality reduction on a matrix using a pre-fit DimReduction object.
#' @param x DimReduction object.
#' @param mat IterableMatrix object.
#' @return IterableMatrix object of the projected data.
#' @details DimReduction subclasses should use this to project new data with the same features, to project into the same latent space.
#' All required information to run a projection should be held in x$fitted_params, including pertinent parameters when construction the DimReduction subclass object.
#' If there are no rownames, assume that the matrix is in the same order as the original matrix, assuming that they have the same number of features.
#' If there are rownames, reorder the matrix to match the order of the original matrix
#' @export
project <- function(x, mat, ...) {
  UseMethod("project")
}
#' @export 
project.default <- function(x, mat, ...) {
  rlang::abort("project method not implemented for BPCells objects.")
}
#' @export
project.DimReduction <- function(x, mat, ...) {
  rlang::abort("project method not implemented for base DimReduction object.")
}

#################
# LSI Implementation
#################


#' Perform latent semantic indexing (LSI) on a matrix.
#' 
#' Given a `(features x cells)` matrix, perform LSI to perform tf-idf, z-score normalization, and PCA to create a latent space representation of the matrix of shape `(n_dimensions, ncol(mat))`.
#' @param mat (IterableMatrix) dimensions features x cells.
#' @param n_dimensions (integer) Number of dimensions to keep during PCA.
#' @param corr_cutoff (numeric) Numeric filter for the correlation of a PC to the sequencing depth.  If the PC has a correlation that is great or equal to
#' the corr_cutoff, it will be excluded from the final PCA matrix.
#' @param normalize (function) Normalize matrix using a given function. If `NULL`, no normalization is performed.
#' @param threads (integer) Number of threads to use.
#' @returns An object of class `c("LSI", "DimReduction")` with the following attributes:
#' - `cell_embeddings`: The projected data
#' - `fitted_params`: A tibble of the parameters used for iterative LSI, with rows as iterations. Columns include the following:
#'   - `normalization_method`: The normalization method used
#'   - `feature_means`: The means of the features used for normalization
#'   - `pcs_to_keep`: The PCs that were kept after filtering by correlation to sequencing depth
#'   - `svd_params`: The matrix calculated for SVD
#' - `feature_names`: The names of the features in the matrix
#' @details Compute LSI through first doing a log(tf-idf) transform, z-score normalization, then PCA.  Tf-idf implementation is from Stuart & Butler et al. 2019.
#' 
#' Running on a 2600 cell dataset with 50000 peaks and 4 threads, as an example:
#' - 17.1 MB memory usage, 25.1 seconds runtime
#' @seealso `project()` `DimReduction()` `normalize_tfidf()`
#' @export
LSI <- function(
  mat, n_dimensions = 50L, corr_cutoff = 1, normalize = normalize_tfidf,
  threads = 1L, verbose = FALSE
) {
  if (rlang::is_missing(mat)) {
    return(
      create_partial(
        missing_args = list(
          n_dimensions = missing(n_dimensions),
          corr_cutoff = missing(corr_cutoff), 
          normalize = missing(normalize), 
          threads = missing(threads),
          verbose = missing(verbose)
        )
      )
    )
  }
  assert_is(mat, "IterableMatrix")
  assert_is_wholenumber(n_dimensions)
  assert_len(n_dimensions, 1)
  assert_greater_than_zero(n_dimensions)
  assert_true(n_dimensions < min(ncol(mat), nrow(mat)))
  assert_true((corr_cutoff >= 0) && (corr_cutoff <= 1))
  assert_is_wholenumber(threads)
  
  if (verbose) log_progress("Starting LSI")
  # Normalization
  if (verbose) log_progress("Normalizing matrix")
  mat_stats <- matrix_stats(mat, row_stats = c("mean"), col_stats = c("mean"), threads = threads)
  read_depth <- mat_stats$col_stats["mean", ] * nrow(mat)
  if (!is.null(normalize)) mat <- partial_apply(normalize, threads = threads)(mat)
  
  # Save to prevent re-calculation of queued operations
  mat <- write_matrix_dir(
    convert_matrix_type(mat, type = "float"),
    tempfile("mat"), compress = TRUE
  )
  # Run pca
  if (verbose) log_progress("Calculating SVD")
  svd_attr <- svds(mat, k = n_dimensions, threads = threads)
  pca_res <- t(svd_attr$u) %*% mat
  
  # Filter out PCs that are highly correlated with sequencing depth
  pca_corrs <- abs(cor(read_depth, t(pca_res)))
  pca_feats_to_keep <- which(pca_corrs < corr_cutoff)
  if (length(pca_feats_to_keep) != n_dimensions) {
    if (verbose) log_progress(sprintf("Dropping PCs %s due to high correlation with sequencing depth", paste(setdiff(1:n_dimensions, pca_feats_to_keep), collapse = ", ")))
    pca_res <- pca_res[pca_feats_to_keep, ]
  }
  fitted_params <- list(
    normalization_method = normalize,
    feature_means = mat_stats$row_stats["mean", ],
    pcs_to_keep = pca_feats_to_keep,
    svd_params = svd_attr
  )
  res <- DimReduction(
    x = pca_res,
    fitted_params = fitted_params,
    feature_names = rownames(mat)
  )
  class(res) <- c("LSI", class(res))
  return(res)
}


#' @export
project.LSI <- function(x, mat, threads = 1L, ...) {
  assert_is(mat, "IterableMatrix")
  assert_is(x, "LSI")

  fitted_params <- x$fitted_params
  # Do a check to make sure that the number of rows in the matrix is the same as the number of rows in SVD$u
  assert_true(nrow(mat) == nrow(fitted_params$svd_params$u))
  if (!is.null(rownames(mat)) && !is.null(x$feature_names)) {
    assert_true(all(x$feature_names %in% rownames(mat)))
    mat <- mat[x$feature_names, ]
  }

  if (!is.null(fitted_params$normalization_method))  {
    mat <- fitted_params$normalization_method(
      mat,
      feature_means = fitted_params$feature_means,
      threads = threads
    )
    mat <- write_matrix_dir(
      convert_matrix_type(mat, type = "float"),
      tempfile("mat"), compress = TRUE
    )
  }
  pca_attr <- fitted_params$svd_params
  res <- t(pca_attr$u) %*% mat
  if (length(fitted_params$pcs_to_keep) != nrow(res)) {
    res <- res[fitted_params$pcs_to_keep, ]
  }
  return(res)
}


#' Run iterative LSI on a matrix.
#' 
#' Given a `(features x cells)` matrix, Compute an iterative LSI dimensionality reduction, using the method described in [ArchR](https://doi.org/10.1038/s41588-021-00790-6) (Granja et al; 2019).
#' @param mat (IterableMatrix) 
#' @param n_iterations (int) The number of LSI iterations to perform.
#' @param first_feature_selection_method (function) Method to use for selecting features for the first iteration. Current builtin options are `select_features_by_variance`, `select_features_by_dispersion`, `select_features_by_mean`, `select_features_by_binned_dispersion`
#' @param feature_selection_method (function) Method to use for selecting features for each iteration after the first. Current builtin options are `select_features_by_variance`, `select_features_by_dispersion`, `select_features_by_mean`, `select_features_by_binned_dispersion`
#' @param cluster_method (function) Method to use for clustering. Current builtin options are `cluster_graph_{leiden, louvain, seurat}()`
#' @return An object of class `c("LSI", "DimReduction")` with the following attributes:
#' - `cell_embeddings`: The projected data
#' - `fitted_params`: A tibble of the parameters used for iterative LSI, with rows as iterations. Columns include the following:
#' - `first_feature_selection_method`: The method used for selecting features for the first iteration
#' - `lsi_method`: The method used for LSI
#' - `cluster_method`: The method used for clustering
#' - `feature_means`: The means of the features used for normalization
#' - `iterations`: The number of iterations
#' - `iter_info`: A tibble with the following columns:
#'    - `iteration`: The iteration number
#'    - `feature_names`: The names of the features used for the iteration
#'    - `lsi_results`: The results of LSI for the iteration
#'    - `clusters`: The clusters for the iteration.  This is blank for the first iteration
#' @details
#' The iterative LSI method is as follows:
#' - First iteration:
#'    - Select features based on the `first_feature_selection_method` argument
#'    - Perform LSI on the selected features
#'    - If `n_iterations` is 1, return the PCA results
#'    - Else, cluster the LSI results using `cluster_method`
#' - For each subsequent iteration:
#'    - Pseudobulk the clusters and select the top features based on the variance of the pseudobulked clusters
#'    - Perform LSI on the selected features
#'    - If this is the final iteration, return the PCA results
#'    - Else, cluster the LSI results using `cluster_method`
#' @seealso `LSI()`, `top_features()`, `highly_variable_features()`
#' @inheritParams LSI
#' @export
IterativeLSI <- function(
  mat, 
  n_iterations = 2,
  first_feature_selection_method = select_features_by_binned_dispersion,
  feature_selection_method = select_features_by_dispersion,
  lsi_method = LSI,
  cluster_method = cluster_graph_leiden,
  verbose = FALSE, threads = 1L
) {
  assert_is(mat, "IterableMatrix")
  assert_true(n_iterations > 0)
  assert_is_wholenumber(n_iterations)
  assert_is_wholenumber(threads)
  
  fitted_params = list(
    first_feature_selection_method = first_feature_selection_method,
    lsi_method = lsi_method,
    cluster_method = cluster_method,
    feature_means = matrix_stats(mat, row_stats = "mean", threads = threads)$row_stats["mean",],
    iterations = n_iterations,
    iter_info = tibble::tibble(
      iteration = integer(),
      feature_names = list(),
      lsi_results = list(),
      clusters = list()
    )
  )
  if (verbose) log_progress("Starting Iterative LSI")
  for (i in seq_len(n_iterations)) {
    if (verbose) log_progress(sprintf("Starting Iterative LSI iteration %s of %s", i, n_iterations))
    # add a blank row to the iter_info tibble
    fitted_params$iter_info <- tibble::add_row(fitted_params$iter_info, iteration = i)
    # run variable feature selection
    if (verbose) log_progress("Selecting features")
    if (i == 1) {
      variable_features <- first_feature_selection_method(mat, threads = threads)
    } else {
      variable_features <- feature_selection_method(pseudobulk_res, threads = threads)
    }
    fitted_params$iter_info$feature_names[[i]] <- variable_features %>% dplyr::filter(highly_variable) %>% dplyr::pull(names)
    
    if (is.character(fitted_params$iter_info$feature_names[[i]])) {
      mat_indices <- which(rownames(mat) %in% fitted_params$iter_info$feature_names[[i]])
    } else {
      mat_indices <- fitted_params$iter_info$feature_names[[i]]
    }
    # run LSI
    if (verbose) log_progress("Running LSI")
    lsi_res_obj <- lsi_method(mat[mat_indices,], threads = threads)
    fitted_params$iter_info$lsi_results[[i]] <- lsi_res_obj$fitted_params
    # only cluster + pseudobulk if this isn't the last iteration
    if (i == n_iterations) break
    # cluster the LSI results
    if (verbose) log_progress("Clustering LSI results")
    clustering_res <- t(lsi_res_obj$cell_embeddings) %>% knn_hnsw(ef = 500, threads = threads) %>% knn_to_snn_graph() %>% cluster_method()
    fitted_params$iter_info$clusters[[i]] <- clustering_res
    # pseudobulk and pass onto next iteration
    if (verbose) log_progress("Pseudobulking matrix")
    pseudobulk_res <- pseudobulk_matrix(mat, clustering_res, threads = as.integer(threads)) %>% as("dgCMatrix") %>% as("IterableMatrix")
    rownames(pseudobulk_res) <- rownames(mat)
  }
  if (verbose) log_progress("Finished running LSI")
  res <- DimReduction(
    x = lsi_res_obj$cell_embeddings,
    fitted_params = fitted_params,
    feature_names = rownames(mat)
  )
  class(res) <- c("IterativeLSI", class(res))
  return(res)
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
  assert_is_wholenumber(threads)

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