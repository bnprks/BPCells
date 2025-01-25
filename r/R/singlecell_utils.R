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