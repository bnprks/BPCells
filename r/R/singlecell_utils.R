# Copyright 2025 BPCells contributors
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
#' Apply a feature selection method to a non-normalized `(features x cells)` matrix.  We recommend using counts matrices as input and to
#' apply any normalizations prior to feature selection via the normalize argument (if available).  
#' Instead of directly subsetting the input matrix,
#' an output dataframe is returned, indicating which features are highly variable, and the scoring of each feature.
#' @rdname feature_selection
#' @param mat (IterableMatrix) Counts matrix with dimensions `(features x cells)`.
#' @param num_feats (float) Number of features to mark as highly_variable. If 0 < `num_feats` < 1, then interpret as a fraction of features.
#' @param normalize_method (function) Used to normalize the matrix prior to feature selection by calling `normalize_method(mat)` if it is not `NULL`. 
#' For example, pass `normalize_log()` or `normalize_tfidf()`. 
#' @param threads (integer) Number of threads to use. Also overrides the threads argument in `normalize_method`
#' @returns
#' Return a dataframe with the following columns:
#' - `feature`: Feature name.
#' - `score`: Scoring of the feature, depending on the method used.  
#' - `highly_variable`: Logical vector of whether the feature is highly variable.
#'
#' Each different feature selection method will have a different scoring method.
#' Consider a matrix \eqn{X}, where the row index \eqn{i} refers to each feature
#' and the column index \eqn{j} refers to each cell. For each feature \eqn{x_{i} \in X}, we define the following feature-selection scores:
#' - `select_features_variance`: \eqn{\mathrm{Score}(x_i) = \frac{1}{n - 1} \sum_{j=1}^{n} \bigl(x_{ij} - \bar{x}_i\bigr)^2}
#' @examples
#' ## Prep data
#' set.seed(12345)
#' mat <- matrix(rpois(5*8, lambda=1), nrow=5, ncol=8)
#' rownames(mat) <- paste0("gene", seq_len(nrow(mat)))
#' mat
#' mat <- as(mat, "IterableMatrix")
#' 
#' 
#' #######################################################################
#' ## select_features_variance() examples
#' #######################################################################
#' select_features_variance(
#'     mat, 
#'     num_feats=2, 
#'     normalize_method=normalize_log
#' )
#' 
#' # Because of how the BPCells `normalize` functions behave when the matrix 
#' # argument is missing, we can also customize the normalization parameters 
#' # using partial arguments:
#' variable_features <- select_features_variance(
#'     mat,
#'     num_feats=2,
#'     normalize_method=normalize_log(scale_factor=20)
#' ) 
#' # One can then filter to only variable features using the subset operator:
#' mat[variable_features$feature[variable_features$highly_variable],]
#' 
#' 
#' @seealso `normalize_tfidf()` `normalize_log()`
#' @export
select_features_variance <- function(
  mat, num_feats = 0.05, 
  normalize_method = NULL,
  threads = 1L
) {
  assert_greater_than_zero(num_feats)
  assert_len(num_feats, 1)
  assert_is_wholenumber(threads)
  if (rlang::is_missing(mat)) return(create_partial())
  assert_is_mat(mat)
  
  if (num_feats < 1 && num_feats > 0) num_feats <- floor(nrow(mat) * num_feats)
  if (min(max(num_feats, 0), nrow(mat)) != num_feats) {
    rlang::warn(add_timestamp(sprintf("Number of features asked for (%s) is greater than the number of features in the matrix (%s).", num_feats, nrow(mat))))
  }
  num_feats <- min(max(num_feats, 0), nrow(mat))
  if (!is.null(normalize_method)) mat <- partial_apply(normalize_method, threads = threads, .missing_args_error = FALSE)(mat)
  features_df <- tibble::tibble(
    feature = rownames(mat),
    score = matrix_stats(mat, row_stats = "variance", threads = threads)$row_stats["variance",]
  ) %>% 
    dplyr::mutate(highly_variable = dplyr::row_number(dplyr::desc(score)) <= num_feats)
  return(features_df)
}


#' @rdname feature_selection
#' @returns
#' - `select_features_dispersion`: \eqn{\mathrm{Score}(x_i) = \frac{\frac{1}{n - 1} \sum_{j=1}^{n} \bigl(x_{ij} - \bar{x}_i\bigr)^2}{\bar{x}_i}}
#' @examples
#' #######################################################################
#' ## select_features_dispersion() example
#' #######################################################################
#' select_features_dispersion(
#'   mat,
#'   num_feats = 2,
#'   normalize_method = normalize_log
#' )
#' 
#' 
#' @export
select_features_dispersion <- function(
  mat, num_feats = 0.05, 
  normalize_method = NULL,
  threads = 1L
) {
  assert_greater_than_zero(num_feats)
  assert_len(num_feats, 1)
  assert_is_wholenumber(threads)
  if (rlang::is_missing(mat)) return(create_partial())
  if (num_feats < 1 && num_feats > 0) num_feats <- floor(nrow(mat) * num_feats)
  if (min(max(num_feats, 0), nrow(mat)) != num_feats) {
    rlang::warn(add_timestamp(sprintf("Number of features asked for (%s) is greater than the number of features in the matrix (%s).", num_feats, nrow(mat))))
  }
  num_feats <- min(max(num_feats, 0), nrow(mat))
  if (!is(mat, "IterableMatrix") && canCoerce(mat, "IterableMatrix")) mat <- as(mat, "IterableMatrix")
  assert_is(mat, "IterableMatrix")
  if (!is.null(normalize_method)) mat <- partial_apply(normalize_method, threads = threads, .missing_args_error = FALSE)(mat)
  mat_stats <- matrix_stats(mat, row_stats = "variance", threads = threads)
  features_df <- tibble::tibble(
    feature = rownames(mat),
    score = mat_stats$row_stats["variance", ] / mat_stats$row_stats["mean", ]
  ) %>% 
    dplyr::mutate(highly_variable = dplyr::row_number(dplyr::desc(score)) <= num_feats)
  return(features_df)
}

#' @rdname feature_selection
#' @returns
#' - `select_features_mean`: \eqn{\mathrm{Score}(x_i) = \frac{\sum_{j=1}^{n}\bigl(x_{ij}\bigr)}{n}}
#' @examples
#' #######################################################################
#' ## select_features_mean() example
#' #######################################################################
#' select_features_mean(
#'   mat,
#'   num_feats = 1,
#'   normalize_method = normalize_log
#' )
#' 
#' 
#' @export
select_features_mean <- function(mat, num_feats = 0.05, normalize_method = NULL, threads = 1L) {
  assert_greater_than_zero(num_feats)
  assert_is_wholenumber(threads)
  if (rlang::is_missing(mat)) return(create_partial())
  assert_is_mat(mat)
  if (num_feats < 1 && num_feats > 0) num_feats <- floor(nrow(mat) * num_feats)
  if (min(max(num_feats, 0), nrow(mat)) != num_feats) {
    rlang::warn(add_timestamp(sprintf("Number of features asked for (%s) is greater than the number of features in the matrix (%s).", num_feats, nrow(mat))))
  }
  num_feats <- min(max(num_feats, 0), nrow(mat))
  if (!is.null(normalize_method)) mat <- partial_apply(normalize_method, threads = threads, .missing_args_error = FALSE)(mat)
  # get the sum of each feature, binarized
  features_df <- tibble::tibble(
    feature = rownames(mat),
    score = matrix_stats(mat, row_stats = "mean", threads = threads)$row_stats["mean", ]
  ) %>%
    dplyr::mutate(highly_variable = dplyr::row_number(dplyr::desc(score)) <= num_feats)
  return(features_df)
}

#' @rdname feature_selection
#' @param n_bins (integer) Number of bins to split features into in order to control for the relationship between mean expression and dispersion (see details).
#' @returns
#'  - `select_features_binned_dispersion`: Process described in `details`.
#' @details 
#' `select_features_binned_dispersion` implements the approach from Satija et al. 2015:
#'  1. Bin features into equal-width bins by `log1p(mean)`
#'  2. Calculate dispersion of each feature as `log(variance / mean)`
#'  3. Z-score normalize dispersion within each bin, and select highest normalized dispersion across all bins
#' 
#' If the number of features within a bin is equal to 1, then the mean dispersion for that bin is set to 1.
#' 
#' This should be equivalent to `Seurat::FindVariableFeatures()` with `selection.method="mean.var.plot"`
#'  and `scanpy.pp.highly_variable_genes()` with `flavor="seurat"`.
#' @examples
#' #######################################################################
#' ## select_features_binned_dispersion() example
#' #######################################################################
#' select_features_binned_dispersion(
#'   mat,
#'   num_feats = 2,
#'   n_bins = 2
#' )
#' 
#' 
#' @export
select_features_binned_dispersion <- function(
  mat, num_feats = 0.05, n_bins = 20,
  threads = 1L
) {
  assert_greater_than_zero(num_feats)
  assert_len(num_feats, 1)
  assert_is_wholenumber(n_bins)
  assert_len(n_bins, 1)
  assert_greater_than_zero(n_bins)
  assert_is_wholenumber(threads)
  if (rlang::is_missing(mat)) return(create_partial())
  assert_is_mat(mat)
  if (num_feats < 1 && num_feats > 0) num_feats <- floor(nrow(mat) * num_feats)
  if (min(max(num_feats, 0), nrow(mat)) != num_feats) {
    rlang::warn(add_timestamp(sprintf("Number of features asked for (%s) is greater than the number of features in the matrix (%s).", num_feats, nrow(mat))))
  }
  num_feats <- min(max(num_feats, 0), nrow(mat))
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
  features_df <- tibble::tibble(
    feature = names(feature_means),
    var = feature_vars,
    mean = feature_means,
    dispersion = feature_dispersion
  )
  # Bin by mean, and normalize dispersion with each bin
  features_df <- features_df %>%
    dplyr::mutate(bin = cut(mean, n_bins, labels = FALSE)) %>% 
    dplyr::group_by(bin) %>%
    dplyr::mutate(
      score = (dispersion - mean(dispersion)) / sd(dispersion),
      score = if (dplyr::n() == 1) {1} else {score} # Set feats that are in bins with only one feat to have a norm dispersion of 1  
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(highly_variable = dplyr::row_number(dplyr::desc(score)) <= num_feats) %>% 
    dplyr::select(c("feature", "score", "highly_variable", "dispersion", "bin")) %>%
    dplyr::rename("raw_log_dispersion" = "dispersion")
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
#' @rdname DimReduction
#' @field cell_embeddings (IterableMatrix, dgCMatrix, matrix) Projected data of shape `(cells x n_dimensions)` of the original matrix after a dimensionality reduction.
#' @field fitted_params (list) A list of parameters used for the transformation of a matrix.  This should include all necessary information to project new data with the same features.
#' @field feature_names (character) The names of the features that this DimReduction object was fit on.  Matrices to be projected should have the same feature names.
#' @return - `DimReduction()`: DimReduction object.
#' @export
DimReduction <- function(cell_embeddings, fitted_params = list(), feature_names = character(0), ...) {
  assert_is(cell_embeddings, c("IterableMatrix", "dgCMatrix", "matrix"))
  assert_is(fitted_params, "list")
  structure(list(
    cell_embeddings = cell_embeddings,
    fitted_params = fitted_params,
    feature_names = feature_names,
    ...
    ),
    class = "DimReduction"
  )
}

#' @rdname DimReduction
#' @param (IterableMatrix) Input matrix of shape `(features x cells)`.
#' @param x DimReduction object.
#' @return - `project()`: IterableMatrix object of the projected data.
#' @details 
#' `project()`: Perform a dimensionality reduction on a matrix using a pre-fit DimReduction object.
#' 
#' DimReduction subclasses should use the `project` method on new data with the same features, to project into the same latent space.
#' All required information to run a projection should be held in `x$fitted_params`, including pertinent parameters when constructing the DimReduction subclass object.
#' 
#' If there are no rownames, assume that the matrix is in the same order as the original matrix, assuming that they have the same number of features.
#' If there are rownames, reorder the matrix to match the order of the original matrix
#' @inheritParams DimReduction
#' @export
project <- function(x, mat, ...) {
  UseMethod("project")
}
#' @export 
project.default <- function(x, mat, ...) {
  rlang::abort("project method not implemented for objects that are not a fitted DimReduction")
}

#' @export
print.DimReduction <- function(x, ...) {
  cat(sprintf("Fitted <%s> dimensionality reduction\n\n", class(x)[1]))
  
  # Print feature info
  cat("Number of features:", length(x$feature_names), "\n") 
  cat("Input feature names:", pretty_print_vector(x$feature_names), "\n")

  # Print embedding info
  dim_embeddings <- dim(x$cell_embeddings)
  cat(sprintf("cell_embeddings: %d cells x %d embedding dimensions\n", dim_embeddings[1], dim_embeddings[2]))

  # Print param info
  # params_str <- paste(names(x$fitted_params), collapse = ", ")
  # wrapped_params <- strwrap(params_str, width = getOption("width") - 6)
  # cat("Fitted_params:\n")
  # for (ln in wrapped_params) {
  #   cat("  ", ln, "\n", sep = "")
  # }
  cat("fitted_params: ", stringr::str_wrap(paste(names(x$fitted_params), collapse = ", "), exdent = 2, width = 60), "\n")
}

#################
# LSI Implementation
#################


#' Perform latent semantic indexing (LSI) on a matrix.
#' 
#' Given a `(features x cells)` counts matrix, perform LSI, which sequentially executes tf-idf normalization and PCA to create a latent space representation 
#' of the matrix of shape `(n_dimensions, ncol(mat))`.
#' Returns a DimReduction object, which allows for projection of new matrices with the same features into the same latent space.
#' @rdname LSI
#' @param mat (IterableMatrix) Counts matrix of shape `(features x cells)`.
#' @param n_dimensions (integer) Number of dimensions to keep during PCA.
#' @param corr_cutoff (numeric) Numeric filter for the correlation of a PC to the sequencing depth.  If the PC has a correlation that is greater or equal to
#' the corr_cutoff, it will be excluded from the final PCA matrix.
#' @param scale_factor (numeric) Scaling factor to multiply matrix by during tf-idf normalization.
#' @param threads (integer) Number of threads to use.
#' @returns 
#' `LSI()` An object of class `c("LSI", "DimReduction")` with the following attributes:
#' - `cell_embeddings`: The projected data as a matrix of shape `(cells, n_dimensions)`
#' - `fitted_params`: A named list of the parameters used for LSI, including the following:
#'   - `scale_factor`: The scale factor used for tf-idf normalization
#'   - `feature_means`: The means of the features used for normalization
#'   - `pcs_to_keep`: The PCs that were kept after filtering by correlation to sequencing depth
#'   - `feature_loadings`: SVD component with dimension `(n_variable_features, n_dimensions)` used to linearly transform input data into cell embeddings
#' - `feature_names`: The names of the features in the matrix
#' @details Compute LSI through first doing a log(tf-idf) transform, z-score normalization, then PCA.  Tf-idf implementation is from Stuart & Butler et al. 2019.
#' 
#' Running on a 2600 cell dataset with 50000 peaks and 4 threads, as an example:
#' - 17.1 MB memory usage, 25.1 seconds runtime
#' @seealso `project()` `DimReduction()` `normalize_tfidf()` `normalize_log()` `svds()`
#' @examples
## Prep data
#' nrows <- 50
#' ncols <- 1000
#' mat <- matrix(1:(nrows*ncols), nrow = nrows) %>% as("IterableMatrix")
#' rownames(mat) <- paste0("feat", seq(nrows))
#' colnames(mat) <- paste0("cell", seq(ncols))
#' 
#' 
#' #######################################################################
#' ## LSI() example 
#' #######################################################################
#' lsi_result <- LSI(mat, n_dimensions = 10)
#' lsi_result
#' 
#' 
#' @export
LSI <- function(
  mat, n_dimensions = 50L, corr_cutoff = 1, scale_factor = 1e4,
  threads = 1L, verbose = FALSE
) {
  if (rlang::is_missing(mat)) {
    return(create_partial())
  }
  assert_is_mat(mat)
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
  mat <- normalize_tfidf(
    mat = mat,
    feature_means = mat_stats$row_stats["mean", ],
    scale_factor = scale_factor,
    threads = threads
  )
  # Save to prevent repeating queued calculations during SVD
  temp_mat_dir <- tempfile("mat")
  mat <- write_matrix_dir(
    convert_matrix_type(mat, type = "float"),
    temp_mat_dir, compress = TRUE
  )
  on.exit(unlink(temp_mat_dir, recursive = TRUE, expand = FALSE))
  # Run pca
  if (verbose) log_progress("Calculating SVD")
  svd_attr <- svds(mat, k = n_dimensions, threads = threads)
  pca_res <- svd_attr$v %*% diag(svd_attr$d)
  rownames(pca_res) <- colnames(mat)
  
  # Filter out PCs that are highly correlated with sequencing depth
  pca_corrs <- abs(cor(read_depth, pca_res))
  pca_feats_to_keep <- which(pca_corrs < corr_cutoff)
  if (length(pca_feats_to_keep) != n_dimensions) {
    if (verbose) log_progress(sprintf("Dropping PCs %s due to high correlation with sequencing depth", paste(setdiff(1:n_dimensions, pca_feats_to_keep), collapse = ", ")))
    pca_res <- pca_res[, pca_feats_to_keep, drop = FALSE]
  }
  fitted_params <- list(
    scale_factor = scale_factor,
    feature_means = mat_stats$row_stats["mean", ],
    pcs_to_keep = pca_feats_to_keep,
    feature_loadings = svd_attr$u
  )
  res <- DimReduction(
    cell_embeddings = pca_res,
    fitted_params = fitted_params,
    feature_names = rownames(mat)
  )
  class(res) <- c("LSI", class(res))
  return(res)
}
#' @rdname LSI
#' @return 
#' `project()` IterableMatrix of the projected data of shape `(n_dimensions, ncol(mat))`.
#' @inheritParams project
#' @examples
#' #######################################################################
#' ## project(<LSI>) example 
#' #######################################################################
#' dim(project(lsi_result, mat))
#' 
#' 
#' @export
project.LSI <- function(x, mat, threads = 1L, ...) {
  assert_is_mat(mat)
  assert_is(x, "LSI")

  fitted_params <- x$fitted_params
  # Do a check to make sure that the number of rows in the matrix is the same as the number of rows in SVD$u
  assert_true(nrow(mat) == nrow(fitted_params$feature_loadings))
  if (!is.null(rownames(mat)) && !is.null(x$feature_names)) {
    assert_true(all(x$feature_names %in% rownames(mat)))
    mat <- mat[x$feature_names, ]
  }
  mat <- normalize_tfidf(
    mat,
    feature_means = fitted_params$feature_means,
    scale_factor = fitted_params$scale_factor,
    threads = threads
  )
  feature_loadings <- fitted_params$feature_loadings
  res <- t(mat) %*% feature_loadings
  if (length(fitted_params$pcs_to_keep) != ncol(res)) {
    res <- res[, fitted_params$pcs_to_keep, drop = FALSE]
  }
  return(res)
}



#' Run Iterative LSI on a matrix.
#' 
#' Given a `(features x cells)` counts matrix, perform IterativeLSI to create a latent space representation of the matrix of shape `(n_dimensions, ncol(mat))`.  
#' This uses the method described in [ArchR](https://doi.org/10.1038/s41588-021-00790-6) (Granja et al; 2019). 
#' See details for more specific information.  Returns a DimReduction object, which allows for projection of matrices with the same features into the same latent space.
#' @rdname IterativeLSI
#' @param mat (IterableMatrix) Counts matrix of shape `(features x cells)`.
#' @param n_iterations (int) The number of LSI iterations to perform.
#' @param feature_selection_method (function) Method to use for selecting features for LSI. 
#' Current builtin options are `select_features_variance`, `select_features_dispersion`, `select_features_mean`, `select_features_binned_dispersion`
#' @param cluster_method (function) Method to use for clustering the post-SVD matrix. 
#' The user can pass in partial parameters to the cluster method, such as by passing 
#' `cluster_cells_graph(mat, graph_to_cluster_method = cluster_graph_louvain(resolution = 0.5))` into `cluster_method`.
#' @param threads (integer) Number of threads to use.  Also gets passed down into `feature_selection_method` and `cluster_method`
#' @return 
#' `IterativeLSI()` An object of class `c("IterativeLSI", "DimReduction")` with the following attributes:
#' - `cell_embeddings`: The projected data as a matrix of shape `(cells, n_dimensions)`
#' - `fitted_params`: A list of the parameters used for iterative LSI.  Includes the following:
#'    - `lsi_method`: The method used for LSI
#'    - `feature_selection_method`: The method used for selecting features
#'    - `cluster_method`: The method used for clustering
#'    - `feature_means`: The means of the features used for tf-idf normalization
#'    - `iterations`: The number of LSI iterations ran
#'    - `iter_info`: A tibble of iteration info with rows as iterations. Columns include the following:
#'        - `iteration`: The iteration number
#'        - `feature_names`: The names of the features used for the iteration
#'        - `feature_loadings`: SVD component `u` with dimension `(n_dimensions, n_variable_features)`
#'        - `pcs_to_keep`: The PCs that were kept after filtering by correlation to sequencing depth
#'        - `clusters`: The clusters for the iteration.  This is blank for the first iteration
#' @details
#' The Iterative LSI method is as follows:
#' - First iteration:
#'    - Select features using `feature_selection_method`
#'    - Perform LSI on the selected features
#'    - If `n_iterations` is 1, return the projected data from the first PCA projection
#'    - Else, cluster the LSI results using `cluster_method`
#' - For each subsequent iteration:
#'    - Pseudobulk the clusters from the previous iteration and select the top features based on the pseudobulked clusters
#'    - Perform LSI on the selected features
#'    - If this is the final iteration, return the projected data from this PCA projection
#'    - Else, cluster the LSI results using `cluster_method`
#' 
#' There are some minor differences when compared to the ArchR implementation:
#' - ArchR binarizes data prior to feature selection.  To replicate this, the user can pass `select_features_variance(normalize=binarize)` for their `feature_selection_method`.
#' - `IterativeLSI()` currently does not support utilization of different feature selection methods across each iteration. 
#' If one desires to use a different feature selection method for each iteration, they can take the cluster assignments from the previous 
#' iteration and use them to select features and run LSI.
#' - ArchR uses a default of 25000 features picked during feature selection. As the number of input features is dependent on the input matrix fed into `IterativeLSI()`, 
#' the default for `select_features_variance()` instead picks the number of variable features as a proportion of the total features provided.  To mimic the ArchR implementation,
#' `feature_selection_method` can be set to `select_features_variance(num_feats = 25000)`.
#' - ArchR calculates LSI during non-terminal iterations using a default subset of 10000 cells.  ArchR does this to prevent a memory bottleneck,
#' which BPCells does not encounter even with a non-subsetted matrix. Therefore, IterativeLSI will run LSI on the entire matrix for each iteration.
#' - ArchR defaults on using Seurat clustering for default, which utilizes the Louvain algorithm (See `Seurat::FindClusters()`).  In constrast, `IterativeLSI()` utilizes
#' leiden, which should provide the same clustering results while being faster.
#' - ArchR also plots and calculates a umap of every iteration's dimensionality reduction.  While this is not implemented in `IterativeLSI()`,
#' one can use the `project()` method with the `iteration` argument set to the desired iteration to get projected data.  This can then be fed into `uwot::umap()`
#' - ArchR filters out PCs with a correlation to sequencing depth greater than 0.75. 
#' While corr_cutoff is provided as an argument in `IterativeLSI()`, it is set to not removing any PCs by default.
#' - ArchR filters out outliers dependent on number of accesible regions of cells, by the bottom and top quantiles.  This is not implemented in `IterativeLSI()`, 
#' but can be done as a preprocessing step.
#' @seealso `LSI()` `DimReduction()` `svds()`
#' `cluster_cells_graph()` `select_features_variance()` `select_features_dispersion()` 
#' `select_features_mean()` `select_features_binned_dispersion()`
#' @inheritParams LSI
#' @examples
#' ## Prep data
#' nrows <- 350
#' ncols <- 2000
#' mat <- matrix(1:(nrows*ncols), nrow = nrows) %>% as("IterableMatrix")
#' rownames(mat) <- paste0("feat", seq(nrows))
#' colnames(mat) <- paste0("cell", seq(ncols))
#' 
#' 
#' #######################################################################
#' ## IterativeLSI() examples
#' #######################################################################
#' dim_reduction <- IterativeLSI(mat, n_dimensions = 5)
#' 
#' ## Can customize parameters using partialization
#' dim_reduction <- IterativeLSI(
#'   mat,
#'   n_dimensions = 10,
#'   feature_selection_method = select_features_variance(
#'     num_feats = 0.5,
#'     normalize_method = normalize_tfidf(scale_factor = 5000)
#'   ),
#'   cluster_method = cluster_cells_graph(
#'     graph_to_cluster_method = cluster_graph_louvain(resolution = 0.5),
#'     knn_to_graph_method = knn_to_snn_graph
#'   )
#' )
#' dim_reduction
#' 
#' 
#' @export
IterativeLSI <- function(
  mat, 
  n_iterations = 2,
  feature_selection_method = select_features_variance,
  scale_factor = 1e4, 
  n_dimensions = 50L, 
  corr_cutoff = 1,
  cluster_method = cluster_cells_graph,
  threads = 1L, verbose = FALSE
) {
  assert_has_package("RcppHNSW")
  assert_is_mat(mat)
  assert_greater_than_zero(n_iterations)
  assert_is_wholenumber(n_iterations)
  assert_is_wholenumber(threads)
  assert_is(verbose, "logical")
  fitted_params <- list(
    lsi_method = LSI(n_dimensions = n_dimensions, corr_cutoff = corr_cutoff, scale_factor = scale_factor, threads = threads),
    cluster_method = cluster_method,
    feature_selection_method = feature_selection_method,
    feature_means = matrix_stats(mat, row_stats = "mean", threads = threads)$row_stats["mean",],
    iterations = n_iterations,
    iter_info = tibble::tibble(
      iteration = integer(),
      feature_names = list(),
      feature_loadings = list(),
      pcs_to_keep = list(),
      clusters = list()
    )
  )
  feature_selection_method <- partial_apply(feature_selection_method, threads = threads, .missing_args_error = FALSE)
  if (verbose) log_progress("Starting Iterative LSI")
  for (i in seq_len(n_iterations)) {
    if (verbose) log_progress(sprintf("Starting Iterative LSI iteration %s of %s", i, n_iterations))
    # run variable feature selection
    if (verbose) log_progress("Selecting features")
    if (i == 1) {
      var_feature_table <- feature_selection_method(mat)
    } else {
      var_feature_table <- feature_selection_method(pseudobulk_res)
    }
    variable_features <- var_feature_table %>% dplyr::filter(highly_variable) %>% dplyr::pull(feature)
    if (is.character(variable_features)) {
      mat_indices <- which(rownames(mat) %in% variable_features)
    } else {
      mat_indices <- variable_features
    }
    if (length(mat_indices) < n_dimensions) {
      rlang::abort(
        sprintf(paste0(
          "Number of variable features selected (%s) is smaller than number of dimensions requested (%s). \n",
          "Try setting the num_feats arg in feature_selection_method as a larger value."), 
          length(mat_indices), n_dimensions
        )
      )
    }
    # run LSI
    if (verbose) log_progress("Running LSI")

    lsi_res_obj <- LSI(
      mat = mat[mat_indices, ],
      n_dimensions = n_dimensions,
      corr_cutoff = corr_cutoff,
      scale_factor = scale_factor,
      threads = threads,
      verbose = verbose
    )
    fitted_params$iter_info <- tibble::add_row(
      fitted_params$iter_info, iteration = i,
      feature_names = list(variable_features),
      feature_loadings = list(lsi_res_obj$fitted_params$feature_loadings),
      pcs_to_keep = list(lsi_res_obj$fitted_params$pcs_to_keep)
    )
    # only cluster + pseudobulk if this isn't the last iteration
    if (i == n_iterations) break
    
    # cluster the LSI results
    if (verbose) log_progress("Clustering LSI results")
    clustering_res <- lsi_res_obj$cell_embeddings[, lsi_res_obj$fitted_params$pcs_to_keep, drop = FALSE] %>%
      partial_apply(cluster_method, threads = threads, .missing_args_error = FALSE)()
    fitted_params$iter_info$clusters[[i]] <- clustering_res
    # pseudobulk and pass onto next iteration
    if (verbose) log_progress("Pseudobulking matrix")
    pseudobulk_res <- pseudobulk_matrix(mat, clustering_res, threads = threads)
    rownames(pseudobulk_res) <- rownames(mat)
  }
  if (verbose) log_progress("Finished running Iterative LSI")
  res <- DimReduction(
    cell_embeddings = lsi_res_obj$cell_embeddings,
    fitted_params = fitted_params,
    feature_names = rownames(mat)
  )
  class(res) <- c("IterativeLSI", class(res))
  return(res)
}
#' @rdname IterativeLSI
#' @param iteration (integer) Which iteration of `IterativeLSI`'s features and loadings to use for projection.
#' @return
#' `project()` Matrix of the projected data of shape `(cells, n_dimensions)`.
#' @inheritParams project
#' @examples
#' #######################################################################
#' ## project(<IterativeLSI>) example
#' #######################################################################
#' dim(project(dim_reduction, mat))
#' 
#' 
#' @export
project.IterativeLSI <- function(x, mat, iteration = x$fitted_params$iterations, threads = 1L, ...) {
  assert_is_mat(mat)
  assert_is_wholenumber(threads)
  assert_is_wholenumber(iteration)
  fitted_params <- x$fitted_params
  # Get the desired iteration row of iter_info tibble
  assert_true(iteration <= x$fitted_params$iterations)
  iter_info <- fitted_params$iter_info[iteration, ]

  # Do a check to make sure that the fitted features all exist in input matrix
  if (!is.null(rownames(mat)) && !is.null(x$feature_names)) {
    assert_true(all(x$feature_names %in% rownames(mat)))
  }
  # Subset to variable features
  if (is.character(iter_info$feature_names[[1]])) {
    mat_indices <- which(rownames(mat) %in% iter_info$feature_names[[1]])
  } else {
    mat_indices <- iter_info$feature_names[[1]]
  }
  mat <- mat[mat_indices,]
  # Run LSI
  # since we don't hold the LSI object, copy the internal logic from `project.LSI()`
  lsi_attr <- attr(x$fitted_params$lsi_method, "args")

  mat <- normalize_tfidf(
    mat = mat,
    feature_means = fitted_params$feature_means,
    scale_factor = lsi_attr$scale_factor,
    threads = threads
  )

  feature_loadings <- iter_info$feature_loadings[[1]]
  res <-  t(mat) %*% feature_loadings
  if (length(iter_info$pcs_to_keep[[1]]) != ncol(res)) {
    res <- res[, iter_info$pcs_to_keep[[1]], drop = FALSE]
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