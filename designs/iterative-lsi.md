# Iterative LSI Design docs

This document describes design choices and required functionality for creating the interface and logic for performing iterative LSI on a dataset.  
Iterative LSI takes in a dataset, then performs feature selection by top-accesibility/variability, then executes a loop of LSI, clustering, and pseudobulked feature selection.  
The last iteration uses the information from the previous iteration's feature selection to perform LSI on the dataset.

To provide the functionality to BPCells, we can utilize the previous work on LSI and highly variable feature selection by binning, as well as creation of a couple new functions to streamline design.  
We also discuss the features and wants of our iterative LSI interface.

## Normalizations

As normalizations during the dimensionality reduction, and feature selection process have been shown to be highly variable (haha!),
we decided to allow for modularity/parameterization in which normalizations are being used.  This would also allow for cleanup of the previous LSI function as well,
so we can stored the tf-idf logic sepearately from the SVD logic.  We also want to allow for the user to pass in pre-computed statistics, such as means, variance, and inverse document frequencies,
to allow for similar operations across multiple datasets.

For the purposes of this document, we will only create `normalize_log()` and `normalize_tfidf()`, as these are the only normalizations that we will use in the iterative LSI process.  We can add more normalizations as needed.

These normalizations should be able to be passed in as a function.  In the case for LSI, we will pass in the `normalize_tfidf()`, which would perform a tf-idf normalization prior to SVD.  We will also want to pass in feature means, so we can have the same normalization across multiple datasets.

### Log normalization

```R
#' Normalize a matrix using log normalization
#' @param mat (IterableMatrix) Matrix to normalize
#' @param scale_factor (numeric) Scale factor for log normalization
#' @returns log normalized matrix.
#' @export
normalize_log <- function(mat, scale_factor = 1e4)
```

### TF-IDF normalization

```R
#' Normalize a matrix using term frequency-inverse document frequency
#' @param mat (IterableMatrix) Matrix to normalize
#' @param feature_means (numeric) Means of the features to normalize by
#' @returns tf-idf normalized matrix.
#' @export
normalize_tfidf <- function(mat, feature_means = NULL)
```

## Feature Selection

Several variable feature selection methods are required for iterative LSI, on top of the highly variable feature selection by binning approach that we have created.  
These methods include feature selection by feature means, variance, and dispersion

We seek to have a shared interface for these methods, so they can be swapped out/parameterized easily.  
Data input should always be a IterableMatrix, and the output should be a tibble with at least the columns of `name`, `score`, and `highly_variable`.  
The `score` column should be the metric used to determine the variable features, and `highly_variable` should be a logical vector of whether the feature is highly variable.  The `name` column should be the name of the feature.  Doing so would allow a user to do the following to filter their original matrix:

Additionally, there is intentional effort to modularize normalizations as functions.  This is developed with the intent to allow for partial arguments to be passed in, and for the function as a whole to be passed into the iterative LSI function as described below.

```R
filtered_df <- features_df %>% dplyr::filter(highly_variable)
feats_of_interest <- which(rownames(mat) %in% filtered_df$name)
mat[feats_of_interest,]
```

### Feature selection by mean

Is self explanatory, and currently used in the first iteration of ArchR.  Uses acessibility scores via sum of non-zeros for each feature.

```R
#' Get the top features from a matrix, based on the mean of each feature.
#' @param num_feats Number of features to deem as highly accessible.  If the number is higher than the number of features in the matrix,
#' all features will be returned.
#' @inheritParams highly_variable_features
#' @returns
#' Return a dataframe with the following columns, sorted descending by variance:
#'   - `names`: Feature name
#'   - `score`: Total accessibility of each feature
#'   - `highly_variable`: Logical vector of whether the feature is highly variable
#' @export
select_features_by_mean <- function(mat, num_feats, threads = 1L)
```

### Variable feature selection by variance

```R
#' Get the most variable features within a matrix.
#' @param num_feats (integer) Number of features to return.  If the number is higher than the number of features in the matrix, 
#' all features will be returned.
#' @param normalize (function) Normalize matrix using a given function.
#' @returns
#' Return a dataframe with the following columns, sorted descending by variance:
#'   - `names`: Feature name
#'   - `score`: Variance of the feature
#'   - `highly_variable`: Logical vector of whether the feature is highly variable
#' @inheritParams highly_variable_features
#' @details 
#' Calculate using the following process:
#'  1. Normalize by term frequency for each pseudobulk
#'  2. Do an optional log normalization of the pseudobulk matrix
#'  3. Find `num_feats` features with the highest variance across clusters
select_features_by_variance <- function(
  mat, num_feats = 25000, 
  normalize = normalize_log, return_subsetted_matrix = FALSE,
  threads = 1L, verbose = TRUE
)
```

### Variable feature selection by dispersion

```R
#' Get the features with the highest dispersion within a matrix
#' @param num_feats (integer) Number of features to return.  If the number is higher than the number of features in the matrix, 
#' all features will be returned.
#' @param normalize (function) Normalize matrix pseudobulked matrix using a given function.
#' @returns
#' Return a dataframe with the following columns, sorted descending by variance:
#'   - `names`: Feature name
#'   - `score`: Variance of the feature
#'   - `highly_variable`: Logical vector of whether the feature is highly variable
#' @inheritParams highly_variable_features
#' @details 
#' Calculate using the following process:
#'  1. Do an optional log normalization of the pseudobulk matrix
#'  2. Find the dispersion (variance/mean) of each feature
#'  3. Find `num_feats` features with the highest variance across clusters
select_features_by_dispersion <- function(
  mat, num_feats = 25000, 
  normalize = normalize_log, return_subsetted_matrix = FALSE,
  threads = 1L, verbose = TRUE
)
```

## Iterative LSI

This iterative LSI approach will deviate from the previous ArchR LSI in a couple ways.  

Logic wise,  In a sense, the first iteration can be just seen as the feature selection and SVD, where the primary goal is to return the post-SVD matrix.  For every iteration after, they will take the previous iteration's feature selection, and perform clustering, pseudobulking, feature selection, and SVD. The iterative LSI function should essentially act as a wrapper, with the logic primarily being the iteration loop.  Within the iteration loop, this function will also perform pseudobulking post clustering, to be fed into variable feature selection.

Styling wise, many of the arguments are removed for the sake of being able to shove them in as functions with partial operations done to them.  This allows for the parameters of iterative LSI to be more clear in terms of which parameters are being used for which operation.  We want to allow for the user to project their data using a pre-used iterative LSI object as well.

The return type of iterative LSI will be a c("IterativeLSI", "DimReduction") object, with the following attributes:

- `cell_embeddings`: The projected data
- `fitted_params`: A tibble of the parameters used for iterative LSI, with rows as iterations. Columns include the following:
  - `iteration`: The iteration number
  - `cluster`: The cluster assignments for each cell
  - `feature_selection`: The tibble of feature selection returned by the feature selection method
  - `feature_means`: The means of the features used for normalization
  - `svd_params`: The matrix calculated for SVD
These feature loadings allow for other datasets to be passed in using the same interface.

```R
#' Run iterative LSI on a matrix.
#' 
#' This function will compute an iterative LSI dimensionality reduction on an ArchRProject.
#' @param mat (IterableMatrix) Matrix to perform LSI on.
#' @param n_iterations The number of LSI iterations to perform.
#' @param first_feature_selection_method (function) First iteration selection method for features to use for LSI. 
#' Refer to `variable_features_by_top_accessibility()` or `variable_features_by_bin_variance()`.
#' @param lsi_method (function) LSI function to use. Refer to `lsi()`.
#' @param cluster_method (function) Clustering function to use. Refer to `cluster_graph_leiden()`, `cluster_graph_louvain()`, or `cluster_graph_seurat()`.
#' @returns An object of class `c("IterativeLSI", "DimReduction")` with the following attributes:
#' - `cell.embeddings`: The projected data
#' - `fitted_params`: A tibble of the parameters used for iterative LSI, with rows as iterations. Columns include the following:
#'   - `iteration`: The iteration number
#'   - `cluster`: The cluster assignments for each cell
#'   - `feature_selection`: The tibble of feature selection returned by the feature selection method
#'   - `feature_means`: The means of the features used for normalization
#'   - `svd_params`: The matrix calculated for SVD
#' @seealso `lsi()`, `variable_features_by_top_accessibility()`, `variable_features_by_bin_variance()`
#' @inheritParams lsi
iterative_lsi <- function(
    mat, 
    n_iterations = 2,
    first_feature_selection_method = variable_features_by_top_accessibility,
    lsi_method = lsi,
    cluster_method = cluster_graph_leiden, 
    threads = 1L
) 
```

### Utilizing functions for arguments

There is a lot of focus of passing functions in rather than traditional parameters into these methods.  The most important reason for this is to allow for partial arguments to be passed into these functions.

As a short demonstration, a user can do this:

```R
iterative_lsi(
    variable_feature = purrr::partial(
        variable_features_variance,
        normalization = purrr::partial(
            normalize_log, 
            scale_factor=100000
        )
    )
)
```

Or even utilize implicit partials if we stuff the purrr::partial into the function:
This allows to hold all the logic of the called functions in one place, allowing for changes in called functions to not break the iterative LSI function from performing as a wrapper.

```R
iterative_lsi(
    variable_feature = variable_features_variance(
        normalization=normalize_log(scale_factor=10000)
    )
)
```
