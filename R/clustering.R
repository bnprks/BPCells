
#' K Nearest Neighbor (KNN) Graph
#'
#' Convert a KNN object (e.g. returned by `knn_hnsw()` or `knn_annoy()`) into
#' a graph. The graph is represented as a sparse adjacency matrix.
#'
#' @rdname knn_graph
#' @details **knn_to_graph**
#' Create a knn graph
#' @param knn List of 2 matrices -- idx for cell x K neighbor indices,
#'  dist for cell x K neighbor distances
#' @param use_weights boolean for whether to replace all distance weights with 1
#' @param self_loops boolean for whether to allow cells to count themselves as neighbors
#' @return **knn_to_graph**
#'   Sparse matrix (dgCMatrix) where `mat[i,j]` = distance from cell `i` to
#'   cell `j`, or 0 if cell `j` is not in the K nearest neighbors of `i`
knn_to_graph <- function(knn, use_weights = FALSE, self_loops = TRUE) {
  if (use_weights) {
    weights <- knn$dist
  } else {
    weights <- 1
  }
  mat <- Matrix::sparseMatrix(
    i = rep.int(seq_len(nrow(knn$idx)), ncol(knn$idx)),
    j = knn$idx,
    x = weights
  )
  if (!self_loops) {
    diag(mat) <- 0
    mat <- Matrix::drop0(mat)
  }
  rownames(mat) <- rownames(knn$idx)
  colnames(mat) <- rownames(knn$idx)
  mat
}


#' @rdname knn_graph
#' @details **knn_to_snn_graph**
#'  Convert a knn object into a shared nearest neighbors adjacency matrix.
#'  This follows the algorithm that Seurat uses to compute SNN graphs
#' @param min_val minimum jaccard index between neighbors. Values below this will
#'  round to 0
#' @param self_loops Whether to allow self-loops in the output graph
#' @param return_type Whether to return a sparse adjacency matrix or an edge list
#' @return **knn_to_snn_graph**
#' - `return_type == "matrix"`:
#'      Sparse matrix (dgCMatrix) where `mat[i,j]` = jaccard index of the overlap
#'      in nearest neigbors between cell `i` and cell `j`, or 0 if the jaccard index
#'      is < `min_val`. Only the lower triangle is filled in, which is compatible with
#'      the BPCells clustering methods
#' - `return_type == "list"`:
#'      List of 3 equal-length vectors `i`, `j`, and `weight`, along with an integer `dim`. 
#'      These correspond to the rows, cols, and values of non-zero entries in the lower triangle
#'      adjacency matrix. `dim` is the total number of vertices (cells) in the graph 
#' @export
knn_to_snn_graph <- function(knn, min_val = 1 / 15, self_loops = FALSE, return_type=c("matrix", "list")) {
  return_type <- match.arg(return_type)
  # Solve x / (2*K - x) >= min_val --> x >= 2*K*min_val / (1 + min_val)
  min_int <- ceiling(2*min_val*ncol(knn$idx) / (1 + min_val))
  snn <- build_snn_graph_cpp(knn$idx, min_neighbors = min_int)

  snn$snn_count <- snn$snn_count / (2*ncol(knn$idx) - snn$snn_count)
  
  if (!self_loops) {
    keepers <- snn$i != snn$j
    snn$i <- snn$i[keepers]
    snn$j <- snn$j[keepers]
    snn$snn_count <- snn$snn_count[keepers]
  }

  names(snn)[names(snn) == "snn_count"] <- "weight"
  snn$dim <- nrow(knn$idx)
  if (return_type == "list") {
    return(snn)
  }
  
  # Return as a sparse matrix
  Matrix::sparseMatrix(
    i = snn$i + 1L, j = snn$j + 1L, x = snn$weight,
    dims = c(snn$dim, snn$dim),
    repr="C"
  )
}

#' @rdname knn_graph
#' @details **knn_to_geodesic_graph**
#'  Convert a knn object into an undirected weighted graph, using the same 
#'  geodesic distance estimation method as the UMAP package. 
#'  This matches the output of `umap._umap.fuzzy_simplicial_set`
#'  from the `umap-learn` python package, used by default in `scanpy.pp.neighbors`.
#'  Because this only re-weights and symmetrizes the KNN graph, it will usually use
#'  less memory and return a sparser graph than `knn_to_snn_graph` which computes
#'  2nd-order neighbors. Note: when cells don't have themselves listed as the nearest
#'  neighbor, results may differ slightly from `umap._umap.fuzzy_simplicial_set`, which
#'  assumes self is always successfully found in the approximate nearest neighbor search.
#'  
#' @param threads Number of threads to use during calculations
#' @return **knn_to_geodesic_graph**
#' - `return_type == "matrix"`:
#'      Sparse matrix (dgCMatrix) where `mat[i,j]` = normalized similarity between cell `i` and cell `j`.
#'      Only the lower triangle is filled in, which is compatible with the BPCells clustering methods
#' - `return_type == "list"`:
#'      List of 3 equal-length vectors `i`, `j`, and `weight`, along with an integer `dim`. 
#'      These correspond to the rows, cols, and values of non-zero entries in the lower triangle
#'      adjacency matrix. `dim` is the total number of vertices (cells) in the graph 
#' @export
knn_to_geodesic_graph <- function(knn, return_type=c("matrix", "list"), threads=0L) {
  return_type <- match.arg(return_type)
  graph <- build_umap_graph_cpp(knn$dist, knn$idx)
  
  graph$dim <- nrow(knn$idx)
  if (return_type == "list") {
    return(graph)
  }
  
  # Return as a sparse matrix
  Matrix::sparseMatrix(
    i = graph$i + 1L, j = graph$j + 1L, x = graph$weight,
    dims = c(graph$dim, graph$dim),
    repr="C"
  )
}

#' Cluster an adjacency matrix
#' @rdname cluster
#' @details **cluster_graph_leiden**: Leiden graph clustering algorithm `igraph::cluster_leiden()`
#' @param snn Symmetric adjacency matrix (dgCMatrix) output from e.g. knn_to_snn_graph. Only the lower triangle is used
#' @param resolution Resolution parameter. Higher values result in more clusters
#' @param seed Random seed for clustering initialization
#' @param ... Additional arguments to underlying clustering function
#' @return Factor vector containing the cluster assignment for each cell.
#' @export
cluster_graph_leiden <- function(snn, resolution = 1e-3, seed = 12531, ...) {
  # Set seed without permanently changing seed state
  prev_seed <- get_seed()
  on.exit(restore_seed(prev_seed), add = TRUE)
  set.seed(seed)

  igraph::graph_from_adjacency_matrix(snn, weighted = TRUE, diag = FALSE, mode = "lower") %>%
    igraph::cluster_leiden(resolution_parameter = resolution, ...) %>%
    igraph::membership() %>%
    as.factor()
}


#' @rdname cluster
#' @details **cluster_graph_louvain**: Louvain graph clustering algorithm `igraph::cluster_louvain()`
#' @export
cluster_graph_louvain <- function(snn, resolution = 1, seed = 12531) {
  # Set seed without permanently changing seed state
  prev_seed <- get_seed()
  on.exit(restore_seed(prev_seed), add = TRUE)
  set.seed(seed)

  igraph::graph_from_adjacency_matrix(snn, weighted = TRUE, diag = FALSE, mode = "lower") %>%
    igraph::cluster_louvain(resolution = resolution) %>%
    igraph::membership() %>%
    as.factor()
}

#' @rdname cluster
#' @details **cluster_graph_seurat**: Seurat's clustering algorithm `Seurat::FindClusters()`
#' @export
cluster_graph_seurat <- function(snn, resolution = 0.8, ...) {
  assert_has_package("Seurat")
  Seurat::as.Graph(snn) %>%
    Seurat::FindClusters(resolution = resolution, ...) %>%
    .[[1]]
}

#' Convert grouping vector to sparse matrix
#'
#' Converts a vector of membership IDs into a sparse matrix
#'
#' @param groups Vector with one entry per cell, specifying the cell's group
#' @param group_order Optional vector listing ordering of groups
#' @return cell x group matrix where an entry is 1 when a cell is in a given group
cluster_membership_matrix <- function(groups, group_order = NULL) {
  if (is.null(group_order)) {
    group_order <- levels(as.factor(groups))
  }
  membership <- match(as.character(groups), group_order)

  if (anyNA(membership)) {
    rlang::warn(sprintf(
      "Could not match %d groups to group_order: %s",
      sum(is.na(membership)),
      pretty_print_vector(groups[is.na(membership)])
    ))
  }
  Matrix::sparseMatrix(
    i = seq_along(groups)[!is.na(membership)],
    j = membership[!is.na(membership)],
    x = 1,
    dims = c(length(groups), length(group_order)),
    dimnames = list(names(groups), group_order)
  )
}


#' Get a knn matrix from reduced dimensions
#'
#' Search for approximate nearest neighbors between cells in the reduced
#' dimensions (e.g. PCA), and return the k nearest neighbors (knn) for each
#' cell. Optionally, we can find neighbors between two separate sets of cells by
#' utilizing both data and query.
#'
#' @rdname knn
#' @details **knn_hnsw**: Use RcppHNSW as knn engine
#' @param data cell x dims matrix for reference dataset
#' @param query cell x dims matrix for query dataset (optional)
#' @param k number of neighbors to calculate
#' @param metric distance metric to use
#' @param threads Number of threads to use. Note that result is non-deterministic
#'          if threads > 1
#' @param ef ef parameter for RccppHNSW::hnsw_search. Increase for slower search but
#'          improved accuracy
#' @param verbose whether to print progress information during search
#' @return List of 2 matrices -- idx for cell x K neighbor indices,
#'         dist for cell x K neighbor distances.
#'         If no query is given, nearest neighbors are found mapping
#'         the data matrix to itself, prohibiting self-neighbors
#' @export
knn_hnsw <- function(data, query = NULL, k = 10, metric = c("euclidean", "cosine"), verbose = TRUE, threads = 1, ef = 100) {
  metric <- match.arg(metric)
  index <- RcppHNSW::hnsw_build(
    data,
    distance = metric,
    verbose = verbose,
    n_threads = threads
  )
  no_query <- is.null(query)
  if (no_query) query <- data
  res <- RcppHNSW::hnsw_search(
    query,
    index,
    k,
    ef = ef,
    verbose = verbose,
    n_threads = threads
  )
  if (no_query) {
    missed_searches <- sum(res$idx[, 1] != seq_len(nrow(data)))
    if (missed_searches > 0) {
      warning(sprintf("KNN search didn't find self-neighbor for %d datapoints. Try higher ef value", missed_searches))
    }
  }

  rownames(res$idx) <- rownames(data)
  rownames(res$dist) <- rownames(data)
  return(res)
}

#' @rdname knn
#' @details **knn_annoy**: Use RcppAnnoy as knn engine
#' @param n_trees Number of trees during index build time. More trees gives higher accuracy
#' @param search_k Number of nodes to inspect during the query, or -1 for default value. Higher number gives higher accuracy
#' @export
knn_annoy <- function(data, query = data, k = 10, metric = c("euclidean", "cosine", "manhattan", "hamming"), n_trees = 50, search_k = -1) {
  metric <- match.arg(metric)
  annoy <- switch(metric,
    "euclidean" = new(RcppAnnoy::AnnoyEuclidean, ncol(data)),
    "cosine" = new(RcppAnnoy::AnnoyAngular, ncol(data)),
    "manhattan" = new(RcppAnnoy::AnnoyManhattan, ncol(data)),
    "hamming" = new(RcppAnnoy::AnnoyHamming, ncol(data)),
  )
  for (i in seq_len(nrow(data))) {
    annoy$addItem(i - 1, data[i, ])
  }
  annoy$build(n_trees)

  idx <- matrix(nrow = nrow(query), ncol = k)
  dist <- matrix(nrow = nrow(query), ncol = k)
  rownames(idx) <- rownames(query)
  rownames(dist) <- rownames(query)
  for (i in seq_len(nrow(query))) {
    res <- annoy$getNNsByVectorList(query[i, ], k, search_k, include_distances = TRUE)
    idx[i, ] <- res$item + 1
    dist[i, ] <- res$dist
  }
  if (metric == "cosine") dist <- 0.5 * (dist * dist)
  list(idx = idx, dist = dist)
}
