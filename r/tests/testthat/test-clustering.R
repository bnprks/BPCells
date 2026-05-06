# Copyright 2025 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

test_that("C++ SNN calculation works",{
    
    k <- 10

    for (cells in c(100, 1000)) {
        neighbor_sim <- matrix(sample.int(cells, cells*k, replace=TRUE), nrow=cells)

        # Make sure no cell is listed twice as a nearest neighbor
        for (i in seq_len(cells)) {
            while (length(unique(neighbor_sim[i,])) != k) {
                neighbor_sim[i,] <- sample.int(cells, k)
            }
        }
        input <- list(idx = neighbor_sim, dist = matrix(runif(cells*k), nrow=cells))
        min_val <- 1/15
        snn <- knn_to_snn_graph(input, min_val=min_val, self_loops=TRUE)
        mat <- knn_to_graph(input, use_weights=FALSE)

        mat <- mat %*% t(mat)
        mat <- mat / (2 * k - mat)
        mat@x[mat@x < min_val] <- 0
        # Prune the explicit 0 entries from storage
        mat <- Matrix::drop0(mat) 
        mat <- Matrix::tril(mat)
        expect_identical(
            snn,
            as(mat, "dgCMatrix")
        )
        snn2 <- knn_to_snn_graph(input, min_val=min_val, self_loops=FALSE)
        diag(mat) <- 0
        mat <- Matrix::drop0(mat)
        expect_identical(
            snn2,
            as(mat, "dgCMatrix")
        )
    }

})

test_that("C++ UMAP graph calculation works", {
    # This uses a pre-calculated graph from the umap-learn python package.
    # See ../data/generate_iris_geodesic_graph.R for the generation code 
    test_data <- readRDS("../data/iris_geodesic_graph.rds")
    knn <- test_data$knn
    res <- knn_to_geodesic_graph(knn)
    expect_equal(as.matrix(res + t(res)), as.matrix(test_data$graph), tolerance=1e-6)
})

test_that("igraph clustering doesn't crash", {
    skip_if_not_installed("igraph")
    test_data <- readRDS("../data/iris_geodesic_graph.rds")
    knn <- test_data$knn
    graph <- knn_to_geodesic_graph(knn)
    graph_list <- knn_to_geodesic_graph(knn, return_type="list")

    # The `resolution_parameter` param in igraph `cluster_leiden()` is deprecated,
    # causing `expect_no_condition()` to fail. This workaround avoids test failures from 
    # the deprecation warning, but is confirmed to fail if a call to `warning()` or `stop()` is
    # inserted into `cluster_graph_leiden()`
    suppressWarnings({
        expect_no_warning(cluster_graph_leiden(graph))
        expect_no_warning(cluster_graph_leiden(graph, objective_function="CPM"))
        expect_no_error(cluster_graph_leiden(graph))
        expect_no_error(cluster_graph_leiden(graph, objective_function="CPM"))
    })
    expect_identical(
        suppressWarnings(cluster_graph_leiden(graph)),
        suppressWarnings(cluster_graph_leiden(graph_list))
    )

    expect_no_condition(cluster_graph_louvain(graph))
    expect_no_condition(cluster_graph_louvain(graph_list))
    expect_identical(cluster_graph_louvain(graph), cluster_graph_louvain(graph_list))
})

test_that("knn_hnsw rownames come from query", {
    skip_if_not_installed("RcppHNSW")
    data <- matrix(rnorm(200), nrow=20, dimnames=list(paste0("data_", 1:20)))
    query <- matrix(rnorm(50), nrow=5, dimnames=list(paste0("query_", 1:5)))
    res <- knn_hnsw(data, query=query, k=3, verbose=FALSE)
    expect_equal(rownames(res$idx), rownames(query))
    expect_equal(rownames(res$dist), rownames(query))
})

test_that("cluster_cells_graph works", {
    skip_if_not_installed("RcppAnnoy")
    skip_if_not_installed("RcppHNSW")
    mat <- matrix(sample.int(1000, 10000, replace=TRUE), nrow=1000)
    # check with default params
    res <- expect_no_error(cluster_cells_graph(mat))
    # check with threads, method partialization
    expect_true(class(res) == "factor")
    expect_equal(nrow(mat), length(res))
    res_partialized <- expect_no_error(
        cluster_cells_graph(
            mat, knn_method = knn_annoy(k = 9),
            knn_to_graph_method = knn_to_snn_graph(min_val = 1/10),
            cluster_graph_method = cluster_graph_louvain(resolution = 0.8),
        ))
    expect_true(class(res) == "factor")
    expect_equal(nrow(mat), length(res))
})
