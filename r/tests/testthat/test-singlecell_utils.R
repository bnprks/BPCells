# Copyright 2023 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

generate_sparse_matrix <- function(nrow, ncol, fraction_nonzero = 0.5, max_val = 10) {
  m <- matrix(rbinom(nrow * ncol, 1, fraction_nonzero) * sample.int(max_val, nrow * ncol, replace = TRUE), nrow = nrow)
  rownames(m) <- paste0("row", seq_len(nrow(m)))
  colnames(m) <- paste0("col", seq_len(ncol(m)))
  as(m, "dgCMatrix")
}

test_that("Wilcoxon rank sum works (small)", {
    x <- c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46)
    y <- c(1.15, 0.88, 0.90, 0.74, 1.21)

    med_xy <- median(c(x,y))
    x <- x - med_xy
    y <- y - med_xy

    res_r <- wilcox.test(x, y, exact=FALSE, correct=TRUE)

    vals <- c(x,y)
    X <- matrix(c(vals - median(vals), vals), ncol=2)
    colnames(X) <- c("centered", "uncentered")
    Y <- factor(c(rep.int("X", length(x)), rep.int("Y", length(y))))

    m <- as(X, "dgCMatrix") %>% 
        as("IterableMatrix") %>%
        t()

    res_bp_1 <- marker_features(m, Y) %>%
        dplyr::filter(foreground == "X", feature == "centered")
    res_bp_2 <- marker_features(m, Y) %>%
        dplyr::filter(foreground == "Y", feature == "uncentered")
    
    expect_equal(res_bp_1$p_val_raw, res_r$p.value)
    expect_equal(res_bp_1$p_val_raw, res_bp_2$p_val_raw)

    expect_equal(res_bp_1$foreground_mean, mean(x))
    expect_equal(res_bp_1$background_mean, mean(y))

    expect_equal(res_bp_2$foreground_mean, mean(y))
    expect_equal(res_bp_2$background_mean, mean(x))
    
    expect_equal(res_bp_1$feature, "centered")
})

test_that("Wilcoxon rank sum matches immunogenomics::presto", {
    skip_if_not_installed("pbmc3k.SeuratData")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("presto")
    mat <- pbmc3k.SeuratData::pbmc3k.final@assays$RNA@data
    clusters <- pbmc3k.SeuratData::pbmc3k.final@active.ident
    res_presto <- presto::wilcoxauc(mat, clusters) %>%
        tibble::as_tibble()
    
    mat_bpcells <- mat %>% 
        t() %>% 
        as("IterableMatrix") %>% 
        t()
    res_bpcells <- marker_features(mat_bpcells, clusters) %>%
        dplyr::arrange(match(foreground, levels(clusters)), match(feature, rownames(mat)))
    
    expect_equal(res_bpcells$feature, res_presto$feature)
    expect_equal(res_bpcells$foreground, res_presto$group)
    expect_equal(res_bpcells$p_val_raw, res_presto$pval)

    expect_equal(res_bpcells$foreground_mean, res_presto$avgExpr)
    expect_equal(res_bpcells$background_mean, res_presto$avgExpr - res_presto$logFC)

})

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

        min_val <- 1/15
        snn <- knn_to_snn_graph(list(idx=neighbor_sim), min_val=min_val, self_loops=TRUE)

        mat <- knn_to_graph(list(idx=neighbor_sim), use_weights=FALSE)

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
        snn2 <- knn_to_snn_graph(list(idx=neighbor_sim), min_val=min_val, self_loops=FALSE)
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

non_zeros <- function(vec) {
  sum(vec > 0)
}
create_pseudobulk_r <- function(mat, cell_group, method) {
  res <- suppressWarnings(
    mat %>%
    t() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(group = cell_group) %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), !!as.symbol(method))) %>%
    tibble::column_to_rownames("group") %>%
    t()
  )
  return(as.matrix(res))
}

test_that("Pseudobulk aggregation works", {
  m0 <- generate_sparse_matrix(20, 10, max_val = 10)
  m1 <- m0 |> as("IterableMatrix")
  m0 <- as.matrix(m0)
  m0_t <- t(m0)
  m1_t <- t(m1)
  for (matrices in list(list(m0_t, m1_t), list(m0, m1))) {
    # Check across two equal groups, one group, and groups of equal length
    m <- matrices[[1]]
    m_bpcells <- matrices[[2]]
    groups <- c(rep.int(1, ncol(m)/2), rep.int(2, ncol(m)/2))
    groups_one_type <- c(rep.int(1, ncol(m)))
    groups_equal_length <- seq(ncol(m))
    for (cell_group in list(groups, groups_one_type, groups_equal_length)) {
      m_sum <- create_pseudobulk_r(m, cell_group, "sum")
      m_mean <- create_pseudobulk_r(m, cell_group, "mean")
      m_non_zeros <- create_pseudobulk_r(m, cell_group, "non_zeros")
      m_bpcells_sum <- pseudobulk_matrix(m_bpcells, cell_group, method = "sum")
      m_bpcells_mean <- pseudobulk_matrix(m_bpcells, cell_group, method = "mean")
      m_bpcells_non_zeros <- pseudobulk_matrix(m_bpcells, cell_group, method = "nonzeros")
      expect_equal(m_sum, m_bpcells_sum)
      # check for clipping if matrix isn't transposed
      expect_equal(m_mean, m_bpcells_mean)
      expect_equal(m_non_zeros, m_bpcells_non_zeros)
      # make sure that we dont check for variances if we have number of groups == number of cells
      if (length(unique(cell_group)) < ncol(m)) {
        m_var <- create_pseudobulk_r(m, cell_group, "var")
        m_bpcells_var <- pseudobulk_matrix(m_bpcells, cell_group, method = "variance")
        expect_equal(m_var, m_bpcells_var, info = paste("Transposed:", m_bpcells@transpose))
      }
    }
  }
})

test_that("Pseudobulk aggregation works with multiple return types", {
  m0 <- generate_sparse_matrix(20, 10, max_val = 10)
  m1 <- m0 |> as("IterableMatrix")
  m0 <- as.matrix(m0)
  m0_t <- t(m0)
  m1_t <- t(m1)
  methods <- c("nonzeros", "sum", "mean", "variance")
  for (matrices in list(list(m0_t, m1_t), list(m0, m1))) {
    # Check across two equal groups, one group, and groups of equal length
    m <- matrices[[1]]
    m_bpcells <- matrices[[2]]
    cell_group <- c(rep.int(1, ncol(m)/2), rep.int(2, ncol(m)/2))
    m_sum <- create_pseudobulk_r(m, cell_group, "sum")
    m_mean <- create_pseudobulk_r(m, cell_group, "mean")
    m_non_zeros <- create_pseudobulk_r(m, cell_group, "non_zeros")
    m_var <- create_pseudobulk_r(m, cell_group, "var")
    for (first_method_idx in seq(1, length(methods)-1)) {
      for (second_method_idx in seq(first_method_idx+1, length(methods))) {
        methods_of_interest <- c(methods[first_method_idx], methods[second_method_idx])
        m_bpcells_res <- pseudobulk_matrix(m_bpcells, cell_group, method = methods_of_interest)
        if ("variance" %in% methods_of_interest) {
          expect_equal(m_var, m_bpcells_res$var)
        }
        if ("mean" %in% methods_of_interest) {
          expect_equal(m_mean, m_bpcells_res$mean)
        }
        if ("sum" %in% methods_of_interest) {
          expect_equal(m_sum, m_bpcells_res$sum)
        }
        if ("nonzeros" %in% methods_of_interest) {
          expect_equal(m_non_zeros, m_bpcells_res$nonzeros)
        }
      }
    }
  }
})