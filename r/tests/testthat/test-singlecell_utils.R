# Copyright 2025 BPCells contributors
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

generate_dense_matrix <- function(nrow, ncol) {
  m <- matrix(runif(nrow * ncol), nrow = nrow)
}

generate_dense_matrix <- function(nrow, ncol) {
  m <- matrix(runif(nrow * ncol), nrow = nrow)
}

test_that("select_features works general case", {
  m1 <- generate_sparse_matrix(100, 50) %>% as("IterableMatrix")
  for (fn in c("select_features_variance", "select_features_dispersion", "select_features_mean")) {
    res <- do.call(fn, list(m1, num_feats = 5))
    expect_equal(nrow(res), nrow(m1)) # Check that dataframe has correct features we're expecting
    expect_equal(sum(res$highly_variable), 5) # Only 10 features marked as highly variable
    expect_setequal(res$feature, rownames(m1))
    res_more_feats_than_rows <- suppressWarnings(do.call(fn, list(m1, num_feats = 10000))) # more features than rows
    res_feats_equal_rows <- do.call(fn, list(m1, num_feats = 100))
    res_feats_partial <- get(fn)(num_feats = 100)(m1)
    expect_identical(res_feats_equal_rows, res_feats_partial)
    expect_identical(res_more_feats_than_rows, res_feats_equal_rows)
    if (fn == "select_features_variance") {
      # Check that normalization actually does something
      res_no_norm <- do.call(fn, list(m1, num_feats = 10, normalize_method = NULL))
      # Check that we can do partial functions on normalization too
      res_norm_partial <- do.call(fn, list(m1, num_feats = 10, normalize_method = normalize_log(scale = 1e3, threads = 1L)))
      res_norm_implicit_partial <- select_features_variance(normalize_method = normalize_log(scale_factor = 1e3), num_feats = 10)(m1)
      expect_identical(res_norm_partial, res_norm_implicit_partial)
      expect_true(!all((res_no_norm %>% dplyr::arrange(feature))$score == (res_norm_partial %>% dplyr::arrange(feature))$score))
    }
  }
})


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
  m1_int <- convert_matrix_type(m1, "uint32_t")
  for (matrices in list(list(m0_t, m1_t), list(m0, m1), list(m0, m1_int))) {
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



test_that("Feature selection by bin variance works", {
    mat <- generate_sparse_matrix(500, 26, fraction_nonzero = 0.1) %>% as("IterableMatrix")
    # Test only that outputs are reasonable.  There is a full comparison in `tests/real_data/` that compares implementation to Seurat
    res_table <- select_features_binned_dispersion(mat, num_feats = 10, n_bins = 5, threads = 1)
    res_table_t <- select_features_binned_dispersion(t(mat), num_feats = 10, n_bins = 5, threads = 1)
    res_feats <- res_table %>% dplyr::filter(highly_variable) %>% dplyr::pull(feature) 
    res <- mat[res_feats,]
    res_feats_t <- res_table_t %>% dplyr::filter(highly_variable) %>% dplyr::pull(feature)
    res_t <- t(mat[,res_feats_t])
    
    expect_equal(nrow(res), 10)
    expect_equal(ncol(res), 26)
    expect_equal(nrow(res_t), 10)
    expect_equal(ncol(res_t), 500)
})

test_that("LSI works", {
    set.seed(12345)
    mat <- matrix(runif(240), nrow=10) %>% as("dgCMatrix") %>% as("IterableMatrix")
    rownames(mat) <- paste0("feat", seq_len(nrow(mat)))
    colnames(mat) <- paste0("cell", seq_len(ncol(mat)))
    # Test only that outputs are reasonable.  There is a full comparison in `tests/real_data/` that compares implementation to ArchR
    n_dimensions <- 5
    lsi_res_obj <- LSI(mat, n_dimensions = n_dimensions)
    lsi_res_t_obj <- LSI(t(mat), n_dimensions = n_dimensions)
    lsi_res <- lsi_res_obj$cell_embeddings
    lsi_res_t <- lsi_res_t_obj$cell_embeddings
    # Check that projection results in the same output if used on the same input matrix
    lsi_res_proj <- project(lsi_res_obj, mat)
    # Check setting pca correlations to non-1 value
    lsi_res_obj_corr <- LSI(mat, n_dimensions = n_dimensions, corr_cutoff = 0.2)

    expect_equal(ncol(lsi_res), 5)
    expect_equal(nrow(lsi_res), ncol(mat))
    expect_equal(ncol(lsi_res_t), 5)
    expect_equal(nrow(lsi_res_t), nrow(mat))
    expect_equal(nrow(lsi_res_proj), ncol(mat))
    expect_lt(ncol(lsi_res_obj_corr$cell_embeddings), n_dimensions)
    expect_equal(lsi_res, lsi_res_proj, tolerance = 1e-7)
})

test_that("Iterative LSI works", {
  skip_if_not_installed("RcppHNSW")
  mat <- matrix(data = runif(50000, 0, 1), nrow=500, ncol = 100) %>% as("dgCMatrix") %>% as("IterableMatrix")
  rownames(mat) <- paste0("feat", seq_len(nrow(mat)))
  colnames(mat) <- paste0("cell", seq_len(ncol(mat)))
  lsi_res_obj <- expect_no_error(IterativeLSI(mat, n_iterations = 2, n_dimensions = 10L, cluster_method = cluster_cells_graph(knn_method = knn_annoy)))
  lsi_res_proj <- project(lsi_res_obj, mat)
  lsi_res_proj_iter_1 <- expect_no_error(project(lsi_res_obj, mat, iteration = 1L))
  lsi_res_embedding <- lsi_res_obj$cell_embeddings
  expect_equal(nrow(lsi_res_embedding), ncol(mat))
  expect_equal(ncol(lsi_res_embedding), 10)
  expect_equal(nrow(lsi_res_proj_iter_1), ncol(mat))
  expect_equal(ncol(lsi_res_proj_iter_1), 10)
  expect_equal(lsi_res_embedding, lsi_res_proj, tolerance = 1e-7)
})

test_that("Iterative LSI works with parameterized clustering", {
  skip_if_not_installed("RcppAnnoy")
  mat <- matrix(data = runif(50000, 0, 1), nrow=500, ncol = 100) %>% as("dgCMatrix") %>% as("IterableMatrix")
  rownames(mat) <- paste0("feat", seq_len(nrow(mat)))
  colnames(mat) <- paste0("cell", seq_len(ncol(mat)))
  lsi_res_obj <- expect_no_error(
    IterativeLSI(
      mat, n_dimensions = 10L,
      cluster_method = cluster_cells_graph(
        knn_method = knn_annoy(k = 12),
        knn_to_graph_method = knn_to_snn_graph(min_val = 0.1)
      )
    )
  )
  lsi_res_proj <- project(lsi_res_obj, mat)
  lsi_res_embedding <- lsi_res_obj$cell_embeddings
  expect_equal(nrow(lsi_res_embedding), ncol(mat))
  expect_equal(ncol(lsi_res_embedding), 10)
  expect_equal(lsi_res_embedding, lsi_res_proj, tolerance = 1e7)
})