generate_sparse_matrix <- function(nrow, ncol, fraction_nonzero = 0.5, max_val = 10) {
  m <- matrix(rbinom(nrow * ncol, 1, fraction_nonzero) * sample.int(max_val, nrow * ncol, replace = TRUE), nrow = nrow)
  rownames(m) <- paste0("row", seq_len(nrow(m)))
  colnames(m) <- paste0("col", seq_len(ncol(m)))
  as(m, "dgCMatrix")
}
generate_dense_matrix <- function(nrow, ncol) {
  matrix(runif(nrow * ncol), nrow = nrow)
}
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

test_that("Matrix col quantiles works", {
  m0 <- generate_sparse_matrix(20, 10)
  m1 <- m0 |> as("IterableMatrix")
  m0 <- as.matrix(m0)
  # For single quantiles
  for (quantile in c(0, 0.25, 0.5, 0.75, 0.99)) {
    for (type in c(4, 5, 6, 7, 8, 9)) {
      m_quantile <- m0 %>% matrixStats::colQuantiles(probs = quantile, type = type)
      m_quantile_bpcells <- colQuantiles(m1, probs = quantile, type = type)
      expect_equal(m_quantile, m_quantile_bpcells, info = paste("Quantile:", quantile, "Type:", type))
    }
  }
  # for multiple quantiles
  for (quantiles in list(c(0, 0.25, 0.5, 0.75, 0.99), c(0.25, 0.5, 0.75))) {
    for (type in c(4, 5, 6, 7, 8, 9)) {
      m_quantile <- m0 %>% matrixStats::colQuantiles(probs = quantiles, type = type)
      m_quantile_bpcells <- colQuantiles(m1, probs = quantiles, type = type)
      expect_equal(m_quantile, m_quantile_bpcells, info = paste("Quantiles:", quantiles, "Type:", type))
    }
  }
})

test_that("Matrix row quantiles works with transposed matrices", {
  m0 <- generate_sparse_matrix(20, 10)
  m1 <- m0 |> as("IterableMatrix") |> t()
  m0 <- as.matrix(m0) |> t()
  # For single quantile
  m_quantile <-  m0 %>% matrixStats::rowQuantiles(probs = 0.75, type = 7)
  m_quantile_bpcells <- rowQuantiles(m1, probs = 0.75, type = 7)
  expect_equal(m_quantile, m_quantile_bpcells)
  # With multiple quantiles
  m_quantiles <-  m0 %>% matrixStats::rowQuantiles(probs = c(0.25, 0.75), type = 7)
  m_quantiles_bpcells <- rowQuantiles(m1, probs = c(0.25, 0.75), type = 7)
  expect_equal(m_quantiles, m_quantiles_bpcells)
})