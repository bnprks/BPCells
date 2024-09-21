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
  for (matrices in list(list(m0, m1), list(m0_t, m1_t))) {
    # Check across two equal groups, one group, and groups of equal length
    m <- matrices[[1]]
    m_bpcells <- matrices[[2]]
    groups <- c(rep.int(1, ncol(m)/2), rep.int(2, ncol(m)/2))
    groups_one_type <- c(rep.int(1, ncol(m)))
    groups_equal_length <- seq(ncol(m))
    for (cell_group in list(groups, groups_one_type, groups_equal_length)) {
      # create clipping
      min_vals  <- m %>% apply(2, quantile, 0.99, type = 1)
      #pmin doesn't work correctly when comparing matrix to vec, so we do it iteratively
      m_sum <- create_pseudobulk_r(m, cell_group, "sum")
      m_mean <- create_pseudobulk_r(m, cell_group, "mean")
      m_non_zeros <- create_pseudobulk_r(m, cell_group, "non_zeros")
      m_bpcells_sum <- pseudobulk_matrix(m_bpcells, cell_group, method = "sum")$sum
      m_bpcells_mean <- pseudobulk_matrix(m_bpcells, cell_group, method = "mean")$mean
      m_bpcells_non_zeros <- pseudobulk_matrix(m_bpcells, cell_group, method = "non-zeros")$non_zeros
      expect_equal(m_sum, m_bpcells_sum)
      # check for clipping if matrix isn't transposed
      if (m_bpcells@transpose == FALSE) {
        m_clipped <- m
        m_sum_clipped <- create_pseudobulk_r(m, cell_group, "sum")
        m_bpcells_clipped <- pseudobulk_matrix(m_bpcells, cell_group, method = "sum", clip_values = TRUE)$sum
        expect_equal(m_sum_clipped, m_bpcells_clipped)
      }
      expect_equal(m_mean, m_bpcells_mean)
      expect_equal(m_non_zeros, m_bpcells_non_zeros)
      # make sure that we dont check for variances if we have number of groups == number of cells
      if (length(unique(cell_group)) < ncol(m)) {
        m_clipped <- m
        for (idx in seq_along(min_vals)) {
          m_clipped[, idx] <- pmin(m[, idx], min_vals[[idx]])
        }
        m_var <- create_pseudobulk_r(m, cell_group, "var")
        m_bpcells_var <- pseudobulk_matrix(m_bpcells, cell_group, method = "var")
        expect_equal(m_var, m_bpcells_var$var)
      }
    }
  }
})

test_that("Matrix quantiles work", {
  m0 <- generate_sparse_matrix(20, 10)
  m1 <- m0 |> as("IterableMatrix")
  m0 <- as.matrix(m0)
  for (quantile in c(0, 0.25, 0.5, 0.75, 0.99)) {
    m_quantile <- tibble::as_tibble(m0) %>% apply(2, quantile, probs = quantile, type = 1)
    m_quantile_bpcells <- matrix_quantile_per_cell(m1, quantile)
    expect_equal(m_quantile, m_quantile_bpcells)
  }
})