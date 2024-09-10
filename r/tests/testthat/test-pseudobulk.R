generate_sparse_matrix <- function(nrow, ncol, fraction_nonzero = 0.5, max_val = 10) {
  m <- matrix(rbinom(nrow * ncol, 1, fraction_nonzero) * sample.int(max_val, nrow * ncol, replace = TRUE), nrow = nrow)
  rownames(m) <- paste0("row", seq_len(nrow(m)))
  colnames(m) <- paste0("col", seq_len(ncol(m)))
  as(m, "dgCMatrix")
}
generate_dense_matrix <- function(nrow, ncol) {
  matrix(runif(nrow * ncol), nrow = nrow)
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
      # Test with mean and sum
      m_sum <- suppressWarnings(m %>% t() %>% tibble::as_tibble() %>% dplyr::mutate(group = cell_group) %>% dplyr::group_by(group) %>% dplyr::summarise_each(dplyr::funs(sum)) %>% t())
      m_mean <- suppressWarnings(m %>% t() %>% tibble::as_tibble() %>% dplyr::mutate(group = cell_group) %>% dplyr::group_by(group) %>% dplyr::summarise_each(dplyr::funs(mean)) %>% t())
      m_sum <- as.data.frame(m_sum)
      m_mean <- as.data.frame(m_mean)
      # using this aggregation results in a row indicating group, so we remove it
      colnames(m_sum) <- m_sum[1,]
      colnames(m_mean) <- m_mean[1,]
      m_sum <- tail(m_sum, -1)
      m_mean <- tail(m_mean, -1)
      m_bpcells_sum <- pseudobulk_counts_matrix_multiply(m_bpcells, cell_group, method = "sum")
      m_bpcells_mean <- pseudobulk_counts_matrix_multiply(m_bpcells, cell_group, method = "mean")
      expect_equal(m_sum, m_bpcells_sum)
      expect_equal(m_mean, m_bpcells_mean)
    }
  }
})

test_that("Matrix percentiles work", {
  m0 <- generate_sparse_matrix(20, 10)
  m1 <- m0 |> as("IterableMatrix")
  m0 <- as.matrix(m0)
  m0_t <- t(m0)
  m1_t <- t(m1)
  for (matrices in list(list(m0, m1), list(m0_t, m1_t))) {
    m <- matrices[[1]]
    m_bpcells <- matrices[[2]]
    for (percentile in c(0, 0.25, 0.5, 0.75, 0.99)) {
      m_percentile <- tibble::as_tibble(m) %>% apply(2, quantile, probs = percentile, type = 1)
      m_percentile_bpcells <- matrix_percentile_per_cell(m_bpcells, percentile)
      expect_equal(m_percentile, m_percentile_bpcells)
    }
  }
})