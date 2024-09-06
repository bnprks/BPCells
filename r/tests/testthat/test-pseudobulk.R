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
  # Check across two equal groups, one group, and groups of equal length
  groups <- c(rep.int(1, 5), rep.int(2, 5))
  groups_one_type <- c(rep.int(1, 10))
  group_equal_length <- seq(10)
  m0 <- as.matrix(m0)
  #m0 <- tibble::as_tibble(m0)
  
  for (cell_group in list(groups, groups_one_type, group_equal_length)) {
    # Test with mean and sum
    m0_sum <- suppressWarnings(m0 %>% t() %>% tibble::as_tibble() %>% dplyr::mutate(group = cell_group) %>% dplyr::group_by(group) %>% dplyr::summarise_each(dplyr::funs(sum)) %>% t())
    m0_mean <- suppressWarnings(m0 %>% t() %>% tibble::as_tibble() %>% dplyr::mutate(group = cell_group) %>% dplyr::group_by(group) %>% dplyr::summarise_each(dplyr::funs(mean)) %>% t())
    m0_sum <- as.data.frame(m0_sum)
    m0_mean <- as.data.frame(m0_mean)
    # using this aggregation results in a row indicating group, so we remove it
    colnames(m0_sum) <- m0_sum[1,]
    colnames(m0_mean) <- m0_mean[1,]
    m0_sum <- tail(m0_sum, -1)
    m0_mean <- tail(m0_mean, -1)
    m1_sum <- pseudobulk_counts(m1, cell_group, approach = "sum")
    m1_mean <- pseudobulk_counts(m1, cell_group, approach = "mean")
    expect_equal(m0_sum, m1_sum)
    expect_equal(m0_mean, m1_mean)
  }
})