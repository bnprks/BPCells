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
  # Test with mean and sum
  m0 <- generate_sparse_matrix(20, 10, max_val = 10)
  m1 <- m0 |> as("IterableMatrix")
  groups <- c(rep.int(1, 5), rep.int(2, 5))
  m0 <- as.matrix(m0)
  #m0 <- tibble::as_tibble(m0)
  m0_sum <- suppressWarnings(m0 %>% t() %>% tibble::as_tibble() %>% dplyr::mutate(group = groups) %>% dplyr::group_by(group) %>% dplyr::summarise_each(dplyr::funs(sum)) %>% t())
  colnames(m0_sum) <- m0_sum[1,]
  m0_sum <- m0_sum[-1,]
  browser()
}
)