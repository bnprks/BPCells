
generate_sparse_matrix <- function(nrow, ncol, fraction_nonzero = 0.5, max_val = 10) {
  m <- matrix(rbinom(nrow * ncol, 1, fraction_nonzero) * sample.int(max_val, nrow * ncol, replace = TRUE), nrow = nrow)
  as(m, "dgCMatrix")
}
generate_dense_matrix <- function(nrow, ncol) {
  m <- matrix(runif(nrow * ncol), nrow = nrow)
}

to_matrix <- function(x) {
  x <- as.matrix(x)
  attr(x, "dimnames") <- NULL
  x
}
to_vector <- function(x) {
  x <- as.numeric(x)
  attributes(x) <- NULL
  x
}


test_that("MatrixStats basic test", {
  withr::local_seed(195123)

  m1 <- generate_sparse_matrix(5, 1000)
  m2 <- t(m1)

  i1 <- as(m1, "IterableMatrix")
  i2 <- t(i1)


  stats1 <- matrix_stats(i1, "variance", "variance")
  stats2 <- matrix_stats(i2, "variance", "variance")

  expect_identical(c("nonzero", "mean", "variance"), rownames(stats1$row_stats))
  expect_identical(c("nonzero", "mean", "variance"), rownames(stats1$col_stats))

  expect_identical(c("nonzero", "mean", "variance"), rownames(stats2$row_stats))
  expect_identical(c("nonzero", "mean", "variance"), rownames(stats2$col_stats))

  expect_equal(stats1$row_stats["nonzero", ], rowSums(as.matrix(m1 != 0)), tolerance = testthat_tolerance())
  expect_equal(stats1$row_stats["mean", ], rowMeans(as.matrix(m1)), tolerance = testthat_tolerance())

  expect_equal(stats1$row_stats["variance", ], matrixStats::rowVars(as.matrix(m1)), tolerance = testthat_tolerance())
  expect_equal(stats1$col_stats["nonzero", ], colSums(as.matrix(m1 != 0)), tolerance = testthat_tolerance())

  expect_equal(stats1$col_stats["mean", ], colMeans(as.matrix(m1)), tolerance = testthat_tolerance())
  expect_equal(stats1$col_stats["variance", ], matrixStats::colVars(as.matrix(m1)), tolerance = testthat_tolerance())
})

# Given a dgCMatrix m1 and its equivalent IterableMatrix i1, check
# that all permutations of matrix_stats work as expeected
test_stats_comprehensive <- function(m1, i1, threads=0L) {
  m2 <- t(m1)
  i2 <- t(i1)
  for (row_stats in c("none", "nonzero", "mean", "variance")) {
    for (col_stats in c("none", "nonzero", "mean", "variance")) {
      stats1 <- matrix_stats(i1, row_stats, col_stats, threads=threads)
      stats2 <- matrix_stats(i2, row_stats, col_stats, threads=threads)

      expected_rownames <- c("nonzero", "mean", "variance")[seq_len(match(row_stats, c("none", "nonzero", "mean", "variance")) - 1)]
      expected_colnames <- c("nonzero", "mean", "variance")[seq_len(match(col_stats, c("none", "nonzero", "mean", "variance")) - 1)]

      if (row_stats == "none") expected_rownames <- NULL
      if (col_stats == "none") expected_colnames <- NULL

      expect_identical(expected_rownames, rownames(stats1$row_stats))
      expect_identical(expected_colnames, rownames(stats1$col_stats))

      expect_identical(expected_rownames, rownames(stats2$row_stats))
      expect_identical(expected_colnames, rownames(stats2$col_stats))

      if (row_stats %in% c("nonzero", "mean", "variance")) {
        expect_equal(stats1$row_stats["nonzero", ], rowSums(as.matrix(m1 != 0)), tolerance = testthat_tolerance())
      }
      if (row_stats %in% c("mean", "variance")) {
        expect_equal(stats1$row_stats["mean", ], rowMeans(as.matrix(m1)), tolerance = testthat_tolerance())
      }
      if (row_stats %in% c("variance")) {
        expect_equal(stats1$row_stats["variance", ], matrixStats::rowVars(as.matrix(m1)), tolerance = testthat_tolerance())
      }

      if (col_stats %in% c("nonzero", "mean", "variance")) {
        expect_equal(stats1$col_stats["nonzero", ], colSums(as.matrix(m1 != 0)), tolerance = testthat_tolerance())
      }
      if (col_stats %in% c("mean", "variance")) {
        expect_equal(stats1$col_stats["mean", ], colMeans(as.matrix(m1)), tolerance = testthat_tolerance())
      }
      if (col_stats %in% c("variance")) {
        expect_equal(stats1$col_stats["variance", ], matrixStats::colVars(as.matrix(m1)), tolerance = testthat_tolerance())
      }
    }
  }
}

test_that("MatrixStats comprehensive tests", {
  withr::local_seed(195123)
  m1 <- generate_sparse_matrix(5, 1000)
  i1 <- as(m1, "IterableMatrix")

  test_stats_comprehensive(m1, i1)
})


test_that("rbind MatrixStats comprehensive tests", {
  withr::local_seed(195123)

  m1 <- generate_sparse_matrix(1000, 10)

  mat_list <- lapply(1:10, function(i) {
    as(m1, "IterableMatrix")[(i-1)*100 + 1:100,]
  })

  i1 <- do.call(rbind, mat_list)
  test_stats_comprehensive(m1, i1)
})

test_that("MatrixStats multithreaded works", {
  withr::local_seed(195123)
  m1 <- generate_sparse_matrix(5, 1000)
  i1 <- as(m1, "IterableMatrix")

  test_stats_comprehensive(m1, i1, threads=2)
})

test_that("svds works", {
  withr::local_seed(195123)
  m1 <- matrix(runif(240), nrow=10)
  m1_t <- t(m1)
  i1 <- as(m1, "dgCMatrix") %>% as("IterableMatrix")
  i1_t <- as(m1_t, "dgCMatrix") %>% as("IterableMatrix")

  ans <- svd(m1)
  ans_t <- svd(m1_t)
  mine <- svds(i1, k=5)
  mine_t <- svds(i1_t, k=5)

  mine_threaded <- svds(i1, k=5, threads=2)
  mine_t_threaded <- svds(i1_t, k=5, threads=2)

  equal_svds <- function(ans, mine) {
    k <- length(mine$d)
    sign_flip <- sign(ans$u[1,1:k]) * sign(mine$u[1,])

    expect_equal(ans$d[1:k], mine$d)
    expect_equal(ans$u[,1:k], multiply_cols(mine$u, sign_flip))
    expect_equal(ans$v[,1:k], multiply_cols(mine$v, sign_flip))
  }
  equal_svds(ans, mine)
  equal_svds(ans_t, mine_t)

  equal_svds(ans, mine_threaded)
  equal_svds(ans_t, mine_t_threaded)

  equal_svds(ans_t, svds(t(i1), k=5))
})