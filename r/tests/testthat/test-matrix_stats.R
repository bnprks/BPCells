# Copyright 2021 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.


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
  skip_if_not_installed("matrixStats")
  withr::local_seed(195123)

  m1 <- generate_sparse_matrix(5, 1000)
  m2 <- t(m1)

  i1 <- as(m1, "IterableMatrix")
  i2 <- t(i1)

  
  stats_m1 <- matrix_stats(as.matrix(m1), "variance", "variance")
  stats_m1_dgc <- matrix_stats(m1, "variance", "variance")
  stats1 <- matrix_stats(i1, "variance", "variance")
  stats2 <- matrix_stats(i2, "variance", "variance")

  expect_error(matrix_stats(list("a"), "variance"), "cannot be converted to an IterableMatrix object")

  expect_identical(stats_m1_dgc, stats_m1)
  expect_identical(stats_m1_dgc, stats1)

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
  skip_if_not_installed("matrixStats")
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

test_that("MatrixStats with row/col shift comprehensive tests", {
  withr::local_seed(195123)
  m1 <- generate_sparse_matrix(100, 50)
  i1 <- as(m1, "IterableMatrix")

  test_stats_comprehensive(add_cols(m1, seq_len(ncol(m1))), add_cols(i1, seq_len(ncol(m1))))
  test_stats_comprehensive(add_rows(m1, seq_len(nrow(m1))), add_rows(i1, seq_len(nrow(m1))))
})

test_that("MatrixStats multithreaded works", {
  withr::local_seed(195123)
  m1 <- generate_sparse_matrix(5, 1000)
  i1 <- as(m1, "IterableMatrix")

  test_stats_comprehensive(m1, i1, threads=2)
})

equal_svds <- function(ans, mine) {
  k <- length(mine$d)
  sign_flip <- sign(ans$u[1,1:k]) * sign(mine$u[1,])

  expect_equal(ans$d[1:k], mine$d)
  expect_equal(ans$u[,1:k], multiply_cols(mine$u, sign_flip))
  expect_equal(ans$v[,1:k], multiply_cols(mine$v, sign_flip))
}

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

  equal_svds(ans, mine)
  equal_svds(ans_t, mine_t)

  equal_svds(ans, mine_threaded)
  equal_svds(ans_t, mine_t_threaded)

  equal_svds(ans_t, svds(t(i1), k=5))
})

test_that("svds registers generic with RSpectra", {
  skip_if_not_installed("RSpectra")

  withr::local_seed(195123)
  m1 <- matrix(runif(240), nrow=10)
  i1 <- as(m1, "dgCMatrix") %>% as("IterableMatrix")

  ans <- svd(m1)
  
  # Ensure RSpectra works with BPCells and BPCells works with RSpectra
  equal_svds(ans, RSpectra::svds(i1, k=5))
  equal_svds(ans, BPCells::svds(m1, k=5))
})

test_that("rowMaxs and colMaxs works comprehensive", {
  skip_if_not_installed("matrixStats")
  withr::local_seed(195123)
  m1_neg <- matrix(runif(12, min = -10, max = -1), nrow = 3, ncol = 4)
  m2_dense <- generate_dense_matrix(4, 64)
  # Make a sparse matrix with half positive and half negative values
  m3_sparse <- generate_sparse_matrix(4, 512)
  m3_sparse@x <- m3_sparse@x * sample(c(-1,1), length(m3_sparse@x), replace=TRUE)

  m_test_cases <- list(m1_neg, m2_dense, m3_sparse)
  for (m in m_test_cases) {
    i_transpose <- t(m) %>% as("dgCMatrix") %>% as("IterableMatrix")
    i <- m %>% as("dgCMatrix") %>% as("IterableMatrix")
    m_tranpose_cases <- list(i_transpose, t(i_transpose), i, t(i))
    for (i in m_tranpose_cases) {
      expect_identical(rowMaxs(i), matrixStats::rowMaxs(as.matrix(i)))
      expect_identical(colMaxs(i), matrixStats::colMaxs(as.matrix(i)))
    }
  }
})

test_that("Matrix col quantiles works", {
  skip_if_not_installed("matrixStats")
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
  skip_if_not_installed("matrixStats")
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