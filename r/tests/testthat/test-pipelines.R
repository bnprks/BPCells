# Copyright 2024 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

generate_sparse_matrix <- function(nrow, ncol, fraction_nonzero = 0.5, max_val = 10) {
  m <- matrix(rbinom(nrow * ncol, 1, fraction_nonzero) * sample.int(max_val, nrow * ncol, replace = TRUE), nrow = nrow)
  rownames(m) <- paste0("feat", seq_len(nrow(m)))
  colnames(m) <- paste0("cell", seq_len(ncol(m)))
  as(m, "dgCMatrix")
}

test_that("Highly variable feature pipeline works", {
  mat <- generate_sparse_matrix(500, 26, fraction_nonzero = 0.1) %>% as("IterableMatrix")
  rownames(mat) <- paste0("feat", seq_len(nrow(mat)))
  colnames(mat) <- paste0("cell", seq_len(ncol(mat)))
  # Test only that outputs are reasonable.  There is a full comparison in `tests/real_data/` that compares implementation to Seurat
  hvf_transformer <- VarFeatSelectorTransformer(num_feats = 10, n_bins = 5)
  res <- fit(hvf_transformer, mat) %>% project(mat)
  expect_equal(nrow(res), 10)
  expect_equal(ncol(res), 26)
})

test_that("LSI Pipeline works", {
  mat <- matrix(runif(240), nrow=10) %>% as("dgCMatrix") %>% as("IterableMatrix")
  rownames(mat) <- paste0("feat", seq_len(nrow(mat)))
  colnames(mat) <- paste0("cell", seq_len(ncol(mat)))
  lsi_transformer <- LSITransformer(n_dimensions = 5)
  res <- fit(lsi_transformer, mat) %>% project(mat)
  expect_equal(nrow(res), 5)
  expect_equal(ncol(res), ncol(mat))
})

test_that("Pipeline with var feat selection and LSI works", {
  mat <- generate_sparse_matrix(500, 26, fraction_nonzero = 0.1) %>% as("IterableMatrix")
  rownames(mat) <- paste0("feat", seq_len(nrow(mat)))
  colnames(mat) <- paste0("cell", seq_len(ncol(mat)))
  # Test only that outputs are reasonable.  There is a full comparison in `tests/real_data/` that compares implementation to Seurat
  hvf_transformer <- VarFeatSelectorTransformer(num_feats = 25, n_bins = 5)
  lsi_transformer <- LSITransformer(n_dimensions = 5)
  # test pipeline creation with both `c()` and `Pipeline()` constructor
  pipeline_1 <- c(hvf_transformer, lsi_transformer)
  pipeline_2 <- Pipeline(
    steps = list(hvf_transformer, lsi_transformer)
  )
  assert_true(all.equal(pipeline_1, pipeline_2))
  res <- fit(pipeline_1, mat) %>% project(mat)
  expect_equal(nrow(res), 5)
  expect_equal(ncol(res), 26)
})