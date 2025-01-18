# Copyright 2024 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

library("BPCells")
library("ArchR")

# Set up temp dir in case it's not already set
create_temp_dir <- function(dir = NULL) {
  if (is.null(dir)) {
    dir <- file.path(tempdir(), "lsi_test")
    if (dir.exists(dir)) unlink(dir, recursive = TRUE)
    dir.create(dir)
  }
  return(dir)
}

#' Perform a dimensionality reduction with tf-idf and SVD (LSI) on a matrix on ArchR and BPCells.
#' As LSI uses an iterative approach on ArchR, we compare by using a single-iteration private function on ArchR.
#' As the SVD implementation is not necessarily the same between the two packages, we take the SVD matrix
#' from both functions and compare the matrix multiplication of the U and SVD matrices, which should give an approximation
#' we can compare between the two packages.
#' @param proj An archr project.
test_lsi_similarity_to_archr <- function(dir = NULL) {
    dir <- create_temp_dir(dir)
    setwd(dir)
    # add iterative lsi for dim reduction
    proj <- getTestProject()
    proj <- addPeakMatrix(proj)
    # Get the peak matrix
    test_mat <- assay(getMatrixFromProject(proj, useMatrix = "PeakMatrix"))
    # Calculate LSI on ArchR
    # running LSI without binarizing, as we don't do this in the BPCells implementation
    # we also don't filter quantile outliers.
    lsi_archr <- ArchR:::.computeLSI(
        mat = test_mat,
        LSIMethod = 2,
        nDimensions = 20,
        binarize = FALSE,
        outlierQuantiles = NULL
    )
    svd_archr <- lsi_archr$svd
    lsi_mat_archr <- t(lsi_archr$matSVD)
    # set rownames to NA, as we don't have rownames in the BPCells implementation
    rownames(lsi_mat_archr) <- NULL
    # PCA Matrix = T(u) * Pre-SVD Matrix
    # u * PCA Matrix = u * T(u) * Pre-SVD Matrix
    # u * PCA Matrix = Pre-SVD Matrix
    pre_svd_mat_approx_archr <- lsi_archr$svd$u %*% lsi_mat_archr
    # Calculate LSI on BPCells
    # Do not use z-score normalization, as this isn't done with ArchR
    lsi_bpcells <- LSI(
        test_mat %>% as("dgCMatrix") %>% as("IterableMatrix"),
        n_dimensions = 20
    )
    pre_svd_mat_approx_bpcells <- lsi_bpcells$fitted_params$svd_params$u %*% lsi_bpcells$cell_embeddings
    testthat::expect_true(all.equal(pre_svd_mat_approx_archr, pre_svd_mat_approx_bpcells, tolerance = 1e-4))
    # convert signs
    lsi_mat_archr <- sweep(lsi_mat_archr, MARGIN = 1,  (2 * (lsi_mat_archr[,1] * lsi_bpcells$cell_embeddings[,1] > 0) - 1), `*`)
    # Check for post-pca matrix similarity
    testthat::expect_true(all.equal(lsi_mat_archr, lsi_bpcells$cell_embeddings, tolerance = 1e-4))
    # also check for correlation between the two matrices in PC space
    testthat::expect_true(cor(as.vector(lsi_mat_archr), as.vector(lsi_bpcells$cell_embeddings)) > 0.999)
}
test_lsi_similarity_to_archr()