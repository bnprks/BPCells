# Copyright 2024 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

source("test_helpers.R")
devtools::load_all(config[["path_bpcells"]])
devtools::load_all(config[["path_archr"]])

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
    lsi_archr <- .computeLSI(
        mat = test_mat,
        LSIMethod = 2,
        nDimensions = 2,
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
    lsi_bpcells <- lsi(
        test_mat %>% as("dgCMatrix") %>% as("IterableMatrix"), 
        z_score_norm = FALSE,
        n_dimensions = 2,
        save_lsi = TRUE
    )
    pre_svd_mat_approx_bpcells <- lsi_bpcells$svd_attr$u %*% lsi_bpcells$pca_res
    testthat::expect_true(all.equal(pre_svd_mat_approx_archr, pre_svd_mat_approx_bpcells, tolerance = 1e-6))
}
test_lsi_similarity_to_archr()