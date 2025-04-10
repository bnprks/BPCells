# Copyright 2024 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

library("BPCells")

# Set up temp dir in case it's not already set
create_temp_dir <- function(dir = NULL) {
  if (is.null(dir)) {
    dir <- file.path(tempdir(), "lsi_test")
    if (dir.exists(dir)) unlink(dir, recursive = TRUE)
    dir.create(dir)
  }
  return(dir)
}

# Compare the feature selection output of BPCells to that of Scanpy.
# Scanpy technically utilizes the Seurat (Satija et al. 2015) method for feature selection, so we should expect similar results of either pkg. 
# This function calls a python script that runs Scanpy feature selection on a test dataset, and writes both input/output to `dir`.
# It then reads in the input/output from the python script, calls the BPCells feature selection function, and compares the output to the Scanpy output.
compare_feat_selection_to_scanpy <- function(dir = NULL) {
    dir <- create_temp_dir(dir)

    # Call python script
    system2("python3", c("Scanpy_variable_feat_selection.py", dir))
    
    # read in input csv
    input_mat_scanpy <- t(read.csv(file.path(dir, "highly_var_genes_scanpy_input.csv"), row.names = 1))
    output_mat_scanpy <- read.csv(file.path(dir, "highly_var_genes_scanpy_output.csv"), row.names = 1)
    # filter output mat to only where "highly_variable" is true
    output_mat_scanpy$highly_variable <- as.logical(output_mat_scanpy$highly_variable)
    output_mat_scanpy <- output_mat_scanpy[output_mat_scanpy$highly_variable,] %>% 
        dplyr::arrange(desc(dispersions_norm)) %>% 
        dplyr::select(-highly_variable) %>% # convert rownames to a column
        tibble::rownames_to_column("name") %>%
        dplyr::as_tibble()

    # Scanpy undoes a log1p transformation on the input matrix, so we do the same here
    input_mat_bpcells <- expm1(input_mat_scanpy)
    
    output_bpcells <- highly_variable_features(
        input_mat_bpcells %>% as("dgCMatrix") %>% as("IterableMatrix"),
        num_feats = 50,
        n_bins = 20,
        save_feat_selection = TRUE
    )
    output_mat_bpcells <- output_bpcells$feature_selection
    expect_true(all.equal(output_mat_bpcells$name, output_mat_scanpy$name))
    expect_true(all.equal(output_mat_bpcells$mean, output_mat_scanpy$means, tolerance = 1e-6))
    expect_true(all.equal(output_mat_bpcells$dispersion, output_mat_scanpy$dispersions, tolerance = 1e-6))
    expect_true(all.equal(output_mat_bpcells$feature_dispersion_norm, output_mat_scanpy$dispersions_norm, tolerance = 1e-6))
}

compare_feat_selection_to_scanpy()
