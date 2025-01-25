# Copyright 2024 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

import sys, tempfile, os
import numpy as np
import pandas as pd
import scanpy as sc


def call_highly_var_genes_single_batch(temp_dir: str) -> None:
    """
    Call highly_variable genes on a single batch of PBMC68k data using their interpreation of 
    the Seurat implementation.
    Write the input anndata object csv at `<temp_dir>/highly_var_genes_scanpy_input.csv`
    Write the output as a csv, at `<temp_dir>/highly_var_genes_scanpy_output.csv`

    Args:
        temp_dir (str): Path to the temporary directory to write the input and output files.
    """
    # Dataset is only (765, 700)
    adata = sc.datasets.pbmc68k_reduced()
    adata.var_names_make_unique()
    res = sc.pp._highly_variable_genes.highly_variable_genes(adata, 
                                                             n_top_genes = 50,
                                                             n_bins = 20, 
                                                             flavor = "seurat",
                                                             inplace = False,
                                                             check_values = False).drop(columns = 'mean_bin')
    # remove mean_bin
    adata.to_df().to_csv(os.path.join(temp_dir, "highly_var_genes_scanpy_input.csv"))
    res.to_csv(os.path.join(temp_dir, "highly_var_genes_scanpy_output.csv"))


if __name__ == "__main__":
    # We use the first argument as the temporary directory
    if len(sys.argv) < 2:
        # If no argument is provided, use the current directory
        call_highly_var_genes_single_batch(".")
    else:
        call_highly_var_genes_single_batch(sys.argv[1])
    