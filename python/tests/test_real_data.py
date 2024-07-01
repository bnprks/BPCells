# Copyright 2023 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.



# %%
import tempfile
import os

import pandas as pd
import numpy as np
import scipy
import h5py

import bpcells.experimental

import pytest

import utils

@utils.slow_data_test
def test_500_pbmc_matrix(tmp_path, fetch_cached_file):
    """ 
    This test takes about 15 seconds to run once data downloads are cached

    1. On a 500 PBMC 10x dataset, calculate 1kb of per-base coverage in all non-overlapping peaks of < 1kb,
       subset to chr20-chr22
    2. Compare the total counts per peak against the (hopefully identical) counts from 10x.
       (Note since each peak is variable width, this does exercise some checking of the per-base information)
    3. Also calculate per-peak coverage for full 1kb peaks, and confirm the sums are consistent
    """
    
    # Download inputs from 10x
    fragments_path = fetch_cached_file("https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_500_nextgem/atac_pbmc_500_nextgem_fragments.tsv.gz")
    peaks_path = fetch_cached_file("https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_500_nextgem/atac_pbmc_500_nextgem_peaks.bed")
    peak_matrix_path = fetch_cached_file("https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_500_nextgem/atac_pbmc_500_nextgem_raw_peak_bc_matrix.h5")
    barcode_metrics_path = fetch_cached_file("https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_500_nextgem/atac_pbmc_500_nextgem_singlecell.csv")
    chrom_sizes_path = fetch_cached_file("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes")


    # Convert BPCells fragments
    bpcells_path = os.path.join(tmp_path, "fragments")
    bpcells.experimental.import_10x_fragments(fragments_path, bpcells_path, shift_end=1)

    # Get a non-overlapping subset of the peaks
    peaks = pd.read_csv(peaks_path, sep="\t", comment="#", names=["chrom", "start", "end"])
    peaks["orig_width"] = peaks.end - peaks.start
    max_width = 1000
    peaks.end = peaks.start + max_width
    too_wide = peaks.orig_width > max_width
    overlaps = np.concatenate(([False], (peaks.end.to_numpy()[:-1] > peaks.start.to_numpy()[1:]) & (peaks.chrom.to_numpy()[:-1] == peaks.chrom.to_numpy()[1:])))

    chr_subset = ["chr20", "chr21", "chr22"]
    peak_subset = peaks[~too_wide & ~overlaps & peaks.chrom.isin(chr_subset)]

    # Choose dummy "clusters" based on the first two letters of the barcode
    cell_metadata = pd.read_csv(barcode_metrics_path).iloc[1:]
    cell_metadata = cell_metadata[cell_metadata.passed_filters >= 1]
    clusts = cell_metadata.barcode.str.slice(0,2)
    group_order = sorted(clusts.unique())

    groups = bpcells.experimental.build_cell_groups(
        bpcells_path,
        cell_ids = cell_metadata.barcode,
        group_ids = clusts,
        group_order = group_order
    )

    # Compute the ground-truth pseudobulk peak counts fro the 10x matrix
    h5 = h5py.File(peak_matrix_path, "r")

    mat_10x = scipy.sparse.csr_matrix(
        (h5["matrix/data"][:], h5["matrix/indices"][:], h5["matrix/indptr"][:])
    )
    cell_clust_mat = np.zeros((groups.max() + 1, cell_metadata.shape[0]))
    cell_group_idx = [group_order.index(b[:2]) for b in cell_metadata.barcode]
    cell_clust_mat[(cell_group_idx, range(cell_metadata.shape[0]))] = 1

    pseudobulk_mat = cell_clust_mat @ mat_10x
    pseudobulk_mat = pseudobulk_mat[:,~too_wide & ~overlaps & peaks.chrom.isin(chr_subset)]

    # Calculate the BPCells matrices
    basepair_matrix = bpcells.experimental.pseudobulk_insertion_counts(
        bpcells_path,
        peak_subset,
        groups,
        bin_size = 1
    )

    peak_matrix = bpcells.experimental.pseudobulk_insertion_counts(
        bpcells_path,
        peak_subset,
        groups,
        bin_size = max_width
    )

    # Check that the results match with the 10x version
    for i in range(pseudobulk_mat.shape[1]):
        peak_width = peak_subset.orig_width.iloc[i]
        assert np.all(basepair_matrix[i, :,:peak_width].sum(axis=-1) == pseudobulk_mat[:,i])

    assert np.all(basepair_matrix.sum(axis=-1) == peak_matrix[:,:,0])

    # Make a precalculated insertion matrix
    precalculated_path = os.path.join(tmp_path, "precalculated_mat")
    chrom_sizes = pd.read_csv(chrom_sizes_path, sep="\t", names=["chrom", "size"])
    chrom_sizes = {t.chrom: t.size for t in chrom_sizes.itertuples() if t.chrom in chr_subset}
    bpcells.experimental.precalculate_insertion_counts(
        bpcells_path,
        precalculated_path,
        groups,
        chrom_sizes,
        threads=4
    )
    m = bpcells.experimental.PrecalculatedInsertionMatrix(precalculated_path)
    assert m.shape[0] == len(group_order)
    basepair_precalculated = m.get_counts(peak_subset)
    assert np.all(basepair_precalculated == basepair_matrix)
