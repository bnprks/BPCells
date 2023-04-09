# BPCells 1.0 Roadmap
- ~~Parallelization~~ (basic support complete. See below)
- Native python library (re-using C++ backend)
- Peak-gene correlations
- MACS peak calling

Contributions welcome :)

# BPCells 0.2.0 (github main branch - in progress)
## Features
- New `svds()` function, based on the excellent Spectra C++ library (used in RSpectra) by Yixuan Qiu.
  This should ensure lower memory usage compared to `irlba`, while achieving similar speed + accuracy.
- Limited parallelization is now supported. This is easiest to use via the `threads` argument to 
  `matrix_stats()` and `svds()`.
    - All normalizations are supported, but a few operations like `marker_features()` and writing a
      matrix to disk remain single-threaded.
    - Running `svds()` with many threads on gene-major matrices can result in high memory usage for now.
      This problem is not present for cell-major matrices.
- Reading text-based MatrixMarket inputs (e.g. from 10x or Parse) is now supported via
  `import_matrix_market()` and the convenience function `import_matrix_market_10x()`. Our
  implementation uses disk-backed sorting to allow importing large files with low memory usage.
- Added `binarize()` function and associated generics `<`, `<=`, `>`, and `>=`.
  This only supports comparison with non-negative numbers currently. (Thanks to 
  contribution from @brgew)
- Added `round()` matrix transformation (Thanks to contributions from @brgew)
- Add experimental internal getter/setter function `all_matrix_inputs()` to help enable relocating
  the underlying storage for BPCells matrix transform objects.

## Improvements
- Merging fragments with `c()` now handles inputs with mismatched chromosome names.
- Merging fragments is now 2-3.5x faster 
- SNN graph construction in `knn_to_snn_graph()` should work more smoothly on large datasets due to C++ implementation
- Reduced memory usage in `marker_features()` for samples with millions of cells and a large number
  of clusters to compare.
- On Windows, increased the maximum number of files that can be simultaneously open. Previously, opening >63 compressed
  counts matrices simultaneously would hit the limit. Now at least 1,000 simultaneous matrices should be possible.
- Subsetting peak or tile matrices with `[` now propagates through so we always avoid computing parts of
  the peak/tile matrix that have been discarded by our subset. Subsetting a tile matrix will automatically
  convert into a peak matrix when possible for improved efficiency.

## Bug-fixes
- Fixed a few fragment transforms where using `chrNames(frags) <- val` or `cellNames(frags) <- val` could cause
  downstream errors.
- Fixed errors in `transpose_storage_order()` for matrices with >4 billion non-zero entries.
- Fixed error in `transpose_storage_order()` for matrices with no non-zero entries.
- Fixed bug writing fragment files with >512 chromosomes.
- Fixed bug when reading fragment files with >4 billion fragments.
- Fixed file permissions errors when using read-only hdf5 files (Issue #26 reported thanks to @ttumkaya)
- Renaming `rownames()` or `colnames()` is now propagated when saving matrices (Issue #29 reported thanks to @realzehuali)
- Fixed 64-bit integer overflow (!) that could cause incorrect p-value calculations in `marker_features()` for features with
  more than 2.6 million zeros.
- Improved robustness of the Windows installation process for setups that do not need the -lsz linker flag to compile hdf5
- Fixed possible memory safety bug where wrapped R objects (such as dgCMatrix) could be potentially garbage collected
  while C++ was still trying to access the data in rare circumstances.

# BPCells 0.1.0

## Features
- ATAC-seq Analysis
    - Reading/writing 10x fragment files on disk
    - Reading/writing compressed fragments on disk (in folder or hdf5 group)
    - Interconversion of fragments objects with GRanges / data.frame
    - Merging of multiple source fragment files transparently at run time
    - Calculation of Cell x Peak matrices, and Cell x Tile matrices
    - ArchR-compatible QC calculations
    - ArchR-compatible gene activity score calculations
    - Filtering fragments by chromosmes, cells, lengths, or genomic region
    - Fast peak calling approximation via overlapping tiles
- Single cell matrices
    - Conversion to/from R sparse matrices
    - Read-write access to 10x hdf5 feature matrices, and read-only access to AnnData files
    - Reading/writing of compressed matrices on disk (in folder or hdf5 group)
    - Support for integer or single/double-precision floating point matrices on disk
    - Fast transposition of storage order, to switch between indexing by cell or
      by gene/feature.
    - Concatenation of multiple source matrix files transparently at run time
    - Single-pass calculation of row/column mean and variance
    - Wilcoxon marker feature calculation
    - Transparent handling of vector `+`, `-`, `*`, `/`, and `log1p` for streaming
      normalization, along with other less common operations. This allows implementation of ATAC-seq LSI and Seurat default
      normalization, along with most published log-based normalizations.
    - SCTransform pearson residual calculation
    - Multiplication of sparse matrices
- Single cell plotting utilities
    - Read count knee cutoffs
    - UMAP embeddings
    - Dot plots
    - Transcription factor footprinting / TSS profile plotting
    - Fragments vs. TSS Enrichment ATAC-seq QC plot
    - Pseudobulk genome track plots, with gene annotation plots
- Additional utility functions
    - Matching gene symbols/IDs to canonical symbols
    - Download transcript annotations from Gencode or GTF files 
    - Download + parse UCSC chromosome sizes
    - Parse peak files BED format; Download ENCODE blacklist region
    - Wrappers for knn graph calculation + clustering


Note: All operations interoperate with all storage formats. For example, all matrix operations can be applied directly to an AnnData or 10x matrix file. In many cases
the bitpacking-compressed formats will provide performance/space advantages, but
are not required to use the computations.