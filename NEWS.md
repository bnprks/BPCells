# BPCells 1.0 Roadmap
- ~~Parallelization~~ (basic support complete. See below)
- Native python library (re-using C++ backend)
- Peak-gene correlations
- MACS peak calling

Contributions welcome :)

# BPCells 0.2.0 (github main branch - in progress)
## Breaking changes
- New slots have been added to 10x matrix objects, so any saved RDS files may need to have
  their 10x matrix inputs re-opened and replaced by calling `all_matrix_inputs()`. Outside of
  loading old RDS files no changes should be needed.

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
- Add getter/setter function `all_matrix_inputs()` to help enable relocating
  the underlying storage for BPCells matrix transform objects.
- All hdf5-writing functions now support a `gzip_level` parameter, which will enable a shuffle + gzip filter for
  compression. This is generally much slower than bitpacking compression, but it adds improved storage options for
  files that must be read by outside programs. Thanks to @ycli1995 for submitting this improvement in pull #42.
- AnnData export now supported via `write_matrix_anndata_hdf5()` (issue #49)
- Re-licensed code base to use dual-licensed Apache V2 or MIT instead of GPLv3
- Assigning to a subset is now supported (e.g. `m1[i,j] <- m2`). Note that this does not modify data on disk. Instead,
  it uses a series of subsetting and concatenation operations to provide the *appearance* of overwriting the appropriate
  entries.
- Added `knn_to_geodesic_graph()`, which matches the Scanpy default construction for
  graph-based clustering

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
- Subsetting RowBindMatrices and ColBindMatrices now propagates through so we avoid touching matrices with no selected indices
- Added logic to help reduce cases where subsetting causes BPCells to fall back to a less efficient matrix-vector multiply algorithm.
  This affects most math transforms. As part of this, the filtering part of a subset will propagate to earlier transformation steps, while the
  reordering will not. Thanks to @nimanouri-nm for raising issue #65 to fix a bug in the initial implementation.
- Additional C++17 filesystem backwards compatibility that should allow slightly older compilers such as GCC 7.5 to 
  build BPCells.
- `as.matrix()` will produce integer matrices when appropriate (Thanks to @Yunuuuu in pull #77)
- 10x HDF5 matrices can now read and write non-integer types when requested (Thanks to @ycli1995 in pull #75)
- Old-style 10x files from cellranger v2 can now read multi-genome files, which are returned as a list (Thanks to @ycli1995 in pull #75)

## Bug-fixes
- Fixed a few fragment transforms where using `chrNames(frags) <- val` or `cellNames(frags) <- val` could cause
  downstream errors.
- Fixed errors in `transpose_storage_order()` for matrices with >4 billion non-zero entries.
- Fixed error in `transpose_storage_order()` for matrices with no non-zero entries.
- Fixed bug writing fragment files with >512 chromosomes.
- Fixed bug when reading fragment files with >4 billion fragments.
- Fixed file permissions errors when using read-only hdf5 files (Issue #26 reported thanks to @ttumkaya)
- Renaming `rownames()` or `colnames()` is now propagated when saving matrices (Issue #29 reported thanks to @realzehuali, with an additional fix after report thanks to @Dario-Rocha)
- Fixed 64-bit integer overflow (!) that could cause incorrect p-value calculations in `marker_features()` for features with
  more than 2.6 million zeros.
- Improved robustness of the Windows installation process for setups that do not need the -lsz linker flag to compile hdf5
- Fixed possible memory safety bug where wrapped R objects (such as dgCMatrix) could be potentially garbage collected
  while C++ was still trying to access the data in rare circumstances.
- Fixed case when dimnames were not preserved when calling `convert_matrix_type()` twice in a row such that it cancels out (e.g. double -> uint32_t -> double). Thanks to @brgrew reporting issue #43
- Caused and fixed issue resulting in unusably slow performance reading matrices from HDF5 files. Broken versions range from commit  21f8dcf until the fix in 3711a40 (October 18-November 3, 2023). Thanks to @abhiachoudhary for reporting this in issue #53
- Fixed error with `svds()` not handling row-major matrices correctly. Thanks to @ycli1995 for reporting this in issue #55
- Fixed error with row/col name handling for AnnData matrices. Thanks to @lisch7 for reporting this in issue #57
- Fixed error with merging matrices of different data types. Thanks to @Yunuuuu for identifying the issue and providing a fix (#68 and #70)
- Fixed issue with losing dimnames on subset assignment `[<-`. Thanks to @Yunuuuu for identifying the issue #67
- Fixed incorrect results with some cases of scaling matrix after shifting. Thanks to @Yunuuuu for identifying the issue #72
- Fixed infinite loop bug when calling `transpose_storage_order()` on a densely-transformed matrix. Thanks to @Yunuuuu for reporting this in issue #71
- h5ad outputs will now subset properly when loaded by the Python anndata package (Thanks to issue described by @ggruenhagen3 in issue #49 and fixed by @ycli1995 in pull #81)

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