# BPCells 0.3.1 (in-progress main branch)

## Breaking changes
- Change first parameter name of `cluster_graph_leiden()`, `cluster_graph_louvain()` and `cluster_graph_seurat()` from `snn` to `mat` to more accurately reflect the input type.  (pull request #189)

## Features
- Add `write_matrix_anndata_hdf5_dense()` which allows writing matrices in AnnData's dense format, most commonly used for `obsm` or `varm` matrices. (Thanks to @ycli1995 for pull request #166)
- Add normalization helper functions `normalize_log()` and `normalize_tfidf()` (pull request #168)
- Add functions `normalize_tfidf()` and `normalize_log()`, which allow for easy normalization of iterable matrices using TF-IDF or log1p(pull request #189)
- Add feature selection functions `select_features_variance()`, and `select_features_{dispersion,mean,binned_dispersion}()`, with parameterization for normalization steps, and number of variable features (pull request #189)
- Add `LSI()` and `IterativeLSI()` dimensionality functions to perform latent semantic indexing on a matrix (pull request #189).
- Add capability to create partial function objects in when excluding the first argument of a function.  This is implemented in normalizations, feature selections, dimensionality reductions, and clustering functions. See `select_features_variance()` for usage.  (pull request #189)
- Create a wrapper function `cluster_cells_graph()` that wraps the steps of knn object creation, graph adjacency creation, and clustering all within a single function (pull request #189)

## Improvements
- Speed up taking large subsets of large concatenated matrices, e.g. selecting 9M cells from a 10M cell matrix composed of ~100 concatenated pieces. (pull request #179)
- `matrix_stats()` now also works with types `matrix` and `dgCMatrix`. (pull request #190)
- Fixed memory errors when running `writeInsertionBed()` and `writeInsertionBedGraph()` (pull request #{118, 134})
- Export `merge_peaks_iterative()`, which helps create non-overlapping peak sets.  (pull request #216)


## Bug-fixes
- Fix error message printing when MACS crashes during `call_peaks_macs()` (pull request #175)
- Fix `gene_score_archr()` and `gene_score_weights_archr()` malfunctioning for non-default `tile_width` settings. (Thanks to @Baboon61 for reporting issue #185)
- Fix `gene_score_archr()` when `chromosome_sizes` argument is not sorted. (Thanks to @Baboon61 for reporting issue #188)
- Fix matrix transpose error when BPCells is loaded via `devtools::load_all()` and `BiocGenerics` has been imported previously. (pull request #191)
- Fix error when using a single group in `write_insertion_bedgraph()` (pull request #214)
- Fix GRanges conversion functions sometimes not being defined if BPCells is built as a binary package prior to GenomicRanges being installed. (pull request #231; thanks to @mfansler for reporting issue #229)
- Fix error in `write_matrix_hdf5()` when overwriting to a `.h5` file that does not exist. (pull request #234)

# BPCells 0.3.0 (12/21/2024)

The BPCells 0.3.0 release covers 6 months of changes and 45 commits from 5 contributors. Notable improvements
this release include support for peak calling with MACS and the addition of pseudobulk matrix and stats calculations.
We also released an initial prototype of a BPCells Python library (more details [here](https://bnprks.github.io/BPCells/python/index.html)).
Full details of changes below.

Thanks to @ycli1995, @Yunuuuu, and @douglasgscofield for pull requests that contributed to this release, as well as to users who
sumitted github issues to help identify and fix bugs. We also added @immanuelazn to the team as a new hire! He is responsible for many
of the new features this release and will continue to help with maintenance and new development moving forwards. 

## Features
- `apply_by_col()` and `apply_by_row()` allow providing custom R functions to compute per row/col summaries.
  In initial tests calculating row/col means using R functions is ~2x slower than the C++-based implementation but memory
  usage remains low.
- Add `rowMaxs()` and `colMaxs()` functions, which return the maximum value in each row or column of a matrix. 
  If `matrixStats` or `MatrixGenerics` packages are installed, `BPCells::rowMaxs()` will fall back to their implementations for non-BPCells objects.
  Thanks to @immanuelazn for their first contribution as a new lab hire!
- Add `regress_out()` to allow removing unwanted sources of variation via least squares linear regression models.
  Thanks to @ycli1995 for pull request #110
- Add `trackplot_genome_annotation()` for plotting peaks, with options for directional arrows, colors, labels, and peak widths. (pull request #113)
- Add MACS2/3 input creation and peak calling through `call_peaks_macs()`(pull request #118). Note, renamed from `call_macs_peaks()` in pull request #143
- Add `rowQuantiles()` and `colQuantiles()` functions, which return the quantiles of each row/column of a matrix. Currently `rowQuantiles()` only works on row-major matrices and `colQuantiles()` only works on col-major matrices.
  If `matrixStats` or `MatrixGenerics` packages are installed, `BPCells::colQuantiles()` will fall back to their implementations for non-BPCells objects. (pull request #128)
- Add `pseudobulk_matrix()` which allows pseudobulk aggregation by `sum` or `mean` and calculation of per-pseudobulk `variance` and `nonzero` statistics for each gene (pull request #128)

## Improvements
- `trackplot_loop()` now accepts discrete color scales
- `trackplot_combine()` now has smarter layout logic for margins, as well as detecting when plots are being combined that cover different genomic regions. (pull request #116)
- `select_cells()` and `select_chromosomes()` now also allow using a logical mask for selection. (pull request #117)
- BPCells installation can now also be configured by setting the `LDFLAGS` or `CFLAGS` as environment variables in addition to setting them in `~/.R/Makevars` (pull request #124)
- `open_matrix_anndata_hdf5()` now supports reading AnnData matrices in the dense format. (pull request #146)
- `cluster_graph_leiden()` now has better defaults that produce reasonable cluster counts regardless of dataset size. (pull request #147) 

## Bug-fixes
- Fixed error message when a matrix is too large to be converted to dgCMatrix. (Thanks to @RookieA1 for reporting issue #95)
- Fixed forgetting dimnames when subsetting after certain sets of
  operations. (Thanks to @Yunuuuu for reporting issues #97 and #100)
- Fixed plotting crashes when running `trackplot_coverage()` with fragments from a single cluster. (Thanks to @sjessa for directly reporting this bug and coming up with a fix)
- Fixed issues with `trackplot_coverage()` when called with ranges less than 500 bp in length (Thanks to @bettybliu for directly reporting this bug.)
- Fix Rcpp warning created when handling compressed matrices with only one non-zero entry (pull request #123)
- Fixed discrepancy between default ArchR and BPCells peak calling insertion method, where BPCells defaulted to only using the start of each fragment as opposed to ArchR's method of using both start and end sites of fragments (pull request #143)
- Fix error in `tile_matrix()` with fragment mode (pull request #141)
- Fix precision bug in `sctransform_pearson()` on ARM architecture (pull request #141) 
- Fix type-confusion error when `pseudobulk_matrix()` gets an integer matrix (pull request #174)

## Deprecations
- `trackplot_coverage()` `legend_label` argument is now ignored, as the color legend is no longer shown by default for coverage plots.

# BPCells 0.2.0 (6/14/2024)

We are finally declaring a new release version, covering a large amount of changes and improvements
over the past year. Among the major features here are parallelization options for `svds()` and 
`matrix_stats()`, improved genomic track plots, and runtime CPU feature detection for SIMD code (enables
higher performance, more portable builds). Full details of changes below.

This version also comes with a new installation path, which is done in preparation for a future 
Python package release. (So we can have one folder for R and one for Python, rather than having all
the R files sit in the root folder). This is a breaking change and requires a slightly
modified installation command.

Thanks to @brgew, @ycli1995, and @Yunuuuu for pull requests that contributed to this release, as
well as all users who submitted github issues to help identify and fix bugs.

## Breaking changes
- Installation location has changed, to make room for a future python package release. New
  installs will have to use `remotes::install_github("bnprks/BPCells/r")` (note the additional `/r`)
  - r-universe mirrors will have to add `"subdir": "r"` to their `packages.json` config.
- New slots have been added to 10x matrix objects, so any saved RDS files may need to have
  their 10x matrix inputs re-opened and replaced by calling `all_matrix_inputs()`. Outside of
  loading old RDS files no changes should be needed.
- `trackplot_gene()` now returns a plot with a facet label to match the new trackplot system.
  This label can be removed by by calling `trackplot_gene(...) + ggplot2::facet_null()` to be 
  equivalent to the old function's output.

## Deprecations
- `draw_trackplot_grid()` deprecated, replaced by `trackplot_combine()` with simplified arguments
- `trackplot_bulk()` has been deprecated, replaced by `trackplot_coverage()` with equivalent functionality
- The old function names will output deprecation warnings, but otherwise work as before. 

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
- Add `checksum()`, which allows for calculating an MD5 checksum of a matrix contents. Thanks to @brgrew for submitting this improvement in pull request #83
- `write_insertion_bedgraph()` allows exporting pseudobulk insertion data to bedgraph format

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
- Trackplots have received several improvements
  - Trackplots now use faceting to provide per-plot labels, leading to an easier-to-use `trackplot_combine()` 
  - `trackplot_gene()` now draws arrows for the direction of transcription
  - `trackplot_loop()` is a new track type allows plotting interactions between genomic regions, for instance peak-gene correlations
    or loop calls from Hi-C
  - `trackplot_scalebar()` is added to show genomic scale
  - All trackplot functions now return ggplot objects with additional metadata stored for the plotting height of each track
  - Labels and heights for trackplots can be adjusted using `set_trackplot_label()` and `set_trackplot_height()`
  - The getting started pbmc 3k vignette now includes the updated trackplot APIs in its final example
- Add `rowVars()` and `colVars()` functions, as convenience wrappers around `matrix_stats()`. 
  If `matrixStats` or `MatrixGenerics` packages are installed, `BPCells::rowVars()` will fall back to
  their implementations for non-BPCells objects. Unfortunately, `matrixStats::rowVars()` is not generic, so either `BPCells::rowVars()` or 
  `BPCells::colVars()`
- Optimize mean and variance calculations for matrices added to a per-row or per-column constant.
- Migrate SIMD code to use [`highway`](https://github.com/google/highway).
  - Adds run-time detection of CPU features to eliminate architecture-specific compilation
  - For now, the `Pow` SIMD implementation is removed, but `Square` gets a new SIMD implementation
  - Empirically, most operations using SIMD math instructions are about 2x faster. This includes `log1p()`, and `sctransform_pearson()`
  - Minor speedups on dense-sparse matrix multiply functions (1.1-1.5x faster)

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
- Disk-backed fragment objects now load via absolute path, matching the behavior of matrices and making it so objects 
  loaded via `readRDS()` can be used from different working directories.
- `footprints()` now respects user interrupts via Ctrl-C

# BPCells 0.1.0 (4/7/2023)

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