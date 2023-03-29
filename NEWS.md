# BPCells 1.0 Roadmap
- Parallelization
- Native python library (re-using C++ backend)
- Peak-gene correlations
- MACS peak calling

Contributions welcome :)

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