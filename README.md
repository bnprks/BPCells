# BPCells

BPCells is a package for high performance single cell analysis. It is designed to
cover the processing pipeline from ATAC fragments or RNA counts matrices through
to normalization, basic QC, and PCA. 

Three of the key distinguishing features that allow for high performance in BPCells are:

1. Bit-packing compression to allow for extremely compact storage of
   fragments or counts matrices on-disk or in memory.
2. C++ code that operates on all data in a streaming fashion to support low
   memory usage and efficient use of CPU cache.
3. Matrix-free SVD solvers in combination with implicit normalization calculations 
  to support computing the PCA of a normalized matrix while only ever storing the 
  original counts matrix in memory.

## Installation

BPCells is easiest to install directly from github:

```R
devtools::install_github("bnprks/BPCells")
```

## Getting started

Two key principles to understand about using BPCells is that all operations are
*streaming* and *lazy*. 

Streaming means that only a minimal amount of data is 
stored in memory while a computation is happening. There is almost no
memory used storing intermediate results. Hence, we can compute operations 
on large matrices without ever loading them into memory.

Lazy means that no real work is performed on matrix or fragment objects until
the result needs to be returned as an R object or written to disk. This helps support
the streaming computation, since otherwise we would be forced to compute intermediate
results and use additional memory.

### Basic usage
We begin with a basic example of loading ATAC fragments from a 10x fragments file,
reading a peak set from a bed file, then calculating a cell x peak matrix.
```R
library("BPCells")

# File reading is lazy, so this is instantaneous
fragments <- open_fragments_10x("atac_fragments.tsv.gz")

# This is when we actually read the file, should take 1-2 minutes to scan
# since we bottleneck on gzip decompression.
packed_fragments <- write_fragments_memory(fragments)

# Important to set compress=FALSE for speed. Should take a few seconds
saveRDS(packed_fragments, "fragments.rds", compress=FALSE)

# Reloading from disk is only a few seconds now.
packed_fragments <- readRDS("fragments.rds")

peaks <- read_bed("peaks.bed")

# This is fast because the peak matrix calculation is lazy.
# It will be computed on-the-fly when we actually need results from it.
peakMatrix <- peakMatrix(packed_fragments, peaks)

# Here is where the peak matrix calculation happens. Should take
# under 10 seconds.
R_matrix <- as(peakMatrix, "dgCMatrix")
```

### Streaming operations

The lazy, stream-oriented design means that we can calculate more complicated
transformations in a single pass. This is both faster and more memory-efficient
than calculating several intermediate results in a sequential manner.

As an example, we will perform the following pipeline:
1. Exclude fragments from non-standard chromosomes
2. Subset our cells
3. Add a Tn5 offset
4. Calculate the peakMatrix
5. Calculate the mean-accessibility per peak

If this were done using e.g. GRanges or sparse matrices, we would need to do 3
passes through the fragments while saving intermediate results, and 2 passes over
the peakMatrix.

With BPCell's streaming operations, this can all be done directly from the fragments in a single pass, and the memory
usage is limited to a few bytes per cell for iterating over the peakMatrix 
and returning the colMeans.
```R
# Here I make use of the pipe operator (%>%) for better readability
library("tidyverse")

# We'll subset to just the standard chromosomes
standard_chr <- which(
  stringr::str_detect(chrNames(packed_fragments), "^chr[0-9XY]+$")
)

# Pick a random subset of 100 cells to consider
set.seed(1337)
keeper_cells <- sample(cellNames(packed_fragments), 100)

# Run the pipeline, and save the average accessibility per peak
peak_accessibility <- packed_fragments %>%
  select_chromosomes(standard_chr) %>%
  select_cells(keeper_cells) %>%
  shift_fragments(shift_start=4, shift_end=-5) %>%
  peakMatrix(peaks) %>%
  colMeans()
```

Note that if we knew the cell names ahead of time, we could even perform this
operation directly on our orignal 10x fragments without ever saving the
fragments into memory. This would be fairly slow because 10x fragment files are
slow to decompress. With upcoming support for storing packed fragments directly
on disk, this can become a much faster operation without ever needing to store
fragments in memory.

## Roadmap

### Current support:
- Fragments
    - Reading/writing 10x fragment files on disk
    - Reading/writing packed fragment objects in memory or directly from disk
    - Interconversion of fragments objects with GRanges
    - Calculation of Cell x Peak matrices
- Matrices
    - Conversion to/from R sparse matrices
    - Read-write access to 10x hdf5 feature matrices, and read-only access to AnnData files, and 
    - Reading/writing of packed sparse matrices in memory or directly from disk
    - Multiplication by dense matrices or vectors
    - Calculation of statistics like rowSums, colSums, rowMeans, and colMeans
    - Transparent handling of vector `+`, `-`, `*`, `/`, and `log1p` for streaming
      normalization. This allows implementation of ATAC-seq LSI and Seurat default
      normalization.

### Upcoming additions:
- Support for additional fragment formats:
    - Read fragments from bam files
    - Support direct download of files from URLs
- Support for additional matrix formats:
    - Write hdf5 AnnData matrices
- Support for additional matrix normalizations:
    - sctransform normalization

### Performance estimates
- Bit-packed storage:
    - Packed fragments are about 2x smaller than ArchR fragments and gzipped 10x fragment files
    - Packed fragments can be decompressed at >5GB/s, and so decoding is disk-limited on
      SSDs or RAID arrays while remaining competitive with reading directly from uncompressed 
      RAM
    - Packed matrices are 4-6x smaller than the equivalent sparse matrices,
      with similarly excellent decompression speed. 
      
      (Note: compared to R dgCMatrices,
      the numbers are more like 6-8x smaller due to dgCMatrices always storing
      double-precision floats)
- Comparison to ArchR:
    - Fragment import speed is about 10-fold faster
    - Peak and tile matrix calculations are about 10-fold faster than ArchR
    - BPCells usually uses less memory than ArchR, though ArchR's memory usage is
      never excessive
- PCA calculation:
    - Compared to Seurat/Scanpy's default normalization + PCA, BPCells can be nearly 100x more
      memory efficient for large datasets (1M cells), even when subsetting to just 2,000 genes.
    - The time for PCA calculation is faster than Seurat, but about 2x slower than Scanpy 
      (BPCells) takes about 5 minutes to run PCA on a 1M cell matrix
    - Compared to Seurat's default normalization + PCA, BPCells will likely be about
      10x more efficient in memory and CPU. It is unclear if this will multiply with
      the 4-6x memory savings from bitpacking counts matrices.
    - LSI is not expected to be substantially faster than ArchR, although the C++
      implementation may provide a several-fold speedup, and bitpacking will provide
      a 4-6x memory savings.