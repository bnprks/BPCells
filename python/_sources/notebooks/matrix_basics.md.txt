---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.15.1
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Matrix slicing tutorial

BPCells has prototype Python bindings that allow for matrix creation and slicing, with optional multithreaded reads.

## Performance estimates
On the scimilarity dataset of 15M human cells, compressed storage is 64GB (2.2 bytes/non-zero)

Read speeds for 10k random cells from the 15M human cells (range of 5 random tests)

| Storage location         | 1 thread     | 4 threads |
|--------------|-----------|------------|
| Memory | 2.8-4.7 seconds      | 1.0-1.1 seconds |
| Local SSD      | 4.5-4.9 seconds  | 1.5-1.7 seconds       |
| Networked FS (warm cache)  | 20-21 seconds  | 5.5-6.2 seconds       |
| Networked FS (cold cache) | 🙁 | 76-115 seconds |

+++

## Demo data setup

```{code-cell} ipython3
import bpcells.experimental

import os
import tempfile

import numpy as np
import scipy.sparse
```

```{code-cell} ipython3
tmp = tempfile.TemporaryDirectory()
os.chdir(tmp.name)
```

```{code-cell} ipython3
mat = scipy.sparse.csc_matrix(np.array([
    [1, 0, 4, 0],
    [0, 0, 5, 7],
    [2, 3, 6, 0]]
))
mat
```

```{code-cell} ipython3
mat.toarray()
```

## Basic usage from scipy.sparse

```{code-cell} ipython3
bp_mat = bpcells.experimental.DirMatrix.from_scipy_sparse(mat, "basic_mat")
bp_mat
```

Slicing the matrix returns a `scipy.sparse` matrix

```{code-cell} ipython3
bp_mat[:,:]
```

```{code-cell} ipython3
bp_mat[:,:].toarray()
```

We can use many of the same slicing options as standard numpy matrices

```{code-cell} ipython3
bp_mat[1:3, [0,2]].toarray()
```

```{code-cell} ipython3
bp_mat[[True, False, True], -2:].toarray()
```

We can also make a transposed view of the matrix similar to numpy. No work is done, we just switch between row-major and col-major representations

```{code-cell} ipython3
bp_mat.T
```

## Reopening the matrix later

+++

The matrix path has 13 files (for compressed integer matrices), which contain data and metadata

```{code-cell} ipython3
!ls -l basic_mat
```

```{code-cell} ipython3
bp_mat = bpcells.experimental.DirMatrix("basic_mat")
```

## Import from h5ad

```{code-cell} ipython3
import anndata
anndata.AnnData(mat).write("mat.h5ad")
bp_mat = bpcells.experimental.DirMatrix.from_h5ad("mat.h5ad", "basic_mat_from_h5ad")
bp_mat[:,:].toarray()
```

## Concatenate multiple matrices

+++

We can concatenate multiple matrices to a single file on disk with low memory usage. This allows importing many samples in parallel, then concatenating them together into a single matrix

```{code-cell} ipython3
bpcells.experimental.DirMatrix.from_hstack(
    [bp_mat, bp_mat], 
    "basic_mat_hstack"
)[:,:].toarray()
```

```{code-cell} ipython3
bpcells.experimental.DirMatrix.from_vstack(
    [bp_mat, bp_mat], 
    "basic_mat_vstack"
)[:,:].toarray()
```

## Multithreaded operation

For larger matrices, it can be desirable to perform matrix reading in a multi-threaded manner. When using multiple threads, BPCells will divide the matrix slice query into chunks that are loaded in parallel, then recombined in memory after all threads are completed.

When performing random slicing along the major storage axis, seek latency is the primary performance bottleneck. Setting a high number of threads (even above the actual core count of the machine) can help mitigate filesystem seek latency.

When slicing across the non-major storage axis, decompression speed can become the performance bottleneck. Setting threads to the number of available cores can help to parallelize the decompression speed. For cell-major RNA-seq matrices, each thread can process compressed input at a rate of about 1 GB/s, so filesystems with >1GB/s sequential read speeds will benefit from some parallelization.

```{code-cell} ipython3
bp_mat.threads = 8
bp_mat[:,:].toarray() # This will be performed with 8 threads now
```

## Compressed in-memory storage

For neural network training use-cases, fast slicing performance may be critical to avoid bottlenecking on data loads. In this case, BPCells supports loading the compressed data in memory, which eliminates seek latency while saving ~4x memory usage compared to an uncompressed scipy sparse matrix.

Loading can only be performed from an existing BPCells matrix directory, and the current version involves re-compressing the data in-memory at load time (avoidable, but a bit trickier to code so direct loading isn't implemented yet)

```{code-cell} ipython3
bp_mat_mem = bpcells.experimental.MemMatrix("basic_mat")
bp_mat_mem.threads = 8
bp_mat_mem[:,:].toarray()
```
