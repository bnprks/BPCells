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

# Basepair insertion counts tutorial

BPCells python bindings can be used to query basepair-level coverage for predefined cell types. 

The way this works is in two steps: 

1. 10x or ArchR arrow files are converted to BPCells format. This is most flexible to do with the BPCells R bindings, though single-sample 10x import is supported through the python bindings.
2. The BPCells python bindings use the input fragment files to create large matrix of dimensions (# cell types, # basepairs in genome). This is where cell type groupings are determined.
3. The BPCells python bindings can slice arbitrary genomic regions, returning a numpy array of dimensions (regions, cell types, basepairs)

## Benchmark estimates
Benchmark dataset: 600K cell subset of the [Catlas paper](https://www.sciencedirect.com/science/article/pii/S0092867421012794), with 2.5 billion fragments
Benchmark task: Load 128 random 501-bp peak regions from 111 cell types at basepair resolution
Storage location: Local SSD. Networked file systems will be slower

|               | BPCells                | BigWigs     |
|---------------|------------------------|-------------|
| Creation time | 4.7 minutes, 8 threads | ?           |
| File size     | 6.2 GB                 | 13 GB       |
| Query time    | 0.37 seconds           | 2.2 seconds |

## Main benefits of BPCells
- Cell type count aggregation can be re-run fully from Python
- Query time is about 6x faster than BigWigs

**Caveat for this prototype**: due to development time limitations, this insertion matrix implementation does not support >=2^32 non-zero entries (4.29 billion). The catlas dataset had about 3.2 billion non-zero entries. This limitation can be removed with additional technical work, or as a workaround multiple matrix objects can be created that each individually have <2^32 non-zero entries.

+++

# Usage Demo

```{code-cell} ipython3
import bpcells

import pandas as pd
```

### Data download
We use a public 500-cell 10x dataset

```{code-cell} ipython3
import os.path
import subprocess
import tempfile
tmpdir = tempfile.TemporaryDirectory()
fragments_10x_path = os.path.join(tmpdir.name, "atac_fragments.tsv.gz")

data_url = "https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_500_nextgem/atac_pbmc_500_nextgem_fragments.tsv.gz"
subprocess.run(["curl", data_url], stdout=open(fragments_10x_path, "w"))
```

```{code-cell} ipython3
metadata_url = "https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_500_nextgem/atac_pbmc_500_nextgem_singlecell.csv"
metadata_path = os.path.join(tmpdir.name, "cell_metadata.csv")
subprocess.run(["curl", metadata_url], stdout=open(metadata_path, "w"))

cell_metadata = pd.read_csv(metadata_path)
cell_metadata = cell_metadata[cell_metadata.is__cell_barcode == 1].reset_index()
cell_metadata
```

```{code-cell} ipython3
cell_metadata.is__cell_barcode.sum()
```

## Convert to BPCells format
Notice that the conversion allows for adjusting the start/end coordinates, as well as subsetting to only the barcodes passing QC. Adding 1 to the end coordinate is necessary for 10x inputs produced by cellranger

```{code-cell} ipython3
%%time
fragments_bpcells_path = os.path.join(tmpdir.name, "bpcells_fragments")
bpcells.import_10x_fragments(
    input = fragments_10x_path, 
    output = fragments_bpcells_path, 
    shift_end=1, 
    keeper_cells=cell_metadata.barcode[cell_metadata.is__cell_barcode == 1]
)
```

## Create the insertion matrix
To calculate the insertion matrix, we first define our cell groups, as well as the ordering of the cell groups we want for our output matrix. Here we use the first two characters of the cell barcode since annotated cell types are not available. Note that it is possible to leave cells out when calling `build_cell_groups`, in which case that data will not be included in the precalculated matrix

Next, we precalculate the insertion counts matrix, which can use parallelization to speed up portions of the work.

```{code-cell} ipython3
%%time
barcodes = cell_metadata.barcode
clusters = cell_metadata.barcode.str.slice(0,2)
cluster_order = sorted(set(clusters))

cell_groups_array = bpcells.build_cell_groups(fragments_bpcells_path, barcodes, clusters, cluster_order)

# We could provide a dict or local file path, but URL is easier
chrom_sizes = "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes"

insertions_matrix_path = os.path.join(tmpdir.name, "bpcells_insertions_matrix")

bpcells.precalculate_insertion_counts(
    fragments_bpcells_path, 
    insertions_matrix_path, 
    cell_groups_array, 
    chrom_sizes, 
    threads=4
)
```

## Querying the insertion matrix

We can load the pre-calculated matrix from its input path.

```{code-cell} ipython3
mat = bpcells.PrecalculatedInsertionMatrix(insertions_matrix_path)
mat
```

```{code-cell} ipython3
mat.shape
```

To query the matrix, we use a pandas DataFrame, with columns (chrom, start, end). All the regions must be the same length

```{code-cell} ipython3
query_regions = pd.DataFrame({
    "chrom": ["chr1", "chr1", "chr6"],
    "start": [1_000_000, 2_000_000, 10_000_000],
})
query_regions["end"] = query_regions.start + 1000
query_regions
```

BPCells returns a numpy array of dimensions (regions, cell types, basepairs), holding the per-base counts for each cell type

```{code-cell} ipython3
x = mat.get_counts(query_regions)
x
```

```{code-cell} ipython3
x.shape
```

```{code-cell} ipython3
x.sum()
```

## Pytorch-compatible dataset

It is simple to wrap this matrix as a pytorch-compatible dataset, given a set of regions for the training set. Note the use of the non-standard `__getitems__()` function which pytorch uses to provide batched loading for higher performance. This dataset object can be directly passed to `torch.utils.data.DataLoader`.

```{code-cell} ipython3
class BPCellsDataset:
    def __init__(self, regions, matrix_dir):
        self.regions = regions[["chrom", "start", "end"]]

        matrix_dir = str(os.path.abspath(os.path.expanduser(matrix_dir)))
        self.mat = bpcells.PrecalculatedInsertionMatrix(matrix_dir)
        
        peak_width = self.regions.end[0] - self.regions.start[0]
        assert (self.regions.end - self.regions.start == peak_width).all()

    def __getitem__(self, i):
        return self.__getitems__([i])[0]

    def __getitems__(self, idx):
        # Adding this function allows for batched loading
        # See: https://github.com/pytorch/pytorch/issues/107218

        # Return tensor of shape (batch_size, n_tasks, basepairs)
        return self.mat.get_counts(
            self.regions.iloc[idx,]
        )

    def __len__(self):
        return self.regions.shape[0]
```
