# BPCells

BPCells is a package for high performance single cell analysis on RNA-seq and ATAC-seq datasets. It can analyze
a 1.3M cell dataset with 2GB of RAM in under 10 minutes. This makes analysis of million-cell datasets practical on a laptop.

BPCells provides:

  - Efficient storage of single cell datasets via bitpacking compression
  - Fast, disk-backed RNA-seq and ATAC-seq data processing powered by C++
  - Downstream analysis such as marker genes, and clustering
  - Interoperability with AnnData, 10x datasets, R sparse matrices, and GRanges

Additionally, BPCells exposes its optimized data processing infrastructure for use in scaling 3rd party single cell tools (e.g. Seurat)

## [Learn more at our website](https://bnprks.github.io/BPCells/)

- [Benchmarks](https://bnprks.github.io/BPCells/articles/web-only/benchmarks.html)
- [Multiomic analysis example](https://bnprks.github.io/BPCells/articles/pbmc3k.html)
- [How BPCells works](https://bnprks.github.io/BPCells/articles/web-only/how-it-works.html)
- [Function documentation](https://bnprks.github.io/BPCells/reference/index.html)
- [News](https://bnprks.github.io/BPCells/news/index.html)

## Installation
BPCells is easiest to install directly from github:

```R
remotes::install_github("bnprks/BPCells")
```
Before installing, you must have the HDF5 library installed and accessible on your system.
HDF5 can be installed from your choice of package manager:

- conda: `conda install -c anaconda hdf5` 
- apt: `sudo apt-get install libhdf5-dev` 
- yum: `yum install hdf5-devel`

## Contributing
BPCells is an open source project, and we welcome quality contributions. If you
are interested in contributing and have experience with C++, along with Python
or R, feel free to reach out with ideas you would like to implement yourself.
I'm happy to provide pointers for how to get started, my time permitting.

If you are unfamiliar with C++ it will be difficult for you to contribute code,
but detailed bug reports with
[reproducible examples](https://reprex.tidyverse.org/articles/reprex-dos-and-donts.html)
are still a useful way to help out. Github issues are the best forum for this.

If you maintain a single cell analysis package and want to use BPCells to
improve your scalability, I'm happy to provide advice. We have had a couple of labs
try this so far, with promising success. Email is the best way to get in touch
for this (look in the `DESCRIPTION` file on github for contact info). Python
developers welcome, though the full python package will likely not be
available until mid-summer 2023.

AnnData maintainers: would love to talk about putting bitpacking compression in
AnnData. The [benchmarks](https://bnprks.github.io/BPCells/articles/web-only/benchmarks.html#counts-matrices-rna-or-atac) look promising.