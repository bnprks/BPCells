# BPCells

BPCells is a package for high performance single cell analysis on RNA-seq and ATAC-seq datasets. It can analyze
a 1.3M cell dataset with 2GB of RAM in around 10 minutes ([benchmarks](articles/web-only/benchmarks.html)). This makes analysis of million-cell datasets practical on a laptop.

BPCells provides:

  - Efficient storage of single cell datasets via bitpacking compression
  - Fast, disk-backed RNA-seq and ATAC-seq data processing powered by C++
  - Downstream analysis such as marker genes, and clustering
  - Interoperability with AnnData, 10x datasets, R sparse matrices, and GRanges

Additionally, BPCells exposes its optimized data processing infrastructure for use in scaling 3rd party single cell tools (e.g. Seurat)

## Learn more

- [Benchmarks](articles/web-only/benchmarks.html)
- [Multiomic analysis example](articles/pbmc3k.html)
- [How BPCells works](articles/web-only/how-it-works.html)
- [Additional articles](articles/index.html)
- [Function documentation](reference/index.html)
- [News](news/index.html)

## Installation
BPCells is easiest to install directly from github:

```R
remotes::install_github("bnprks/BPCells")
```
Before installing, you must have the HDF5 library installed and accessible on your system.
HDF5 can be installed from your choice of package manager:

- apt: `sudo apt-get install libhdf5-dev` 
- yum: `sudo yum install hdf5-devel`
- conda: `conda install -c anaconda hdf5` 
  - Note: Linux users should prefer their distro's package manager (e.g. `apt` or `yum`) when possible,
    as it appears to give a slightly more reliable installation experience.

You will also need a C/C++ compiler either gcc >=8.0 (>=9.1 recommended), or clang >= 7.0 (>= 9.0 recommended).
This corresponds to versions from late-2018 and newer.

**Installation troubleshooting** -- see these github issues:

- Windows ([Issue #9](https://github.com/bnprks/BPCells/issues/9#issuecomment-1489741180))
- M1 Macs ([Issue #6](https://github.com/bnprks/BPCells/issues/6#issuecomment-1476976677))
- Macs running MacOS 10.14 (Mojave) or older ([Issue #3](https://github.com/bnprks/BPCells/issues/3#issuecomment-1375238635))

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