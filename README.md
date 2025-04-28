# BPCells

BPCells is a package for high performance single cell analysis of large RNA-seq and ATAC-seq datasets. 
It can run normalization and PCA of a 1.3M cell dataset in 4 minutes with 2GB of RAM, or create scATAC-seq peak matrices from fragment coordinates with 50x less CPU time than ArchR or SnapATAC2. 
BPCells can even handle the full CELLxGENE census human dataset, running full precision PCA on a 44M cell x 60k gene matrix in 6 hours on a laptop or <1 hour on a server. See our [benchmarks page](https://bnprks.github.io/BPCells/articles/web-only/benchmarks.html) for details.

BPCells provides:

  - Efficient storage of single cell datasets via bitpacking compression
  - Fast, disk-backed RNA-seq and ATAC-seq data processing powered by C++
  - Downstream analysis such as marker genes, and clustering
  - Interoperability with AnnData, 10x datasets, R sparse matrices, and GRanges
  - Demonstrated scalability to 44M cells on a laptop

Additionally, BPCells exposes its optimized data processing infrastructure for use in scaling 3rd party single cell tools (e.g. Seurat)

## [Learn more at our website](https://bnprks.github.io/BPCells/)

- [BioRxiv preprint](https://www.biorxiv.org/content/10.1101/2025.03.27.645853v1)
- [Python docs](https://bnprks.github.io/BPCells/python/index.html)
- [Benchmarks](https://bnprks.github.io/BPCells/articles/web-only/benchmarks.html)
- [Multiomic analysis example](https://bnprks.github.io/BPCells/articles/pbmc3k.html)
- [How BPCells works](https://bnprks.github.io/BPCells/articles/web-only/how-it-works.html)
- [Additional articles](https://bnprks.github.io/BPCells/articles/index.html)
- [Function documentation](https://bnprks.github.io/BPCells/reference/index.html)
- [News](https://bnprks.github.io/BPCells/news/index.html)

## R Installation
We recommend installing BPCells directly from github:

```R
remotes::install_github("bnprks/BPCells/r")
```
Before installing, you must have the HDF5 library installed and accessible on your system.
HDF5 can be installed from your choice of package manager. See operating system specific instructions below.


For Mac and Windows users having trouble installing from github, check our [R-universe](https://bnprks.r-universe.dev/BPCells)
page for instructions to install pre-built binary packages. These binary packages automatically track the latest github main branch.

BPCells is available via conda thanks to [@mfansler](https://github.com/mfansler) from the Conda Forge R team (see [issue #241](https://github.com/bnprks/BPCells/issues/241) for details).
Any issues with the bioconda package should be reported at [bioconda-recipes](https://github.com/bioconda/bioconda-recipes/). Version updates are managed by the bioconda team.

<details>
<summary>Click here for operating system specific installation information for github-based installs</summary>
<div>

### Linux
Obtaining the HDF5 dependency is usually pretty straightforward on Linux

- apt: `sudo apt-get install libhdf5-dev` 
- yum: `sudo yum install hdf5-devel`
- conda: `conda install -c conda-forge hdf5` 
  - Note: Linux users should prefer their distro's package manager (e.g. `apt` or `yum`) when possible,
    as it appears to give a slightly more reliable installation experience.

### Windows
Compiling R packages from source on Windows requires installing [R tools for Windows](https://cran.r-project.org/bin/windows/Rtools/). See [Issue #9](https://github.com/bnprks/BPCells/issues/9) for more discussion.

### MacOS
For MacOS, installing HDF5 through homebrew seems to be most reliable: `brew install hdf5`.

**Mac-specific troubleshooting**:

- **Macs with ARM CPUs**: a common error is to have an ARM-based HDF5 install but an x86-based 
  R install. This will cause errors when BPCells tries to access HDF5 during installation. 
    - Check your R installation
  by running `sessionInfo()`, and seeing if it lists ARM or x86 under "Platform". 
    - The easiest option is to use
  ARM R because homebrew will default to an ARM hdf5 installation
    - It is [possible](https://codetinkering.com/switch-homebrew-arm-x86/) (though tricky) to install an x86 copy of homebrew in order to access an x86 version of hdf5
- **Older Macs (10.14 Mojave or older)**: The default compiler on old Macs does not support needed
  C++17 filesystem features. See [issue #3](https://github.com/bnprks/BPCells/issues/3#issuecomment-1375238635) for
  tips getting a newer compiler set up via homebrew.

### Supported compilers
In most cases, you will already have an appropriate compiler. BPCells recommends
gcc >=9.1, or clang >= 9.0. 
This corresponds to versions from late-2018 and newer. 
Older versions may work in some cases so long as they
have basic C++17 support, but they are not officially supported.

</div>
</details>
<details>
<summary>Click here for troubleshooting information for github-based installs</summary>

### General Installation troubleshooting
BPCells tries to print informative error messages during compilation to help diagnose the problem. For a more
verbose set of information, run `Sys.setenv(BPCELLS_DEBUG_INSTALL="true")` prior to `remotes::install_github("bnprks/BPCells/r")`. If you still can't solve the issue with that additional information, feel free to file a Github issue, being
sure to use a [collapsible section](https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/organizing-information-with-collapsed-sections) for the verbose installation log.

</details>

## Python Installation

BPCells can be directly installed via pip:

```shell
python -m pip install bpcells
```

## Contributing
BPCells is an open source project, and we welcome quality contributions. If you
are interested in contributing and have experience with C++, along with Python
or R, feel free to reach out with ideas you would like to implement yourself.
We're happy to provide pointers for how to get started, time permitting.

If you are unfamiliar with C++ it will be difficult for you to contribute code,
but detailed bug reports with
[reproducible examples](https://reprex.tidyverse.org/articles/reprex-dos-and-donts.html)
are still a great way to help out. Github issues are the best forum for this.

If you maintain a single cell analysis package and want to use BPCells to
improve your scalability, we're happy to provide advice. We have had a couple of labs
try this so far, with promising success. Email is the best way to get in touch
for this (look in the `DESCRIPTION` file on github for contact info). Python
developers welcome, though the current python package is still in experimental status.