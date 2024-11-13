# BPCells

The BPCells python bindings are still experimental and the API is subject to change.

The existing functionality is mainly focused on allowing read/write access to BPCells
file formats for integer matrices and scATAC fragments. Future updates will add the
data-processing functions present in the R interface (e.g. streaming normalization, PCA,
or ATAC-seq peak/tile matrix creation). This will provide Python access to the shared
C++ core code.

Notably, plotting functionality is not currently planned for implementation, as it is 
written primarily in R and relies on R plotting libraries not present in Python. There
are a few other helper functions in R BPCells that are implemented in pure R and thus
are unlikely to be added in Python in the near future. If any of this functionality
is of interest to you, we would welcome your contributions -- you would be able
to write most of the code in pure Python. Reach out via github/email if interested.


## Installation

BPCells can be directly installed via pip:

```shell
python -m pip install bpcells
```

## Tutorials

- [Matrix slicing](notebooks/matrix_basics)
- [Basepair insertion dataloading](notebooks/fragment_basics)

## API Reference

- [Fragment functions](api/fragments)
- [Matrix functions](api/matrix)


:::{toctree}
:maxdepth: 2

python
R Docs <https://bnprks.github.io/BPCells>
:::