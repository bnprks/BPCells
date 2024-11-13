# BPCells

BPCells is a package for high performance single cell analysis on RNA-seq and ATAC-seq datasets. It can analyze a 1.3M cell dataset with 2GB of RAM in around 10 minutes (benchmarks). This makes analysis of million-cell datasets practical on a laptop.

The main BPCells interface is in R and the python bindings are still in an 
experimental phase, so the API is subject to change.

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

The current BPCells python bindings provide limited access to BPCells functionality. This
mainly focuses on integer matrix slicing / import, and pseudobulk fragment insertion counts.
The API is subject to change somewhat once a full-featured port of the BPCells R 
functionality is completed.

For more information, see our [github](https://github.com/bnprks/BPCells), [python docs](https://bnprks.github.io/BPCells/python), and [R docs](https://bnprks.github.io/BPCells/).
