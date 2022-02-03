devtools::load_all("~/github/bnprks/BPCells")
# Test reading from matrix files on ARM

cat("read_10x_gene_h5\n")
mat <- readRDS("~/data/pbmc_10k_rna_mat.rds")
dimnames(mat) <- NULL

# Time 0.58s user, 0.917 elapsed, 63.8Mb
cat("h5_p\n")
h5_p <- open_matrix_hdf5("~/data/arm/matrix_test_packed.h5", "packed")
all.equal(mat, as(h5_p, "dgCMatrix"))

# Time 0.95s user, 2.237s elapsed, 251.6Mb
cat("h5_u\n")
h5_u <- open_matrix_hdf5("~/data/arm/matrix_test_unpacked.h5", "unpacked")
all.equal(mat, as(h5_u, "dgCMatrix"))

# Time 3.92s user, 4.24s elapsed, 61Mb
# If I don't turn off buffering it goes in 0.68s elapsed
cat("d_p\n")
d_p <- open_matrix_dir("~/data/arm/matrix_test_packed")
all.equal(mat, as(d_p, "dgCMatrix"))

# Time 32s user, 34s elapsed, 260Mb
# If I don't turn off buffering it goes in 1.2s elapsed
cat("d_u\n")
d_u <- open_matrix_dir("~/data/arm/matrix_test_unpacked")
all.equal(mat, as(d_u, "dgCMatrix"))
