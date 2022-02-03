devtools::load_all("~/github/bnprks/BPCells")
# Test writing to matrix files on ARM

cat("read_10x_gene_h5\n")
# system.time({
#   mat <- read_10x_gene_h5("~/data/pbmc_granulocyte_sorted_10k_raw_feature_bc_matrix.h5")
# })
mat <- readRDS("~/data/pbmc_10k_rna_mat.rds")
dimnames(mat) <- NULL

# Time: 1.9s User, 60Mb
cat("mem_p\n")
system.time({
  mem_p <- write_matrix_memory(mat)
})
all.equal(mat, as(mem_p, "dgCMatrix"))

cat("mem_u\n")
system.time({
  mem_u <- write_matrix_memory(mem_p, compress=FALSE)
})
all.equal(mat, as(mem_u, "dgCMatrix"))

# Time 0.58s user, 0.917 elapsed, 63.8Mb
cat("h5_p\n")
system.time({
  h5_p <- write_matrix_hdf5(mat, "~/data/arm/matrix_test_packed.h5", "packed")
})
all.equal(mat, as(h5_p, "dgCMatrix"))

# Time 0.95s user, 2.237s elapsed, 251.6Mb
cat("h5_u\n")
system.time({
  h5_u <- write_matrix_hdf5(mat, "~/data/arm/matrix_test_unpacked.h5", "unpacked", compress = FALSE)
})
all.equal(mat, as(h5_u, "dgCMatrix"))

# Time 3.92s user, 4.24s elapsed, 61Mb
# If I don't turn off buffering it goes in 0.68s elapsed
cat("d_p\n")
system.time({
  d_p <- write_matrix_dir(mat, "~/data/arm/matrix_test_packed", buffer_size = 8192L)
})
all.equal(mat, as(d_p, "dgCMatrix"))

# Time 32s user, 34s elapsed, 260Mb
# If I don't turn off buffering it goes in 1.2s elapsed
cat("d_u\n")
system.time({
  d_u <- write_matrix_dir(mat, "~/data/arm/matrix_test_unpacked", compress=FALSE)
})
all.equal(mat, as(d_u, "dgCMatrix"))

bench::mark(
  colSums(mat),
  colSums(h5_p),
  colSums(h5_u),
  colSums(d_p),
  colSums(d_u),
  colSums(mem_u),
  colSums(mem_p)
)