devtools::load_all("~/Dropbox/greenleaf/playground/fragment_io/BPCells/")

read_10x_gene_h5 <- function(path) {
  h5 <- hdf5r::H5File$new(path, mode="r")
  rna <- Matrix::sparseMatrix(
    i=h5[["matrix/indices"]][]+1, x=h5[["matrix/data"]][], p=h5[["matrix/indptr"]][], dims = h5[["matrix/shape"]][]
  )
  rownames(rna) <- h5[["matrix/features/id"]][]
  colnames(rna) <- h5[["matrix/barcodes"]][]
  rna <- rna["Gene Expression" == h5[["matrix/features/feature_type"]][], ]
  h5$close()
  rna
}

# Time: 33s User, 48s Elapsed, 419Mb
system.time({
  mat <- read_10x_gene_h5("~/Downloads/pbmc_granulocyte_sorted_10k_raw_feature_bc_matrix.h5")
})
mat <- readRDS("~/Downloads/pbmc_10k_rna_mat.rds")
dimnames(mat) <- NULL

# Time: 1.9s User, 60Mb
system.time({
  mem_p <- write_matrix_memory(mat)
})

system.time({
  mem_u <- write_matrix_memory(mem_p, compress=FALSE)
})

# Time 0.58s user, 0.917 elapsed, 63.8Mb
system.time({
  h5_p <- write_matrix_hdf5(mat, "~/Downloads/matrix_test_packed.h5", "packed")
})

# Time 0.95s user, 2.237s elapsed, 251.6Mb
system.time({
  h5_u <- write_matrix_hdf5(mat, "~/Downloads/matrix_test_unpacked.h5", "unpacked", compress = FALSE)
})

# Time 3.92s user, 4.24s elapsed, 61Mb
# If I don't turn off buffering it goes in 0.68s elapsed
system.time({
  d_p <- write_matrix_dir(mat, "~/Downloads/matrix_test_packed", buffer_size = 8192L)
})

# Time 32s user, 34s elapsed, 260Mb
# If I don't turn off buffering it goes in 1.2s elapsed
system.time({
  d_u <- write_matrix_dir(mat, "~/Downloads/matrix_test_unpacked", compress=FALSE)
})

m_10x <- open_matrix_10x_hdf5("~/Downloads/pbmc_granulocyte_sorted_10k_raw_feature_bc_matrix.h5")

bench::mark(
  colSums(mat),
  colSums(h5_p),
  colSums(h5_u),
  colSums(d_p),
  colSums(d_u),
  colSums(mem_u),
  colSums(mem_p)
)


