devtools::load_all("~/Dropbox/greenleaf/playground/fragment_io/BPCells/")

input <- open_10x_fragments("~/Dropbox/greenleaf/playground/fragment_io/BPCells/tests/data/mini_fragments.tsv.gz")

#m_u <- write_fragments_memory2(input, compress = FALSE)
m_p <- write_fragments_memory2(input)
as(m_p, "GRanges")


# read_10x_gene_h5 <- function(path) {
#   h5 <- hdf5r::H5File$new(path, mode="r")
#   rna <- Matrix::sparseMatrix(
#     i=h5[["matrix/indices"]][]+1, x=h5[["matrix/data"]][], p=h5[["matrix/indptr"]][], dims = h5[["matrix/shape"]][]
#   )
#   rownames(rna) <- h5[["matrix/features/id"]][]
#   colnames(rna) <- h5[["matrix/barcodes"]][]
#   rna <- rna["Gene Expression" == h5[["matrix/features/feature_type"]][], ]
#   h5$close()
#   rna
# }

# Time: 33s User, 48s Elapsed, 419Mb
# system.time({
#   mat <- read_10x_gene_h5("~/Downloads/pbmc_granulocyte_sorted_10k_raw_feature_bc_matrix.h5")
# })
# system.time({
#     mat <- readRDS("~/Downloads/pbmc_10k_rna_mat.rds")
# })

# system.time({
#   d_p <- write_matrix_dir(mat, "~/Downloads/matrix_test_packed")
# })
# Time: 1.9s User, 60Mb
# system.time({
#   mem_p <- write_matrix_memory(mat)
# })

# system.time({
#   mem_u <- write_matrix_memory(mem_p, compress=FALSE)
# })
# cat("Checking equality\n")
# all.equal(mat, as(mem_p, "dgCMatrix"))
