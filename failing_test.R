devtools::load_all("~/Dropbox/greenleaf/playground/fragment_io/BPCells/")

withr::local_seed(1258123)

generate_sparse_matrix <- function(nrow, ncol, fraction_nonzero=0.5, max_val=10) {
    m <- matrix(rbinom(nrow*ncol, 1, fraction_nonzero)*sample.int(max_val, nrow*ncol, replace=TRUE), nrow=nrow)
    as(m, "dgCMatrix")
}

dir <- withr::local_tempdir()
args <- list(
    list(nrow=10, ncol=10, fraction_nonzero=0.2, max_val=10),
    list(nrow=100, ncol=100, fraction_nonzero=0.2, max_val=100),
    list(nrow=1000, ncol=1000, fraction_nonzero=0.1, max_val=1000)
)

m <- generate_sparse_matrix(nrow=10, ncol=10, fraction_nonzero=0.2, max_val=10)
mp <- write_matrix_dir(m, file.path(dir, "packed", 10), compress = TRUE)
mem <- write_matrix_memory(m, compress = FALSE)
mem2 <- write_matrix_memory(mp, compress=FALSE)

vec_search <- function(needle, haystack) {
    initial_matches <- which(needle[1] == haystack)
    full_matches <- purrr::map_lgl(initial_matches, function(i) all(needle == haystack[i:(i+length(needle)-1)]))
    initial_matches[full_matches]
}
# x <- integer(0)
# row <- integer(0)
# col_ptr <- 0L
# for (col in 128:129) {
#     max_bits <- col %% 32
#     vals <- sample.int(2^max_bits, col, replace=TRUE)
#     rows <- seq.int(0, col*2^(col%%8), 2^(col%%8))[seq_len(col)] + 1

#     col_ptr <- c(col_ptr, col)
#     row <- c(row, rows)
#     x <- c(x, vals)
# }
# col_ptr <- cumsum(col_ptr)
# m <- Matrix::sparseMatrix(i=row, x=x, p=col_ptr)

# mp <- write_packed_matrix_memory(as(m, "IterableMatrix"))
# expect_equal(m, as(mp, "dgCMatrix"))


set.seed(231084)
m <- generate_sparse_matrix(nrow=1000, ncol=1000, fraction_nonzero=0.1, max_val=1000)
write_packed_matrix_hdf5(as(m, "IterableMatrix"), "test.hdf5", "d2")
mp <- open_packed_matrix_hdf5("test.hdf5", "d2")
mr <- write_unpacked_matrix_cpp(iterate_matrix(mp))

mp2 <- write_packed_matrix_memory(as(m, "IterableMatrix"))
mr2 <- write_unpacked_matrix_cpp(iterate_matrix(mp2))

mp1 <- write_packed_matrix_memory(as(m, "IterableMatrix"))
mp2 <- write_packed_matrix_memory(mp)
all.equal(m@p, m2@col_ptr)

#all.equal(m, as(m2, "dgCMatrix"))

# generate_sparse_matrix <- function(nrow, ncol, fraction_nonzero=0.5, max_val=10) {
#     m <- matrix(rbinom(nrow*ncol, 1, fraction_nonzero)*sample.int(max_val, nrow*ncol, replace=TRUE), nrow=nrow)
#     as(m, "dgCMatrix")
# }
# set.seed(1234543)
# m <- generate_sparse_matrix(5, 6, fraction_nonzero = 0.25)
# m2 <- write_packed_matrix_memory(as(m, "IterableMatrix"))
# all.equal(m, as(m2, "dgCMatrix"))

m <- matrix(10, nrow=500, ncol=2) %>% as("dgCMatrix")
mp <- write_packed_matrix_memory(as(m, "IterableMatrix"))

# m3 <- iterate_packed_matrix_cpp(m2) %>% write_packed_matrix_cpp()
#m2 <- write_packed_matrix_memory(as(m, "IterableMatrix"))
#m3 <- as(m2, "dgCMatrix")