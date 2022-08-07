generate_sparse_matrix <- function(nrow, ncol, fraction_nonzero=0.5, max_val=10) {
    m <- matrix(rbinom(nrow*ncol, 1, fraction_nonzero)*sample.int(max_val, nrow*ncol, replace=TRUE), nrow=nrow)
    as(m, "dgCMatrix")
}


test_that("Memory Matrix round-trips", {
    args <- list(
        list(nrow=10, ncol=10, fraction_nonzero=0.2, max_val=10),
        list(nrow=100, ncol=100, fraction_nonzero=0.2, max_val=100),
        list(nrow=1000, ncol=1000, fraction_nonzero=0.1, max_val=1000)
    )
    for (a in args) {
        m <- do.call(generate_sparse_matrix, a)
        for (type in c("uint32_t", "float", "double")) {
            m2 <- m %>% as("IterableMatrix") %>% convert_matrix_type(type)
            mp <- write_matrix_memory(m2, compress = TRUE)
            mu <- write_matrix_memory(m2, compress = FALSE)
            expect_equal(m, as(mp, "dgCMatrix"))
            expect_equal(m, as(mu, "dgCMatrix"))
            expect_equal(matrix_type(mp), type)
            expect_equal(matrix_type(mu), type)
            expect_equal(mu, write_matrix_memory(mp, compress=FALSE))
        }
    }
})

test_that("DirPacked Matrix round-trips", {
    dir <- withr::local_tempdir()
    args <- list(
        list(nrow=10, ncol=10, fraction_nonzero=0.2, max_val=10),
        list(nrow=100, ncol=100, fraction_nonzero=0.2, max_val=100),
        list(nrow=1000, ncol=1000, fraction_nonzero=0.1, max_val=1000)
    )
    for (a in args) {
        m <- do.call(generate_sparse_matrix, a)
        for (type in c("uint32_t", "float", "double")) {
            m2 <- m %>% as("IterableMatrix") %>% convert_matrix_type(type)
            mem <- write_matrix_memory(m2, compress = FALSE)
            mp <- write_matrix_dir(m2, file.path(dir, "packed", type, a[["nrow"]]), compress = TRUE)
            mu <- write_matrix_dir(m2, file.path(dir, "unpacked", type, a[["nrow"]]), compress = FALSE)
            expect_equal(m, as(mp, "dgCMatrix"))
            expect_equal(m, as(mu, "dgCMatrix"))
            expect_equal(matrix_type(mp), type)
            expect_equal(matrix_type(mu), type)
            expect_equal(mem, write_matrix_memory(mp, compress = FALSE))
            expect_equal(mem, write_matrix_memory(mu, compress = FALSE))
        }
    }
})


test_that("H5Packed Matrix round-trips", {
    dir <- withr::local_tempdir()
    args <- list(
        list(nrow=10, ncol=10, fraction_nonzero=0.2, max_val=10),
        list(nrow=100, ncol=100, fraction_nonzero=0.2, max_val=100),
        list(nrow=1000, ncol=1000, fraction_nonzero=0.1, max_val=1000)
    )
    for (a in args) {
        m <- do.call(generate_sparse_matrix, a) 
        for (type in c("uint32_t", "float", "double")) {
            m2 <- m %>% as("IterableMatrix") %>% convert_matrix_type(type)
            mem <- write_matrix_memory(m2, compress = FALSE)
            mu <- write_matrix_hdf5(m2, file.path(dir, "subdir", type, "file.h5"), paste0("packed/", as.character(a[["nrow"]])), compress=FALSE)
            mp <- write_matrix_hdf5(m2, file.path(dir, "subdir", type, "file.h5"), paste0("unpacked/", as.character(a[["nrow"]])), compress=TRUE)
            
            expect_equal(m, as(mp, "dgCMatrix"))
            expect_equal(m, as(mu, "dgCMatrix"))
            expect_equal(matrix_type(mp), type)
            expect_equal(matrix_type(mu), type)
            expect_equal(mem, write_matrix_memory(mp, compress = FALSE))
            expect_equal(mem, write_matrix_memory(mu, compress = FALSE))
        }
    }
})

test_that("Packed Matrix works on all bit widths", {
    # Test matrix design
    # - 512 columns, each column has as many entries in it as the column index
    # - Attempt to use col_idx % 32 bits for storing values and col_idx % 8 bits for row offsets
    #   (I had slowness and crashing issues with dgCMatrices that had too many rows)
    # - Test packing the matrix with columns in-order and with shuffled columns

    withr::local_seed(1258123)

    x <- integer(0)
    row <- integer(0)
    col_ptr <- 0L
    for (col in 0:512) {
        max_bits <- col %% 32
        vals <- sample.int(2^max_bits, col, replace=TRUE)
        rows <- seq.int(0, col*2^(col%%8), 2^(col%%8))[seq_len(col)] + 1

        col_ptr <- c(col_ptr, col)
        row <- c(row, rows)
        x <- c(x, vals)
    }
    col_ptr <- cumsum(col_ptr)
    m <- Matrix::sparseMatrix(i=row, x=x, p=col_ptr) 

    mp <- write_matrix_memory(convert_matrix_type(m, "uint32_t"))
    
    expect_equal(m, as(mp, "dgCMatrix"))

    col_reorder <- sample.int(ncol(m))
    m <- m[,col_reorder]
    mp <- write_matrix_memory(convert_matrix_type(m, "uint32_t"))
    expect_equal(m, as(mp, "dgCMatrix"))
})

test_that("Saving matrix tranpose works", {
    dir <- withr::local_tempdir()
    args <- list(
        list(nrow=10, ncol=10, fraction_nonzero=0.2, max_val=10),
        list(nrow=100, ncol=100, fraction_nonzero=0.2, max_val=100),
        list(nrow=1000, ncol=1000, fraction_nonzero=0.1, max_val=1000)
    )
    for (a in args) {
        m <- do.call(generate_sparse_matrix, a) 
        for (type in c("uint32_t", "float", "double")) {
            m2 <- m %>% as("IterableMatrix") %>% convert_matrix_type(type)
            t_mem <- t(m2) %>% write_matrix_memory()
            t_dir <- t(m2) %>% write_matrix_dir(file.path(dir, "packed", type, a[["nrow"]]), compress = TRUE)
            t_hdf5 <- t(m2) %>% write_matrix_hdf5(file.path(dir, "subdir", type, "file.h5"), paste0("unpacked/", as.character(a[["nrow"]])), compress=TRUE)

            t_m <- t(m)
            
            expect_equal(t_m, as(t_mem, "dgCMatrix"))
            expect_equal(t_m, as(t_dir, "dgCMatrix"))
            expect_equal(t_m, as(t_hdf5, "dgCMatrix"))
        }
    }

})