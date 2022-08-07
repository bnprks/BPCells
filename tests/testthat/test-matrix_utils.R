library(Matrix)

generate_sparse_matrix <- function(nrow, ncol, fraction_nonzero=0.5, max_val=10) {
    m <- matrix(rbinom(nrow*ncol, 1, fraction_nonzero)*sample.int(max_val, nrow*ncol, replace=TRUE), nrow=nrow)
    as(m, "dgCMatrix")
}
generate_dense_matrix <- function(nrow, ncol) {
    m <- matrix(runif(nrow*ncol), nrow=nrow)
}

to_matrix <- function(x) {
    x <- as.matrix(x)
    attr(x, "dimnames") <- NULL
    x
}
to_vector <- function(x) {
    x <- as.numeric(x)
    attributes(x) <- NULL
    x
}

test_that("Chained subsetting works", {
    m1 <- generate_dense_matrix(10, 10) %>% as("dgCMatrix")

    m2 <- as(m1, "IterableMatrix")

    # Singleton selections
    expect_equal(m1[c(2,4,3,5),], m2[c(2,4,3,5),] %>% as("dgCMatrix"))
    expect_equal(m1[,c(2,4,3,5)], m2[,c(2,4,3,5)] %>% as("dgCMatrix"))

    # Repeated selections on 1 axis
    expect_equal(m1[c(2,4,3,5),][c(3,1),], m2[c(2,4,3,5),][c(3,1),] %>% as("dgCMatrix"))
    expect_equal(m1[,c(2,4,3,5)][,c(3,1)], m2[,c(2,4,3,5)][,c(3,1)] %>% as("dgCMatrix"))
    
    # Repeated selections on 2 axis
    expect_equal(m1[c(2,4,3,5),c(5,4,2,3)][c(3,1),c(1,3)], m2[c(2,4,3,5),c(5,4,2,3)][c(3,1),c(1,3)] %>% as("dgCMatrix"))

    # Missing one axis in 2nd selection
    expect_equal(m1[c(2,4,3,5),c(5,4,2,3)][c(3,1),], m2[c(2,4,3,5),c(5,4,2,3)][c(3,1),] %>% as("dgCMatrix"))
    expect_equal(m1[c(2,4,3,5),c(5,4,2,3)][,c(1,3)], m2[c(2,4,3,5),c(5,4,2,3)][,c(1,3)] %>% as("dgCMatrix"))

    # Missing one axis in 1st selection
    expect_equal(m1[,c(5,4,2,3)][c(3,1),c(1,3)], m2[,c(5,4,2,3)][c(3,1),c(1,3)] %>% as("dgCMatrix"))
    expect_equal(m1[c(2,4,3,5),][c(3,1),c(1,3)], m2[c(2,4,3,5),][c(3,1),c(1,3)] %>% as("dgCMatrix"))
})

test_that("Dense matrix/vector multiply works", { 
    withr::local_seed(195123)

    m1 <- generate_sparse_matrix(5, 1000)
    m2 <- t(m1)

    i1 <- as(m1, "IterableMatrix")
    i2 <- t(i1)

    b <- generate_dense_matrix(5, 12)

    expect_equal(to_matrix(t(b) %*% m1), t(b) %*% i1)
    expect_equal(to_matrix(m2 %*% b), i2 %*% b)

    expect_equal(to_matrix(t(m1) %*% b), t(i1) %*% b)
    expect_equal(to_matrix(t(b) %*% t(m2)) , t(b) %*% t(i2))    

    y <- as.numeric(generate_dense_matrix(5, 1))
    expect_equal(to_vector(y %*% m1), y %*% i1)
    expect_equal(to_vector(m2 %*% y), i2 %*% y)

    expect_equal(to_vector(t(m1) %*% y), t(i1) %*% y)
    expect_equal(to_vector(y %*% t(m2)) , y %*% t(i2))   

    # Check that everything works fine with integers
    i3 <- convert_matrix_type(i1, "uint32_t")
    i4 <- convert_matrix_type(i2, "uint32_t")

    expect_equal(to_matrix(t(b) %*% m1), t(b) %*% i3)
    expect_equal(to_matrix(m2 %*% b), i4 %*% b)

    expect_equal(to_matrix(t(m1) %*% b), t(i3) %*% b)
    expect_equal(to_matrix(t(b) %*% t(m2)) , t(b) %*% t(i4)) 

    expect_equal(to_vector(y %*% m1), y %*% i3)
    expect_equal(to_vector(m2 %*% y), i4 %*% y)

    expect_equal(to_vector(t(m1) %*% y), t(i3) %*% y)
    expect_equal(to_vector(y %*% t(m2)) , y %*% t(i4))      

})

test_that("Row/Col sum/mean works", {
    m1 <- generate_sparse_matrix(5, 1000)
    m2 <- t(m1)
    
    i1 <- as(m1, "IterableMatrix")
    i2 <- t(i1)

    expect_equal(rowSums(i1), rowSums(m1))
    expect_equal(rowSums(i2), rowSums(m2))
    expect_equal(colSums(i1), colSums(m1))
    expect_equal(colSums(i2), colSums(m2))

    expect_equal(rowMeans(i1), rowMeans(m1))
    expect_equal(rowMeans(i2), rowMeans(m2))
    expect_equal(colMeans(i1), colMeans(m1))
    expect_equal(colMeans(i2), colMeans(m2))

    # Check that everything works fine with integers
    i3 <- convert_matrix_type(i1, "uint32_t")
    i4 <- convert_matrix_type(i2, "uint32_t")

    expect_equal(rowSums(i3), rowSums(m1))
    expect_equal(rowSums(i4), rowSums(m2))
    expect_equal(colSums(i3), colSums(m1))
    expect_equal(colSums(i4), colSums(m2))

    expect_equal(rowMeans(i3), rowMeans(m1))
    expect_equal(rowMeans(i4), rowMeans(m2))
    expect_equal(colMeans(i3), colMeans(m1))
    expect_equal(colMeans(i4), colMeans(m2))
})


test_that("Generic methods work", {
    # Generic methods to test:
    # - dim
    # - dimnames
    # - matrix_type
    # - matrix_is_transform
    # - show
    # - rowSums, colSums, 
    m <- generate_sparse_matrix(5, 1000)
    rownames(m) <- paste0("row", seq_len(nrow(m)))
    colnames(m) <- paste0("col", seq_len(ncol(m)))
    mi <- as(m, "IterableMatrix")
    
    dir <- withr::local_tempdir()

    id_right <- as(Matrix::Diagonal(ncol(m)), "dgCMatrix")
    colnames(id_right) <- colnames(m)
    id_left <- as(Matrix::Diagonal(nrow(m)), "dgCMatrix")
    rownames(id_left) <- rownames(m)

    ident_transforms <- list(
        write_memory_uint32_t = mi %>% convert_matrix_type("uint32_t") %>% write_matrix_memory(compress=TRUE),
        write_memory_unpacked_uint32_t = mi %>% convert_matrix_type("uint32_t") %>% write_matrix_memory(compress=FALSE),
        write_memory_float = mi %>% convert_matrix_type("float") %>% write_matrix_memory(compress=TRUE),
        write_memory_unpacked_float = mi %>% convert_matrix_type("float") %>% write_matrix_memory(compress=FALSE),
        write_memory_double = mi %>% convert_matrix_type("double") %>% write_matrix_memory(compress=TRUE),
        write_memory_unpacked_double = mi %>% convert_matrix_type("double") %>% write_matrix_memory(compress=FALSE),
                
        shift_scale_1 = {t((mi * 1 + 0) * rep_len(1, nrow(m))) * rep_len(1, ncol(m)) + rep_len(0, ncol(m))} %>% t(),
        shift_scale_2 = {t((mi / 1 - 0) / rep_len(1, nrow(m)) - rep_len(0, nrow(m))) / rep_len(1, ncol(m))} %>% t(),

        multiply_right_1 = mi %*% id_right,
        multiply_right_2 = mi %*% as(id_right, "IterableMatrix"),
        multiply_left_1 = id_left %*% mi,
        multiply_types = convert_matrix_type(mi, "uint32_t") %*% convert_matrix_type(as(id_right, "IterableMatrix"), "uint32_t"),
        multiply_types = convert_matrix_type(mi, "float") %*% convert_matrix_type(as(id_right, "IterableMatrix"), "uint32_t"),
        
        subset = mi[seq_len(nrow(m)),][,seq_len(ncol(m))][seq_len(nrow(m)), seq_len(ncol(m))],
        rbind = rbind2(mi[1:2,], mi[3:nrow(m)]),
        cbind = cbind2(mi[,1:2], mi[,3:ncol(m)])
    )

    for (i in names(ident_transforms)) {
        t <- ident_transforms[[i]]
        short_description(t)
        expect_equal(dim(t), dim(m))
        expect_equal(dimnames(t), dimnames(m))
        expect_true(matrix_type(t) %in% c("uint32_t", "float", "double"))
        matrix_is_transform(t)
        
        expect_equal(rowSums(t), rowSums(m))
        expect_equal(colSums(t), colSums(m))
        
        if (i %in% c("shift_scale_1", "shift_scale_2"))
            expect_equal(as.matrix(as(t, "dgCMatrix")), as.matrix(m))
        else
            expect_equal(as(t, "dgCMatrix"), m)
    }
})
