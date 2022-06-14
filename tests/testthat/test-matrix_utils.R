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
    library(Matrix)
    t <- Matrix::t
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
    i3 <- cast_matrix_int(i1)
    i4 <- cast_matrix_int(i2)

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
    library(Matrix)
    t <- Matrix::t
    
    m1 <- generate_sparse_matrix(5, 1000)
    m2 <- Matrix::t(m1)
    
    i1 <- as(m1, "IterableMatrix")
    i2 <- Matrix::t(i1)

    expect_equal(rowSums(i1), rowSums(m1))
    expect_equal(rowSums(i2), rowSums(m2))
    expect_equal(colSums(i1), colSums(m1))
    expect_equal(colSums(i2), colSums(m2))

    expect_equal(rowMeans(i1), rowMeans(m1))
    expect_equal(rowMeans(i2), rowMeans(m2))
    expect_equal(colMeans(i1), colMeans(m1))
    expect_equal(colMeans(i2), colMeans(m2))

    # Check that everything works fine with integers
    i3 <- cast_matrix_int(i1)
    i4 <- cast_matrix_int(i2)

    expect_equal(rowSums(i3), rowSums(m1))
    expect_equal(rowSums(i4), rowSums(m2))
    expect_equal(colSums(i3), colSums(m1))
    expect_equal(colSums(i4), colSums(m2))

    expect_equal(rowMeans(i3), rowMeans(m1))
    expect_equal(rowMeans(i4), rowMeans(m2))
    expect_equal(colMeans(i3), colMeans(m1))
    expect_equal(colMeans(i4), colMeans(m2))
})
