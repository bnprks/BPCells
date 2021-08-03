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


tf_idf_manual <- function(x, scaleTo=1e4) {
    x <- t( t(x/rowMeans(x)) /colSums(x) )
    log(1 + scaleTo * x)
}

test_that("TFIDF normalization works", { 
    withr::local_seed(195123)
    t <- Matrix::t

    m1 <- generate_sparse_matrix(100, 200)
    m2 <- t(m1)

    i1 <- as(m1, "IterableMatrix")
    i2 <- t(i1)

    i1_tf <- normalize_TFIDF(i1)
    i2_tf <- normalize_TFIDF(i2, scaleTo=5)

    m1_tf <- tf_idf_manual(as.matrix(m1))
    m2_tf <- tf_idf_manual(as.matrix(m2), scaleTo=5)

    b <- generate_dense_matrix(200, 3)
    b2 <- generate_dense_matrix(100, 3)

    expect_equal(i1_tf %*% b, m1_tf %*% b)
    expect_equal(t(b) %*% i2_tf, t(b) %*% m2_tf)

    expect_equal(t(b2) %*% i1_tf, t(b2) %*% m1_tf)
    expect_equal(i2_tf %*% b2, m2_tf %*% b2) 
})