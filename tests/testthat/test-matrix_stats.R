
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


test_that("MatrixStats basic test", { 
    withr::local_seed(195123)

    m1 <- generate_sparse_matrix(5, 1000)
    m2 <- t(m1)

    i1 <- as(m1, "IterableMatrix")
    i2 <- t(i1)

    
    stats1 <- matrixStats(i1, "variance", "variance")
    stats2 <- matrixStats(i2, "variance", "variance")

    expect_equal(c("nonzero", "mean", "variance"), rownames(stats1$row_stats))
    expect_equal(c("nonzero", "mean", "variance"), rownames(stats1$col_stats))

    expect_equal(c("nonzero", "mean", "variance"), rownames(stats2$row_stats))
    expect_equal(c("nonzero", "mean", "variance"), rownames(stats2$col_stats))

    expect_equal(stats1$row_stats["nonzero",], rowSums(as.matrix(m1 != 0)))
    expect_equal(stats1$row_stats["mean",], rowMeans(as.matrix(m1)))
    
    expect_equal(stats1$row_stats["variance",], matrixStats::rowVars(as.matrix(m1)))
    expect_equal(stats1$col_stats["nonzero",], colSums(as.matrix(m1 != 0)))
    
    expect_equal(stats1$col_stats["mean",], colMeans(as.matrix(m1)))
    expect_equal(stats1$col_stats["variance",], matrixStats::colVars(as.matrix(m1)))
})

test_that("MatrixStats comprehensive tests", { 
    withr::local_seed(195123)
    
    m1 <- generate_sparse_matrix(5, 1000)
    m2 <- t(m1)

    i1 <- as(m1, "IterableMatrix")
    i2 <- t(i1)

    for (row_stats in c("none", "nonzero", "mean", "variance")) {
        for (col_stats in c("none", "nonzero", "mean", "variance")) {
            stats1 <- matrixStats(i1, row_stats, col_stats)
            stats2 <- matrixStats(i2, row_stats, col_stats)

            expected_rownames <- c("nonzero", "mean", "variance")[seq_len(match(row_stats, c("none", "nonzero", "mean", "variance"))-1)]
            expected_colnames <- c("nonzero", "mean", "variance")[seq_len(match(col_stats, c("none", "nonzero", "mean", "variance"))-1)]
            
            if(row_stats == "none") expected_rownames <- NULL
            if(col_stats == "none") expected_colnames <- NULL

            expect_equal(expected_rownames, rownames(stats1$row_stats))
            expect_equal(expected_colnames, rownames(stats1$col_stats))

            expect_equal(expected_rownames, rownames(stats2$row_stats))
            expect_equal(expected_colnames, rownames(stats2$col_stats))

            if (row_stats %in% c("nonzero", "mean", "variance")) {
                expect_equal(stats1$row_stats["nonzero",], rowSums(as.matrix(m1 != 0)))
            }
            if (row_stats %in% c("mean", "variance")) {
                expect_equal(stats1$row_stats["mean",], rowMeans(as.matrix(m1)))
            }
            if (row_stats %in% c("variance")) {
                expect_equal(stats1$row_stats["variance",], matrixStats::rowVars(as.matrix(m1)))
            }

            if (col_stats %in% c("nonzero", "mean", "variance")) {
                expect_equal(stats1$col_stats["nonzero",], colSums(as.matrix(m1 != 0)))
            }
            if (col_stats %in% c("mean", "variance")) {
                expect_equal(stats1$col_stats["mean",], colMeans(as.matrix(m1)))
            }
            if (col_stats %in% c("variance")) {
                expect_equal(stats1$col_stats["variance",], matrixStats::colVars(as.matrix(m1)))
            }

        }
    }
})