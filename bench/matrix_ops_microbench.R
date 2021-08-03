library(PackedInsertions)
library(Rcpp)
library(Matrix)
library(bench)
library(tidyverse)

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
# About 50s to run
system.time({
res_mat <- bench::press(
    rows = 1000,
    cols = c(500, 1000, 5000, 10000, 20000),
    out_cols = c(1, 10, 30),
    {
        x <- generate_sparse_matrix(rows, cols)
        i <- as(x, "IterableMatrix")

        xt <- t(x)
        it <- as(xt, "IterableMatrix")

        y <- generate_dense_matrix(cols, out_cols)
        yt <- t(y)

        stopifnot(all.equal(to_matrix(x %*% y), i %*% y))
        stopifnot(all.equal(to_matrix(yt %*% xt), yt %*% it))

        bench::mark(
            x %*% y,
            i %*% y,
            yt %*% xt,
            yt %*% it,
            check=FALSE
        )
    }
)
})

ggplot(res, aes(cols, as.numeric(median), color=expression)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    facet_wrap("out_cols") +
    ggtitle("Time to run")

res %>% 
    group_by(cols, out_cols) %>%
    mutate(ratio = as.numeric(median/min(median))) %>%
    ggplot(aes(cols, ratio, color=expression)) +
    geom_point() +
    scale_x_log10() +
    facet_wrap("out_cols") +
    ggtitle("Ratio to best method")

res %>% 
    group_by(expression, cols) %>%
    mutate(ratio = as.numeric(median/min(median))) %>%
    ggplot(aes(out_cols, ratio, color=expression)) +
    geom_point() +
    scale_x_log10() +
    facet_wrap("cols") +
    ggtitle("Scaling of time across output columns")

res %>% 
    group_by(expression, out_cols) %>%
    mutate(ratio = as.numeric(median/min(median))) %>%
    ggplot(aes(cols, ratio, color=expression)) +
    geom_point() +
    scale_x_log10() +
    facet_wrap("out_cols") +
    ggtitle("Scaling of time across input columns")

# 23 seconds
system.time({
    res_vec <- bench::press(
        rows = 1000,
        cols = c(500, 1000, 5000, 10000, 20000),
        {
            x <- generate_sparse_matrix(rows, cols)
            i <- as(x, "IterableMatrix")
            
            xt <- t(x)
            it <- as(xt, "IterableMatrix")
            
            y <- as.numeric(generate_dense_matrix(cols, 1))
            ym <- matrix(y, nrow=cols)
            ymt <- t(ym)
            
            bench::mark(
                x %*% y,
                i %*% ym,
                i %*% y,
                y %*% xt,
                ymt %*% it,
                y %*% it
                check=function(x,y) {
                    x <- as.numeric(x)
                    y <- as.numeric(y)
                    attributes(x) <- NULL
                    attributes(y) <- NULL
                    all.equal(x,y)
                }
            )
        }
    )
})

ggplot(res_vec, aes(cols, as.numeric(median), color=expression)) +
    geom_point()

res_vec %>% 
    group_by(cols) %>%
    mutate(ratio = as.numeric(median/min(median))) %>%
    ggplot(aes(cols, ratio, color=as.character(expression))) +
    geom_point() +
    scale_x_log10() +
    scale_color_brewer(palette="Set1") +
    ggtitle("Ratio to best method")


# x <- generate_sparse_matrix(1000, 20000)
# i <- as(x, "IterableMatrix")
# y <- generate_dense_matrix(20000, 30)
# start_profiler("right_mul.prof")
# for (i in 1:100) {x %*% y}
# stop_profiler()

# About 27s to run
system.time({
    res_sums <- bench::press(
        rows = 1000,
        cols = c(500, 1000, 5000, 10000, 20000),
        {
            x <- generate_sparse_matrix(rows, cols)
            i <- as(x, "IterableMatrix")
            
            xt <- t(x)
            it <- as(xt, "IterableMatrix")
            
            ones <- rep(1, cols)
            
            bench::mark(
                rowSums(x),
                rowSums(i),
                i %*% ones,
                colSums(xt),
                col_sums_manual(PackedInsertions:::iterate_matrix(it)),
                colSums(it),
                ones %*% i,
            )
        }
    )
})
ggplot(res_sums, aes(cols, as.numeric(median), color=expression)) +
    geom_point()

res_sums %>% 
    group_by(cols) %>%
    mutate(ratio = as.numeric(median/min(median))) %>%
    ggplot(aes(cols, ratio, color=as.character(expression))) +
    geom_point() +
    scale_x_log10() +
    scale_color_brewer(palette="Set1") +
    ggtitle("Ratio to best method") +
    ylim(1,2)

