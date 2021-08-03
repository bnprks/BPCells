
context("Fragment IO")

test_that("10x Fragments round-trip", {
    dir <- withr::local_tempdir()
    
    in_path <- "../data/mini_fragments.tsv.gz"
    out_path1 <- file.path(dir, "fragments_copy1.tsv")
    out_path2 <- file.path(dir, "fragments_copy2.tsv.gz")


    input <- open_10x_fragments(in_path)
    write_10x_fragments(input, out_path1)
    write_10x_fragments(input, out_path2)

    res <- system(sprintf("bash -c \"diff -q <(gunzip -c %s) %s\"", in_path, out_path1),
        ignore.stdout = TRUE)
    expect_equal(res, 0)

    res <- system(sprintf("bash -c \"diff -q <(gunzip -c %s) <(gunzip -c %s)\"", in_path, out_path2),
        ignore.stdout = TRUE) 
    expect_equal(res, 0)
})

test_that("10x Fragments detect tsv sorting errors", {
    library(tidyverse)
    dir <- withr::local_tempdir()
    
    in_path <- "../data/mini_fragments.tsv.gz"
    data <- read_tsv(
        in_path, 
        col_names=c("chr", "start", "end", "cell"),
        col_types="ciic",
        progress=FALSE,
        n_max=30000
    )

    error_path1 <- file.path(dir, "fragments_copy1.tsv.gz")
    data %>%
        arrange(start) %>%
        write_tsv(error_path1, col_names=FALSE)

    error_path2 <- file.path(dir, "fragments_copy2.tsv.gz")
    data %>% 
        arrange(chr, cell) %>%
        write_tsv(error_path2, col_names=FALSE)

    expect_error(open_10x_fragments(error_path1) %>% scan_fragments, "TSV not in sorted order")
    expect_error(open_10x_fragments(error_path2) %>% scan_fragments, "TSV not in sorted order")
})


test_that("Packed Fragments short chromosomes round-trip", {
    withr::local_seed(1258123)

    # Aspects for good fragments test file
    # - Include fragments starting at base 0
    # - Include chunks that take both 0 bits and 32 bits to represent
    # - Make chromosomes that have exact multiple of 128 fragments
    # - Even be abusive by including end coordinates less than the start coordinates
    chromosomes <- 129
    cells <- 40
    fragments <- list()
    for (i in seq_len(chromosomes)) {
        max_bits <- min(i, 30)
        fragments[[i]] <- list(
            start = pmin(cumsum(sample.int(2^max_bits, i, replace=TRUE) - 1), .Machine$integer.max),
            end = sample.int(2^max_bits, i, replace=TRUE),
            cell_id = sample.int(cells, i, replace=TRUE) - 1
        )
    }

    raw_fragments <- new("RawFragments", 
        fragments = fragments,
        chr_names = sprintf("chr%d", seq_len(chromosomes)),
        cell_names = sprintf("cell%d", seq_len(cells))
    )

    packed_fragments <- write_packed_fragments(raw_fragments)
    raw_fragments2 <- write_raw_fragments(packed_fragments)
    
    expect_equal(raw_fragments, raw_fragments2)
})


test_that("Packed Fragments example data round-trip", {
    in_path <- "../data/mini_fragments.tsv.gz"
    raw_fragments <- write_raw_fragments(open_10x_fragments(in_path))
    packed_fragments <- write_packed_fragments(raw_fragments)
    raw_fragments2 <- write_raw_fragments(packed_fragments)
    
    expect_equal(raw_fragments, raw_fragments2)
})

test_that("10x Fragments files can skip chromosomes", {
    in_fragments <- open_10x_fragments("../data/mini_fragments.tsv.gz")
    
    raw_fragments <- write_raw_fragments(in_fragments)

    chromosome_selection <- sample.int(length(raw_fragments@chr_names), length(raw_fragments@chr_names)/2)
    chr_names <- raw_fragments@chr_names[chromosome_selection]

    raw_fragments2 <- write_raw_fragments(select_chromosomes(in_fragments, chromosome_selection))
    raw_fragments@fragments <- raw_fragments@fragments[chromosome_selection]
    raw_fragments@chr_names <- chr_names

    expect_equal(raw_fragments2, raw_fragments)
})

test_that("10x Fragments works with comments", {
    dir <- withr::local_tempdir()
    
    tmp <- file.path(dir, "fragments_with_comment.tsv")
    cat("# id=experiment_id\n# description=lots of stuff\n", file=tmp)
    
    in_path <- "../data/mini_fragments.tsv.gz"
    system(sprintf("gunzip -c %s >> %s", in_path, tmp))

    raw_1 <- open_10x_fragments(in_path) %>% write_raw_fragments()
    raw_2 <- open_10x_fragments(tmp) %>% write_raw_fragments()

    expect_equal(raw_1, raw_2)
})

test_that("10x Fragment 5th column output works", {
    dir <- withr::local_tempdir()
    
    in_path <- "../data/mini_fragments.tsv.gz"
    out_path1 <- file.path(dir, "fragments_copy1.tsv")
    out_path2 <- file.path(dir, "fragments_copy2.tsv")

    input <- open_10x_fragments(in_path)
    write_10x_fragments(input, out_path1, append_5th_column=TRUE)
    write_10x_fragments(input, out_path2, append_5th_column=FALSE)

    read_first_line <- function(path, gzip=FALSE) {
        if (gzip) f <- gzfile(path, "r")
        else f <- file(path, "r")
        res <- readLines(f, n=1)
        close(f)
        return(res)
    }
    res_line <- read_first_line(in_path, gzip=TRUE)
    line1 <- read_first_line(out_path1)
    line2 <- read_first_line(out_path2)
    expect_equal(paste0(res_line, "\t0"), line1)
    expect_equal(res_line, line2)
})