
context("Fragment IO v2")

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

    packed_fragments <- write_fragments_memory(raw_fragments)
    raw_fragments2 <- write_raw_fragments(packed_fragments)
    
    expect_equal(raw_fragments, raw_fragments2)
})

test_that("Packed Fragments example data round-trip", {
    in_path <- "../data/mini_fragments.tsv.gz"
    raw_fragments <- write_fragments_memory(open_10x_fragments(in_path), compressed=FALSE)
    packed_fragments <- write_fragments_memory(raw_fragments)
    raw_fragments2 <- write_fragments_memory(packed_fragments, compressed=FALSE)
    
    expect_equal(raw_fragments, raw_fragments2)
})

test_that("Unpacked Fragments example data round-trip", {
    in_path <- "../data/mini_fragments.tsv.gz"
    raw_fragments <- write_fragments_memory(open_10x_fragments(in_path), compressed=FALSE)
    packed_fragments <- write_fragments_memory(raw_fragments, compressed=FALSE)
    raw_fragments2 <- write_fragments_memory(packed_fragments, compressed=FALSE)
    
    expect_equal(raw_fragments, raw_fragments2)
})

test_that("Binary File Fragments example data round-trip", {
    dir <- withr::local_tempdir()

    in_path <- "../data/mini_fragments.tsv.gz"
    raw_fragments <- write_fragments_memory(open_10x_fragments(in_path), compressed=FALSE)
    
    write_fragments_dir(raw_fragments, file.path(dir, "unpacked"), compressed=FALSE)
    write_fragments_dir(raw_fragments, file.path(dir, "packed"), compressed=TRUE)
    
    expect_error(open_fragments_dir(file.path(dir, "packed"), compressed=FALSE))
    expect_error(open_fragments_dir(file.path(dir, "unpacked"), compressed=TRUE))
    
    unpacked <- open_fragments_dir(file.path(dir, "unpacked"), compressed=FALSE)
    packed <- open_fragments_dir(file.path(dir, "packed"), compressed=TRUE)
    
    expect_equal(raw_fragments, write_fragments_memory(unpacked, compressed=FALSE))
    expect_equal(raw_fragments, write_fragments_memory(packed, compressed=FALSE))
})


test_that("HDF5 File Fragments example data round-trip", {
    dir <- withr::local_tempdir()

    in_path <- "../data/mini_fragments.tsv.gz"
    raw_fragments <- write_fragments_memory(open_10x_fragments(in_path), compressed=FALSE)
    
    write_fragments_h5(raw_fragments, file.path(dir, "file.h5"), "unpacked", compressed=FALSE)
    write_fragments_h5(raw_fragments, file.path(dir, "file.h5"), "packed", compressed=TRUE)
    
    expect_error(open_fragments_h5(file.path(dir, "file.h5"), "packed", compressed=FALSE))
    expect_error(open_fragments_h5(file.path(dir, "file.h5"), "unpacked", compressed=TRUE))
    
    unpacked <- open_fragments_h5(file.path(dir, "file.h5"), "unpacked", compressed=FALSE)
    packed <- open_fragments_h5(file.path(dir, "file.h5"), "packed", compressed=TRUE)
    
    expect_equal(raw_fragments, write_fragments_memory(unpacked, compressed=FALSE))
    expect_equal(raw_fragments, write_fragments_memory(packed, compressed=FALSE))
})