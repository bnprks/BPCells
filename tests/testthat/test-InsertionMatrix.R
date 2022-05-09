library(dplyr)

tibble_to_fragments <- function(x, chr_names, cell_names) {
    x %>% mutate(
        chr = factor(chr_names[chr], levels=chr_names),
        cell_id = factor(cell_names[cell_id], levels=cell_names)
    ) %>% convert_to_fragments()
}

test_that("Basic insertion matrix succeeds", { 
    # chr1 tests accuracy of start coordinate overlaps alone
    # chr2 tests accuracy of start+end coordinate overlaps combined 
    # The ordering tests string matching of chromosome name & correct
    # sorting of output rows

    chr1_list <- list()
    for (i in 1:5) {
        chr1_list[[i]] <- tibble::tibble(
            cell_id = i - 1,
            start = rep(i:5, i),
            end = 1001 + i
        )
    }
    chr1_coords <- dplyr::bind_rows(chr1_list) %>%
       dplyr::arrange(start) %>% 
       mutate(chr = 2, cell_id = cell_id+1)

    chr2_coords <- tibble::tribble(
        ~cell_id, ~start, ~end,
        0, 9, 21,
        1, 9, 20,
        2, 10, 21,
        3, 10, 20
    ) %>% dplyr::arrange(start) %>%
        mutate(chr=1, cell_id=cell_id+1)

    raw_fragments <- tibble_to_fragments(
        bind_rows(chr2_coords, chr1_coords),
        chr_names = c("chr2", "chr1"),
        cell_names = sprintf("cell%d", 1:5)
    )

    res <- peakMatrix(
        raw_fragments,
        list(
            chr = c("chr2", "chr1", "chr1", "chr1"), 
            start = c(10, 3, 1002, 1004),
            end = c(20, 5, 1005, 1006))
    ) %>% as("dgCMatrix")

    expect_s4_class(res, "dgCMatrix")
    
    answer <- matrix(c(
        0,1,1,2,0,
        2,4,6,4,0,
        0,8,9,8,0,
        0,0,0,8,5
    ), ncol=4)

    my_answer <- as.matrix(res)
    attr(my_answer, "dimnames") <- NULL

    expect_equal(my_answer, answer)
})


test_that("Out of range peaks work", { 
    # chr1 tests having some peaks before the fragments start
    # chr2 tests having some peaks after the end of the fragments
    # chr3 tests having some fragments but no peaks
    # chr4 tests having some fragments, with a duplicated 1bp-wide peak
    chr1_list <- list()
    for (i in 1:5) {
        chr1_list[[i]] <- tibble::tibble(
            cell_id = i - 1,
            start = i:5 + 400,
            end = 1001 + i
        )
    }
    chr1_coords <- dplyr::bind_rows(chr1_list) %>%
       dplyr::arrange(start) %>% 
       mutate(chr = 2, cell_id = cell_id+1)

    chr2_coords <- tibble::tribble(
        ~cell_id, ~start, ~end,
        0, 9, 21,
        1, 9, 20,
        2, 10, 21,
        3, 10, 20
    ) %>% dplyr::arrange(start) %>%
        mutate(chr=1, cell_id=cell_id+1)

    chr3_coords <- mutate(chr2_coords, chr=3)
    chr4_coords <- mutate(chr2_coords, chr=4)

    raw_fragments <- tibble_to_fragments(
        bind_rows(chr2_coords, chr1_coords, chr3_coords, chr4_coords),
        chr_names = c("chr2", "chr1", "chr3", "chr4"),
        cell_names = sprintf("cell%d", 1:5)
    )

    res <- peakMatrix(
        raw_fragments,
        list(
            chr=c("chr2", "chr2", "chr1", "chr1", "chr4", "chr4"), 
            start=c(10, 1002, 3, 1004, 1, 1), 
            end=c(20, 1005, 5, 1006, 2, 2))
    ) %>% as("dgCMatrix")

    expect_s4_class(res, "dgCMatrix")
    
    answer <- matrix(c(
        0,1,1,2,0,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,2,1,
        0,0,0,0,0,
        0,0,0,0,0
    ), ncol=6)

    my_answer <- as.matrix(res)
    attr(my_answer, "dimnames") <- NULL

    expect_equal(my_answer, answer)
})



test_that("Basic tile matrix works", { 
    # chr1 tests having some peaks before the fragments start
    # chr2 tests having some peaks after the end of the fragments
    # chr3 tests having some fragments but no peaks
    # chr4 tests having some fragments, with a duplicated 1bp-wide peak
    chr1_list <- list()
    for (i in 1:5) {
        chr1_list[[i]] <- tibble::tibble(
            cell_id = i - 1,
            start = i:5 + 400,
            end = 1001 + i
        )
    }
    chr1_coords <- dplyr::bind_rows(chr1_list) %>%
       dplyr::arrange(start) %>% 
       mutate(chr = 2, cell_id = cell_id+1)

    chr2_coords <- tibble::tribble(
        ~cell_id, ~start, ~end,
        0, 9, 21,
        1, 9, 20,
        2, 10, 21,
        3, 10, 20
    ) %>% dplyr::arrange(start) %>%
        mutate(chr=1, cell_id=cell_id+1)

    chr3_coords <- mutate(chr2_coords, chr=3)
    chr4_coords <- mutate(chr2_coords, chr=4)

    raw_fragments <- tibble_to_fragments(
        bind_rows(chr2_coords, chr1_coords, chr3_coords, chr4_coords),
        chr_names = c("chr2", "chr1", "chr3", "chr4"),
        cell_names = sprintf("cell%d", 1:5)
    )

    res <- tileMatrix(
        raw_fragments,
        list(
            chr=c("chr2", "chr1", "chr1", "chr1"), 
            start=c(10, 3, 402, 1004), 
            end=c(20, 5, 405, 1006),
            tile_width = c(4,1,2,2))
    ) %>% as("dgCMatrix")

    expect_s4_class(res, "dgCMatrix")
    
    answer <- matrix(c(
        0,0,1,1,0, # 10-14
        0,0,0,0,0, # 14-18
        0,1,0,1,0, # 18-19
        0,0,0,0,0, # 3
        0,0,0,0,0, # 4
        2,2,1,0,0, # 402-403
        1,1,1,1,0, # 404
        0,0,0,2,1 # 1004-1006
    ), ncol=8)

    my_answer <- as.matrix(res)
    attr(my_answer, "dimnames") <- NULL

    expect_equal(my_answer, answer)
})