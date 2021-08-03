
test_that("Basic insertion matrix succeeds", { 
    # chr1 tests accuracy of start coordinate overlaps alone
    # chr2 tests accuracy of start+end coordinate overlaps combined 
    # The ordering tests string matching of chromosome name & correct
    # sorting of output rows

    chr1_list <- list()
    for (i in 1:5) {
        chr1_list[[i]] <- tibble::tibble(
            cell_id = i - 1,
            starts = rep(i:5, i),
            ends = 1001 + i
        )
    }
    chr1_coords <- dplyr::bind_rows(chr1_list) %>%
       dplyr::arrange(starts) %>% 
       {list(start=.$starts, end=.$ends, cell_id=.$cell_id)}


    chr2_coords <- tibble::tribble(
        ~cell_id, ~starts, ~ends,
        0, 9, 21,
        1, 9, 20,
        2, 10, 21,
        3, 10, 20
    ) %>% dplyr::arrange(starts) %>%
        {list(start=.$starts, end=.$ends, cell_id=.$cell_id)}

    raw_fragments <- new("RawFragments",
        fragments = list(chr2_coords, chr1_coords),
        chr_names = c("chr2", "chr1"),
        cell_names = sprintf("cell%d", 1:5)
    )

    res <- overlapMatrix(
        raw_fragments,
        list(
            chr=c("chr1", "chr2", "chr1", "chr1"), 
            start=c(1004, 10, 1002, 3), 
            end=c(1006, 20, 1005, 5)),
        convert_to_0_based_coords=FALSE
    ) %>% as("dgCMatrix")

    expect_s4_class(res, "dgCMatrix")
    
    answer <- matrix(c(
        0,0,0,8,5,
        0,1,1,2,0,
        0,8,9,8,0,
        2,4,6,4,0
    ), ncol=4)

    my_answer <- as.matrix(res)
    attr(my_answer, "dimnames") <- NULL

    expect_equal(my_answer, answer)
})


test_that("Out of range peaks work", { 
    # chr1 tests having some peaks before the fragments start
    # chr2 tests having some peaks after the end of the fragments
    # chr3 tests having some fragments but no peaks
    chr1_list <- list()
    for (i in 1:5) {
        chr1_list[[i]] <- tibble::tibble(
            cell_id = i - 1,
            starts = i:5 + 400,
            ends = 1001 + i
        )
    }
    chr1_coords <- dplyr::bind_rows(chr1_list) %>%
       dplyr::arrange(starts) %>% 
       {list(start=.$starts, end=.$ends, cell_id=.$cell_id)}

    chr2_coords <- tibble::tribble(
        ~cell_id, ~starts, ~ends,
        0, 9, 21,
        1, 9, 20,
        2, 10, 21,
        3, 10, 20
    ) %>% dplyr::arrange(starts) %>%
        {list(start=.$starts, end=.$ends, cell_id=.$cell_id)}

    raw_fragments <- new("RawFragments",
        fragments = list(chr2_coords, chr1_coords, chr2_coords),
        chr_names = c("chr2", "chr1", "chr3"),
        cell_names = sprintf("cell%d", 1:5)
    )

    res <- overlapMatrix(
        raw_fragments,
        list(
            chr=c("chr1", "chr1", "chr2", "chr2"), 
            start=c(1004, 3, 1002, 10), 
            end=c(1006, 5, 1005, 20)),
        convert_to_0_based_coords = FALSE
    ) %>% as("dgCMatrix")

    expect_s4_class(res, "dgCMatrix")
    
    answer <- matrix(c(
        0,0,0,2,1,
        0,0,0,0,0,
        0,0,0,0,0,
        0,1,1,2,0
    ), ncol=4)

    my_answer <- as.matrix(res)
    attr(my_answer, "dimnames") <- NULL

    expect_equal(my_answer, answer)
})
