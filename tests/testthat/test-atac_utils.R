
test_that("range_overlaps works", {
    a <- tibble::tibble(
        chr=c("chr1", "chr2", "chr1"),
        start=c(1, 5, 10),
        end=c(20, 20, 20)
    )
    b <- tibble::tibble(
        chr=c("chr2", "chr2", "chr1", "chr1", "chr1", "chr1"),
        start=c(1, 25, 2, 11, 0, 2),
        end=c(10, 30, 20, 25, 25, 4)
    )

    expected <- tibble::tibble(
        from = c(1,1,1,1,2,3,3,3),
        to = c(3,4,5,6,1,3,4,5)
    )

    expect_equal(range_overlaps(a,b), expected)
})