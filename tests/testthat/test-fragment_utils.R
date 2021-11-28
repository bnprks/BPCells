
raw_fragments_to_tibble <- function(raw_fragments) {
    dplyr::bind_rows(raw_fragments$fragments, .id="chr")
}

tibble_to_raw_fragments <- function(x, chr_names, cell_names) {
    fragments <- lapply(seq_along(chr_names), function(i) {
        dplyr::filter(x, chr == i) %>%
            {list(start=.$start, end=.$end, cell_id=.$cell_id)}
    })
    new("RawFragments",
        fragments = fragments,
        chr_names = chr_names,
        cell_names = cell_names
    )
}

test_that("Chromosome name/index select works", {
    # Test lots of short chromosomes with index and name selection
    withr::local_seed(1258123)

    nfrags <- 10e3
    short_chrs_tibble <- tibble::tibble(
        chr = sample.int(500, nfrags, replace=TRUE),
        start = sample.int(10000, nfrags, replace=TRUE),
        end = start + sample.int(500, nfrags, replace=TRUE),
        cell_id = sample.int(500, nfrags, replace=TRUE)
    ) %>% dplyr::arrange(chr, start)

    short_chrs <- tibble_to_raw_fragments(
        short_chrs_tibble,
        paste0("chr", 1:500),
        paste0("cell", 1:500)
    )

    chr_selection <- sample.int(500, 250, replace=FALSE)

    short_chrs_ans <- short_chrs
    short_chrs_ans@fragments <- short_chrs@fragments[chr_selection]
    short_chrs_ans@chr_names <- short_chrs@chr_names[chr_selection]

    short_chrs_2 <- select_chromosomes(short_chrs, chr_selection) %>%
        write_raw_fragments()

    short_chrs_3 <- select_chromosomes(short_chrs, paste0("chr", chr_selection)) %>%
        write_raw_fragments()

    expect_equal(short_chrs_ans, short_chrs_2)
    expect_equal(short_chrs_ans, short_chrs_3)
})

test_that("Chromosome name/index select works", {
    # Test lots of short chromosomes with index and name selection
    withr::local_seed(1258123)

    nfrags <- 10e3
    short_chrs_tibble <- tibble::tibble(
        chr = sample.int(500, nfrags, replace=TRUE),
        start = sample.int(10000, nfrags, replace=TRUE),
        end = start + sample.int(500, nfrags, replace=TRUE),
        cell_id = sample.int(500, nfrags, replace=TRUE)
    ) %>% dplyr::arrange(chr, start)

    short_chrs <- tibble_to_raw_fragments(
        short_chrs_tibble,
        paste0("chr", 1:500),
        paste0("cell", 1:500)
    )

    chr_selection <- sample.int(500, 250, replace=FALSE)

    short_chrs_ans <- short_chrs
    short_chrs_ans@fragments <- short_chrs@fragments[chr_selection]
    short_chrs_ans@chr_names <- short_chrs@chr_names[chr_selection]

    short_chrs_2 <- select_chromosomes(short_chrs, chr_selection) %>%
        write_raw_fragments()

    short_chrs_3 <- select_chromosomes(short_chrs, paste0("chr", chr_selection)) %>%
        write_raw_fragments()

    expect_equal(short_chrs_ans, short_chrs_2)
    expect_equal(short_chrs_ans, short_chrs_3)
})


test_that("Cell name/index select works", {
    # Test cells with index and name selection, and try explicitly making
    # chr 2 only contain cells that get filtered out, and
    # chr 3 only contain cells that get kept
    withr::local_seed(1258123)

    nfrags <- 1e4
    nchrs <- 20
    ncells <- 100
    cell_selection <- sample.int(ncells, ncells/2, replace=FALSE)-1

    frags_tibble <- tibble::tibble(
        chr = sample.int(nchrs, nfrags, replace=TRUE),
        start = sample.int(10000, nfrags, replace=TRUE),
        end = start + sample.int(500, nfrags, replace=TRUE),
        cell_id = sample.int(ncells, nfrags, replace=TRUE)-1
    ) %>% dplyr::arrange(chr, start) %>%
        dplyr::filter(
            (chr == 2 & !(cell_id %in% cell_selection)) |
            (chr == 3 & cell_id %in% cell_selection) |
            !(chr %in% c(2, 3))
        )

    frags <- tibble_to_raw_fragments(
        frags_tibble,
        paste0("chr", seq_len(nchrs)),
        paste0("cell", seq_len(ncells))
    )

    frags_ans <- tibble_to_raw_fragments(
        frags_tibble %>% 
            dplyr::filter(cell_id %in% cell_selection) %>%
            dplyr::mutate(cell_id = match(cell_id, cell_selection) - 1),
        frags@chr_names,
        frags@cell_names[cell_selection+1]
    )

    frags_2 <- select_cells(frags, cell_selection+1) %>%
        write_raw_fragments()
    frags_3 <- select_cells(frags, paste0("cell", cell_selection+1)) %>%
        write_raw_fragments()

    expect_equal(frags_ans, frags_2)
    expect_equal(frags_ans, frags_3)
})

test_that("GRanges conversion round-trips", {
    frags <- open_10x_fragments("../data/mini_fragments.tsv.gz")
    raw <- write_raw_fragments(frags)
    ranges <- as(raw, "GRanges")
    raw2 <- RawFragments(ranges)
    expect_equal(raw, raw2)
})

test_that("Insertions iterator works", {
    frags <- open_10x_fragments("../data/mini_fragments.tsv.gz")
    raw <- write_raw_fragments(frags)

    scan1 <- scan_fragments_modulo_cpp(iterate_fragments(raw))
    scan2 <- scan_insertions2_cpp(iterate_fragments(raw))

    expect_equal(scan1[4], scan2[3])
    expect_equal(scan1[1]*2, scan2[1])
})
