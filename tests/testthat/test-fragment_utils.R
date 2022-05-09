raw_fragments_to_tibble <- function(raw_fragments) {
    dplyr::bind_rows(raw_fragments$fragments, .id="chr")
}

tibble_to_fragments <- function(x, chr_names, cell_names) {
    x %>% dplyr::mutate(
        chr = factor(chr_names[chr], levels=chr_names),
        cell_id = factor(cell_names[cell_id], levels=cell_names)
    ) %>% convert_to_fragments()
}

test_that("Chromosome name/index select works", {
    # Test lots of short chromosomes with index and name selection
    library(GenomicRanges)

    withr::local_seed(1258123)
    
    nfrags <- 10e3
    short_chrs_tibble <- tibble::tibble(
        chr = sample.int(500, nfrags, replace=TRUE),
        start = sample.int(10000, nfrags, replace=TRUE),
        end = start + sample.int(500, nfrags, replace=TRUE),
        cell_id = sample.int(500, nfrags, replace=TRUE)
    ) %>% dplyr::arrange(chr, start)

    short_chrs <- tibble_to_fragments(
        short_chrs_tibble,
        paste0("chr", 1:500),
        paste0("cell", 1:500)
    )

    chr_selection <- sample.int(500, 250, replace=FALSE)

    short_chrs_ans <- as(short_chrs, "GRanges")
    short_chrs_ans <- short_chrs_ans[as.character(GenomicRanges::seqnames(short_chrs_ans)) %in% paste0("chr", chr_selection)] %>%
        convert_to_fragments() %>% as("GRanges")
    
    short_chrs_2 <- select_chromosomes(short_chrs, chr_selection) %>%
        as("GRanges")

    short_chrs_3 <- select_chromosomes(short_chrs, paste0("chr", chr_selection)) %>%
        as("GRanges")

    expect_equal(short_chrs_ans, short_chrs_2, check.attributes=FALSE)
    expect_equal(short_chrs_ans, short_chrs_3, check.attributes=FALSE)
})

test_that("Cell name/index select works", {
    # Test cells with index and name selection, and try explicitly making
    # chr 2 only contain cells that get filtered out, and
    # chr 3 only contain cells that get kept
    withr::local_seed(1258123)

    nfrags <- 1e4
    nchrs <- 20
    ncells <- 100
    cell_selection <- sample.int(ncells, ncells/2, replace=FALSE)

    frags_tibble <- tibble::tibble(
        chr = sample.int(nchrs, nfrags, replace=TRUE),
        start = sample.int(10000, nfrags, replace=TRUE),
        end = start + sample.int(500, nfrags, replace=TRUE),
        cell_id = sample.int(ncells, nfrags, replace=TRUE)
    ) %>% dplyr::arrange(chr, start) %>%
        dplyr::filter(
            (chr == 2 & !(cell_id %in% cell_selection)) |
            (chr == 3 & cell_id %in% cell_selection) |
            !(chr %in% c(2, 3))
        )

    frags <- tibble_to_fragments(
        frags_tibble,
        paste0("chr", seq_len(nchrs)),
        paste0("cell", seq_len(ncells))
    )

    frags_ans <- tibble_to_fragments(
        frags_tibble %>% 
            dplyr::filter(cell_id %in% cell_selection) %>%
            dplyr::mutate(cell_id = match(cell_id, cell_selection)),
        frags@chr_names,
        frags@cell_names[cell_selection]
    ) %>% write_fragments_memory() %>% as("GRanges")

    frags_2 <- select_cells(frags, cell_selection) %>%
        as("GRanges")
    frags_3 <- select_cells(frags, paste0("cell", cell_selection)) %>%
        as("GRanges")

    expect_equal(frags_ans, frags_2)
    expect_equal(frags_ans, frags_3)
})

test_that("GRanges conversion round-trips", {
    frags <- open_fragments_10x("../data/mini_fragments.tsv.gz")
    raw <- write_fragments_memory(frags)
    ranges <- as(raw, "GRanges")
    ranges2 <- as(frags, "GRanges")
    expect_equal(ranges, ranges2)
    raw2 <- write_fragments_memory(convert_to_fragments(ranges))
    expect_equal(raw, raw2)
})

test_that("Generic methods work", {
    frags <- open_fragments_10x("../data/mini_fragments.tsv.gz") %>%
        write_fragments_memory()

    dir <- withr::local_tempdir()

    ident_transforms <- list(
        write_memory = write_fragments_memory(frags, compress=TRUE),
        write_memory_unpacked = write_fragments_memory(frags, compress=FALSE),
        write_dir = write_fragments_dir(frags, file.path(dir, "fragments_identical")),
        write_h5 = write_fragments_hdf5(frags, file.path(dir, "fragments_identical.h5")),
        #write_bed = write_fragments_10x(frags, file.path(dir, "fragments_identical.tsv")),
        shift = shift_fragments(frags),
        chrSelectName = select_chromosomes(frags, chrNames(frags)),
        chrSelectIdx = select_chromosomes(frags, seq_along(chrNames(frags))),
        cellSelectName = select_cells(frags, cellNames(frags)),
        cellSelectIdx = select_cells(frags, seq_along(cellNames(frags))),
        chrRename = new("ChrRename", fragments=frags, chr_names=chrNames(frags)),
        cellRename = new("CellRename", fragments=frags, cell_names=cellNames(frags))
    )

    for (i in seq_along(ident_transforms)) {
        t <- ident_transforms[[i]]
        short_description(t)
        expect_equal(chrNames(t), chrNames(frags))
        expect_equal(cellNames(t), cellNames(frags))
        expect_equal(write_fragments_memory(t), frags)
        
        chrNames(t) <- paste0("new-", chrNames(frags))
        cellNames(t) <- paste0("new-", cellNames(frags))
        expect_equal(chrNames(t), paste0("new-", chrNames(frags)))
        expect_equal(cellNames(t), paste0("new-", cellNames(frags)))
        n <- write_fragments_memory(t)
        expect_equal(chrNames(n), paste0("new-", chrNames(frags)))
        expect_equal(cellNames(n), paste0("new-", cellNames(frags)))
    }

})