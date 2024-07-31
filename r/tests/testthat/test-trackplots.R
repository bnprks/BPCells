

get_ggplot_labels <- function(plot) {
    p <- ggplot2::ggplot_build(plot)
    if (length(p$layout$facet$params$facets) == 0) {
        labels_df <- tibble::tibble("(all)"="(all)")
    } else {
        labels_df <- p$layout$layout[names(p$layout$facet$params$facets)]
    }
    attr(labels_df, "facet") <- "wrap"
    attr(labels_df, "type") <- "row"
    as.character(do.call(cbind, p$layout$facet$params$labeller(labels_df)))
}

test_that("gene_region works", {
    genes <- tibble::tibble(
        chr=c("chr1", "chr2"), 
        start=c(1000, 2000),
        end=c(1500, 2500),
        strand=c("+", "-"),
        gene_name=c("gene1", "gene2")
    )
    # Test symmetrical and asymmetric flanking lengths with +/- strand
    expect_identical(gene_region(genes, "gene1", extend_bp=10), 
                     list(chr="chr1", start=990, end=1510))
    expect_identical(gene_region(genes, "gene1", extend_bp=c(10, 100)), 
                     list(chr="chr1", start=990, end=1600))
    expect_identical(gene_region(genes, "gene2", extend_bp=10), 
                     list(chr="chr2", start=1990, end=2510))
    expect_identical(gene_region(genes, "gene2", extend_bp=c(10, 100)), 
                     list(chr="chr2", start=1900, end=2510))
})

test_that("set_trackplot_label works", {
    # Check setting facet label on a trackplot input
    p <- trackplot_scalebar("chr1:1-1000")
    p2 <- set_trackplot_label(p, "Scale bar")
    expect_identical(get_ggplot_labels(p2), c("Scale bar"))

    expect_error(set_trackplot_label(p, c("label1", "label2")))

    # Check setting facet label on un-facetted plot
    p <- ggplot2::ggplot(datasets::mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
    p2 <- set_trackplot_label(p, "My Label")
    expect_identical(get_ggplot_labels(p2), c("My Label"))

    # Check changing facet labels on faceted plot
    p1 <- p + ggplot2::facet_wrap(~cyl)
    expect_identical(get_ggplot_labels(p1), c("4", "6", "8"))
    p2 <- set_trackplot_label(p1, c("a", "b", "c"))
    expect_identical(get_ggplot_labels(p2), c("a", "b", "c"))

    # Check for error on incorrect number of labels provided
    p3 <- set_trackplot_label(p1, c("a", "b"))
    expect_error(get_ggplot_labels(p3), "set_trackplot_label.*[0-9]+ labels were provided")
})

test_that("get/set trackplot height works", {
    # Check setting height on non-trackplot
    p <- ggplot2::ggplot(datasets::mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
    expect_identical(get_trackplot_height(p), ggplot2::unit(1L, "null"))
    p <- set_trackplot_height(p, ggplot2::unit(2L, "pt"))
    expect_identical(get_trackplot_height(p), ggplot2::unit(2L, "pt"))
    p <- set_trackplot_height(p, ggplot2::unit(3L, "pt"))
    expect_identical(get_trackplot_height(p), ggplot2::unit(3L, "pt"))
})

test_that("creating continuous trackplot arrows works", {
    data1 <- tibble::tibble(
        start = c(10000, 20000, 10000),
        end = c(15000, 25000, 15000),
        strand = c(TRUE, FALSE, FALSE),
    )
    data2 <- data1 %>% tibble::add_column(
        chr = c("chr1", "chr1", "chr1"),
        label = c("A", "B", "C"))
    data3 <- tibble::tibble(
        start = c(10000, 10050, 20001),
        end = c(10010, 10060, 20005),
        strand = c(TRUE, FALSE, FALSE),
    )
    region <- tibble::tibble(
        chr = "chr1",
        start = 10000,
        end = 20000
    )
    segs <- trackplot_create_arrow_segs(data1, region, 100)
    segs_expected <- tibble::tibble(
        start = c(seq(10000, 14900, 100), seq(10000, 14900, 100)),
        end = c(seq(10100, 15000, 100), seq(10100, 15000, 100)),
        strand = c(rep(TRUE, 50), rep(FALSE, 50))
    )
    expect_identical(segs, segs_expected)

    segs_expected_head_only <- tibble::tibble(
        start = c(seq(10100, 15000, 100) - 1e-4, seq(10000, 14900, 100)),
        end = c(seq(10100, 15000, 100), seq(10000, 14900, 100) + 1e-4),
        strand = c(rep(TRUE, 50), rep(FALSE, 50))
    )
    expect_identical(
        trackplot_create_arrow_segs(data1, region, 100, head_only=TRUE), 
        segs_expected_head_only
    )

    segs_with_metadata <- trackplot_create_arrow_segs(data2, region, 100)
    segs_expected_with_metadata <- segs_expected %>% 
        tibble::add_column(
            chr = c(rep("chr1", 100)),
            label = c(rep("A", 50), rep("C", 50)))
    expect_identical(segs_with_metadata, segs_expected_with_metadata)

    segs_small <- trackplot_create_arrow_segs(data3, region, 100)
    segs_expected_small <- tibble::tibble(
        start = c(10000, 10050),
        end = c(10010, 10060),
        strand = c(TRUE, FALSE)
    )
    expect_identical(segs_small, segs_expected_small)

    segs_expected_small_head_only <- tibble::tibble(
        start = c(10010-1e-4, 10050),
        end = c(10010, 10050+1e-4),
        strand = c(TRUE, FALSE)
    )
    expect_identical(
        trackplot_create_arrow_segs(data3, region, 100, head_only=TRUE),
        segs_expected_small_head_only
    )
})

test_that("trackplot_genome_annotation doesn't crash", {
    test_data <- tibble::tibble(
        chr = rep.int("chr1", 4),
        start = c(1, 4, 5, 7),
        end = start + 2,
        strand = c("+", "-", "+", "-"),
        cont_data = 1:4,
        discrete_data = c("a", "b", "b", "a")
    )
    expect_no_condition(ggplot2::ggplot_build(trackplot_genome_annotation(test_data, "chr1:0-10")))
    
    # Check color variations
    expect_no_condition(ggplot2::ggplot_build(trackplot_genome_annotation(test_data, "chr1:0-10", color_by="discrete_data")))
    expect_no_condition(ggplot2::ggplot_build(trackplot_genome_annotation(test_data, "chr1:0-10", color_by="cont_data")))
    expect_no_condition(ggplot2::ggplot_build(trackplot_genome_annotation(test_data, "chr1:0-10", color_by="strand")))
    expect_no_condition(ggplot2::ggplot_build(trackplot_genome_annotation(test_data, "chr1:0-10", color_by=test_data$discrete_data)))
    expect_no_condition(ggplot2::ggplot_build(trackplot_genome_annotation(test_data, "chr1:0-10", color_by=test_data$cont_data)))
    
    expect_no_condition(ggplot2::ggplot_build(trackplot_genome_annotation(test_data, "chr1:0-10", color_by="cont_data", colors=c("#222", "#ccc"))))
    expect_no_condition(ggplot2::ggplot_build(trackplot_genome_annotation(test_data, "chr1:0-10", color_by="discrete_data", colors=c("#222", "#ccc"))))
    expect_error(trackplot_genome_annotation(test_data, "chr1:0-10", color_by="discrete_data", colors=c("#222")), "Insufficient values in manual color scale")

    # Check label variations
    expect_no_condition(ggplot2::ggplot_build(trackplot_genome_annotation(test_data, "chr1:0-10", label_by="discrete_data")))
    expect_no_condition(ggplot2::ggplot_build(trackplot_genome_annotation(test_data, "chr1:0-10", label_by="cont_data")))
    expect_no_condition(ggplot2::ggplot_build(trackplot_genome_annotation(test_data, "chr1:0-10", label_by="strand")))
    expect_no_condition(ggplot2::ggplot_build(trackplot_genome_annotation(test_data, "chr1:0-10", label_by=test_data$discrete_data)))

    # Check data re-use for multiple metadata types
    expect_no_condition(ggplot2::ggplot_build(trackplot_genome_annotation(test_data, "chr1:0-10", label_by="strand", color_by="strand", show_strand=TRUE)))

    # Check empty plot with/without arrows
    expect_no_condition(ggplot2::ggplot_build(trackplot_genome_annotation(test_data, "chr2:0-10", label_by="strand", color_by="strand", show_strand=TRUE)))
    expect_no_condition(ggplot2::ggplot_build(trackplot_genome_annotation(test_data, "chr2:0-10", label_by="strand", color_by="strand", show_strand=FALSE)))

    # Check return data works with/without show_strand
    for (region in c("chr1:0-10", "chr2:0-10")) {
        for (show_strand in c(TRUE, FALSE)) {
            expect_named(
                trackplot_genome_annotation(test_data, region, label_by="strand", color_by="strand", show_strand=show_strand, return_data=TRUE), 
                c("data", "arrows")
            )
        }
    }

})

test_that("trackplot_loop doesn't crash", {
    test_data <- tibble::tibble(
        chr = rep.int("chr1", 4),
        start = c(1, 4, 5, 7),
        end = start + 2,
        strand = c("+", "-", "+", "-"),
        cont_data = 1:4,
        discrete_data = c("a", "b", "b", "a")
    )

    expect_no_condition(ggplot2::ggplot_build(trackplot_loop(test_data, "chr1:0-10", color_by="discrete_data")))
    expect_no_condition(ggplot2::ggplot_build(trackplot_loop(test_data, "chr1:0-10", color_by="cont_data")))
    expect_no_condition(ggplot2::ggplot_build(trackplot_loop(test_data, "chr1:0-10", color_by="strand")))
    expect_no_condition(ggplot2::ggplot_build(trackplot_loop(test_data, "chr1:0-10", color_by=test_data$discrete_data)))
    expect_no_condition(ggplot2::ggplot_build(trackplot_loop(test_data, "chr1:0-10", color_by=test_data$cont_data)))
    
    expect_no_condition(ggplot2::ggplot_build(trackplot_loop(test_data, "chr1:0-10", color_by="cont_data", colors=c("#222", "#ccc"))))
    expect_no_condition(ggplot2::ggplot_build(trackplot_loop(test_data, "chr1:0-10", color_by="discrete_data", colors=c("#222", "#ccc"))))
    expect_error(trackplot_loop(test_data, "chr1:0-10", color_by="discrete_data", colors=c("#222")), "Insufficient values in manual color scale")

    # Check return data works
    for (region in c("chr1:0-10", "chr2:0-10")) {
        expect_is(
            trackplot_loop(test_data, region, color_by="strand", return_data=TRUE), 
            "tbl_df"
        )
    }
    
})

test_that("trackplot_coverage doesn't crash", {
    frags <- open_fragments_10x("../data/mini_fragments.tsv.gz") %>%
        write_fragments_memory()
    cell_names <- cellNames(frags)
    frags_metadata <- tibble::tibble(
      id = cell_names,
      atac_reads = as.integer(seq(from=0, to=10, length.out=length(cell_names))),
      cell_type = rep(c("NK", "B", "T"), ceiling(length(cell_names)/3))[1:length(cell_names)]
    )
    # general region
    region1 <- tibble::tibble(
        chr = "chr15",
        start = 30000000,
        end = 52806155,
    )
    # for small motifs
    region2 <- tibble::tibble(
        chr = "chr15",
        start = 31283779,
        end = 31283870
    )
    frags_metadata_one_cluster <- frags_metadata[frags_metadata$cell_type == "T",]
    # checks for only one cluster
    frags_one_cluster <- select_cells(frags,  frags_metadata_one_cluster$id)
    for (region in list(region1, region2)) {
        expect_no_condition(
            ggplot2::ggplot_build(trackplot_coverage(
                  frags,
                  region=region,
                  groups=frags_metadata$cell_type,
                  cell_read_counts=frags_metadata$atac_reads,
                  bins=500
                )
            )
        )
        
        expect_no_condition(
            ggplot2::ggplot_build(trackplot_coverage(
                  frags_one_cluster,
                  region=region,
                  groups=frags_metadata_one_cluster$cell_type,
                  cell_read_counts=frags_metadata_one_cluster$atac_reads,
                  bins=500
                )
            )
        )
        # Check return data works
        expect_is(
            trackplot_coverage(
                frags_one_cluster,
                region=region,
                groups=frags_metadata_one_cluster$cell_type,
                cell_read_counts=frags_metadata_one_cluster$atac_reads,
                bins=500,
                return_data=TRUE
            ),
            "tbl_df"
        )
    }
})