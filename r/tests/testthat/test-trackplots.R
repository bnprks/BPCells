

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