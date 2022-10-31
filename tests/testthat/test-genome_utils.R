
test_that("Gene name translation works", {
    # Test each function with a warning and without a warning
    canonical <- expect_warning(
        canonical_gene_symbol(c("CD45", "CD8", "CD4", "CD3")),
        "Could not match"
    )
    expect_equal(
        canonical,
        c(CD45 = "PTPRC", CD8 = "CD8A", CD4 = "CD4", CD3 = NA)
    )
    expect_equal(
        canonical_gene_symbol(c("CD45", "CD8", "CD4")),
        c(CD45 = "PTPRC", CD8 = "CD8A", CD4 = "CD4")
    )

    match_res <- expect_warning(
        match_gene_symbol(
            c("CD8", "CD4", "CD45", "CD3", "SOX2"),
            c("ENSG00000081237.19", "ENSG00000153563.15", "ENSG00000010610.9", "ENSG00000288825")
        ),
        "Could not match"
    )
    expect_equal(
        match_res,
        c(2, 3, 1, NA, NA)
    )
    expect_equal(
        match_gene_symbol(
            c("CD8", "CD4", "CD45"),
            c("ENSG00000081237.19", "ENSG00000153563.15", "ENSG00000010610.9", "ENSG00000288825")
        ),
        c(2, 3, 1)
    )
})