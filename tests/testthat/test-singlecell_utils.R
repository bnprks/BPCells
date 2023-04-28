test_that("Wilcoxon rank sum works (small)", {
    x <- c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46)
    y <- c(1.15, 0.88, 0.90, 0.74, 1.21)

    med_xy <- median(c(x,y))
    x <- x - med_xy
    y <- y - med_xy

    res_r <- wilcox.test(x, y, exact=FALSE, correct=TRUE)

    vals <- c(x,y)
    X <- matrix(c(vals - median(vals), vals), ncol=2)
    colnames(X) <- c("centered", "uncentered")
    Y <- factor(c(rep.int("X", length(x)), rep.int("Y", length(y))))

    m <- as(X, "dgCMatrix") %>% 
        as("IterableMatrix") %>%
        t()

    res_bp_1 <- marker_features(m, Y) %>%
        dplyr::filter(foreground == "X", feature == "centered")
    res_bp_2 <- marker_features(m, Y) %>%
        dplyr::filter(foreground == "Y", feature == "uncentered")
    
    expect_equal(res_bp_1$p_val_raw, res_r$p.value)
    expect_equal(res_bp_1$p_val_raw, res_bp_2$p_val_raw)

    expect_equal(res_bp_1$foreground_mean, mean(x))
    expect_equal(res_bp_1$background_mean, mean(y))

    expect_equal(res_bp_2$foreground_mean, mean(y))
    expect_equal(res_bp_2$background_mean, mean(x))
    
    expect_equal(res_bp_1$feature, "centered")
})

test_that("Wilcoxon rank sum matches immunogenomics::presto", {
    skip_if_not_installed("pbmc3k.SeuratData")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("presto")
    pbmc <- pbmc3k.SeuratData::pbmc3k.final
    res_presto <- presto::wilcoxauc(pbmc@assays$RNA@data, Seurat::Idents(pbmc)) %>%
        tibble::as_tibble()
    
    mat_bpcells <- pbmc@assays$RNA@data %>% 
        t() %>% 
        as("IterableMatrix") %>% 
        t()
    res_bpcells <- marker_features(mat_bpcells, Seurat::Idents(pbmc)) %>%
        dplyr::arrange(match(foreground, levels(Seurat::Idents(pbmc))), match(feature, rownames(pbmc)))
    
    expect_equal(res_bpcells$feature, res_presto$feature)
    expect_equal(res_bpcells$foreground, res_presto$group)
    expect_equal(res_bpcells$p_val_raw, res_presto$pval)

    expect_equal(res_bpcells$foreground_mean, res_presto$avgExpr)
    expect_equal(res_bpcells$background_mean, res_presto$avgExpr - res_presto$logFC)

})
