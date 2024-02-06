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

test_that("C++ SNN calculation works",{
    
    k <- 10

    for (cells in c(100, 1000)) {
        neighbor_sim <- matrix(sample.int(cells, cells*k, replace=TRUE), nrow=cells)

        # Make sure no cell is listed twice as a nearest neighbor
        for (i in seq_len(cells)) {
            while (length(unique(neighbor_sim[i,])) != k) {
                neighbor_sim[i,] <- sample.int(cells, k)
            }
        }

        min_val <- 1/15
        snn <- knn_to_snn_graph(list(idx=neighbor_sim), min_val=min_val, self_loops=TRUE)

        mat <- knn_to_graph(list(idx=neighbor_sim), use_weights=FALSE)

        mat <- mat %*% t(mat)
        mat <- mat / (2 * k - mat)
        mat@x[mat@x < min_val] <- 0
        # Prune the explicit 0 entries from storage
        mat <- Matrix::drop0(mat) 
        mat <- Matrix::tril(mat)
        expect_identical(
            snn,
            as(mat, "dgCMatrix")
        )
        snn2 <- knn_to_snn_graph(list(idx=neighbor_sim), min_val=min_val, self_loops=FALSE)
        diag(mat) <- 0
        mat <- Matrix::drop0(mat)
        expect_identical(
            snn2,
            as(mat, "dgCMatrix")
        )
    }

})

test_that("C++ UMAP graph calculation works", {
    # This uses a pre-calculated graph from the umap-learn python package.
    # See ../data/generate_iris_geodesic_graph.R for the generation code 
    test_data <- readRDS("../data/iris_geodesic_graph.rds")
    knn <- test_data$knn
    res <- knn_to_geodesic_graph(knn)
    expect_equal(as.matrix(res + t(res)), as.matrix(test_data$graph), tolerance=1e-6)
})