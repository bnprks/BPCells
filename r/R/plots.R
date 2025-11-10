# Copyright 2022 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

#' Color palettes
#'
#' These color palettes are derived from the ArchR color palettes, and provide
#' large sets of distinguishable colors
#'
#' If the requested number of colors is too large, a new palette will be constructed
#' via interpolation from the requested palette
#'
#' @param name Name of the color palette. Valid discrete palettes are: `stallion`, `calm`, `kelly`, `bear`,
#'             `ironMan`, `circus`, `paired`, `grove`, `summerNight`, and `captain`. Valid continuous palettes
#'             are `bluePurpleDark`
#' @param n Minimum number of colors needed
#' @return Character vector of hex color codes
#' @rdname palettes
#' @export
discrete_palette <- function(name, n = 1) {
  palettes <- list(
    # 20-colors
    stallion = c(
      "#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#E6C2DC",
      "#90D5E4", "#89C75F", "#F37B7D", "#9983BD", "#D24B27", "#3BBCA8", "#6E4B9E", "#0C727C", "#7E1416", "#D8A767", "#3D3D3D"
    ),
    calm = c(
      "#7DD06F", "#844081", "#688EC1", "#C17E73", "#484125", "#6CD3A7", "#597873", "#7B6FD0", "#CF4A31", "#D0CD47",
      "#722A2D", "#CBC594", "#D19EC4", "#5A7E36", "#D4477D", "#403552", "#76D73C", "#96CED5", "#CE54D1", "#C48736"
    ),
    kelly = c(
      "#FFB300", "#803E75", "#FF6800", "#A6BDD7", "#C10020", "#CEA262", "#817066", "#007D34", "#F6768E", "#00538A",
      "#FF7A5C", "#53377A", "#FF8E00", "#B32851", "#F4C800", "#7F180D", "#93AA00", "#593315", "#F13A13", "#232C16"
    ),

    # 16-colors
    bear = c(
      "#faa818", "#41a30d", "#fbdf72", "#367d7d", "#d33502", "#6ebcbc", "#37526d",
      "#916848", "#f5b390", "#342739", "#bed678", "#a6d9ee", "#0d74b6",
      "#60824f", "#725ca5", "#e0598b"
    ),

    # 15-colors
    ironMan = c(
      "#371377", "#7700FF", "#9E0142", "#FF0080", "#DC494C", "#F88D51", "#FAD510", "#FFFF5F", "#88CFA4",
      "#238B45", "#02401B", "#0AD7D3", "#046C9A", "#A2A475", "grey35"
    ),
    circus = c(
      "#D52126", "#88CCEE", "#FEE52C", "#117733", "#CC61B0", "#99C945", "#2F8AC4", "#332288",
      "#E68316", "#661101", "#F97B72", "#DDCC77", "#11A579", "#89288F", "#E73F74"
    ),

    # 12-colors
    paired = c(
      "#A6CDE2", "#1E78B4", "#74C476", "#34A047", "#F59899", "#E11E26",
      "#FCBF6E", "#F47E1F", "#CAB2D6", "#6A3E98", "#FAF39B", "#B15928"
    ),

    # 11-colors
    grove = c("#1a1334", "#01545a", "#017351", "#03c383", "#aad962", "#fbbf45", "#ef6a32", "#ed0345", "#a12a5e", "#710162", "#3B9AB2"),

    # 10-colors
    
    # Swap Tableau's 2+3 in the ordering for aesthetic reasons
    tableau = c("#4E79A7", "#E15759", "#F28E2B", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"),
    
    # 7-colors
    summerNight = c("#2a7185", "#a64027", "#fbdf72", "#60824f", "#9cdff0", "#022336", "#725ca5"),

    # 5-colors
    captain = c("grey", "#A1CDE1", "#12477C", "#EC9274", "#67001E")
  )
  if (!(name %in% names(palettes))) {
    stop(sprintf("palette name must be one of: %s", toString(names(palettes))))
  }
  palette <- palettes[[name]]
  if (n > length(palette)) {
    palette <- grDevices::colorRampPalette(palette, space = "Lab")(n)
  }
  return(palette)
}

#' @rdname palettes
#' @export
continuous_palette <- function(name) {
  palettes <- list(
    "bluePurpleDark" = c(RColorBrewer::brewer.pal(9, "BuPu")[-(1:2)], "#000000")
  )
  if (!(name %in% names(palettes))) {
    stop(sprintf("palette name must be one of: %s", toString(names(palettes))))
  }
  palette <- palettes[[name]]
  return(palette)
}

#' Knee plot of single cell read counts
#'
#' Plots read count rank vs. number of reads on a log-log scale.
#'
#' Performs logarithmic downsampling to reduce the number of points plotted
#'
#' @param read_counts Vector of read counts per cell
#' @param cutoff (optional) Read cutoff to mark on the plot
#' @param return_data If true, return data from just before plotting rather than a plot.
#' @param apply_styling If false, return a plot without pretty styling applied
#' @return ggplot2 plot object
#' @examples
#' ## Prep data
#' mat <- get_demo_mat(filter_qc = FALSE, subset = FALSE)
#' reads_per_cell <- colSums(mat)
#' 
#' # Render knee plot
#' plot_read_count_knee(reads_per_cell, cutoff = 1e3)
#' @export
plot_read_count_knee <- function(read_counts, cutoff = NULL, return_data = FALSE, apply_styling = TRUE) {
  # Keep ~1k cells per order of magnitude to speed up plotting
  ranks <- unique(floor(10^(1:(1000 * ceiling(log10(length(read_counts)))) / 1000)))

  data <- list(
    data = tibble::tibble(
      ranks = ranks[ranks < length(read_counts)],
      reads = sort(read_counts, decreasing = TRUE)[ranks]
    ),
    cells_passing = NA
  )

  if (!is.null(cutoff)) data$cells_passing <- sum(read_counts > cutoff)

  if (return_data) {
    return(data)
  }


  plot <- ggplot2::ggplot(data$data, ggplot2::aes(log10(ranks), log10(reads))) +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(labels = scales::label_math(), breaks = scales::breaks_width(1)) +
    ggplot2::scale_y_continuous(labels = scales::label_math(), breaks = scales::breaks_width(1))

  plot <- plot + ggplot2::labs(x = "Barcode Rank", y = "Reads")

  if (!apply_styling) {
    return(plot)
  }

  # Add shaded rectangle and passing cells label
  if (!is.null(cutoff)) {
    cell_label <- tibble::tibble(
      label = sprintf("%s cells", scales::label_comma()(data$cells_passing)),
      x = max(log10(data$data$ranks)),
      y = max(log10(data$data$reads))
    )
    rectangle_highlight <- tibble::tibble(
      xmin = -Inf, xmax = Inf,
      ymin = log10(cutoff), ymax = Inf
    )
    plot <- plot +
      ggplot2::geom_hline(yintercept = log10(cutoff), linetype = "dashed") +
      ggplot2::geom_text(data = cell_label, ggplot2::aes(x, y, label = label), hjust = "inward", vjust = "inward") +
      ggplot2::geom_rect(
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, x = NULL, y = NULL),
        data = rectangle_highlight,
        alpha = 0.1
      )
  }

  plot <- plot +
    ggplot2::annotation_logticks() +
    ggplot2::theme_classic()

  return(plot)
}

#' TSS Enrichment vs. Fragment Counts plot
#'
#' Density scatter plot with log10(fragment_count) on the x-axis and TSS
#' enrichment on the y-axis. This plot is most useful to select which cell barcodes
#' in an experiment correspond to high-quality cells
#'
#' @param atac_qc Tibble as returned by `qc_scATAC()`. Must have columns `nFrags` and `TSSEnrichment`
#' @param min_frags Minimum fragment count cutoff
#' @param min_tss Minimum TSS Enrichment cutoff
#' @param bins Number of bins for density calculation
#' @inheritParams plot_embedding
#' @examples
#' ## Prep data
#' frags <- get_demo_frags(filter_qc = FALSE, subset = FALSE)
#' genes <- read_gencode_transcripts(
#'   file.path(tempdir(), "references"), release = "42",
#'   annotation_set = "basic",
#'   features = "transcript"
#' )
#' blacklist <- read_encode_blacklist(file.path(tempdir(), "references"), genome="hg38")
#' atac_qc <- qc_scATAC(frags, genes, blacklist)
#' 
#' 
#' ## Render tss enrichment vs fragment plot
#' plot_tss_scatter(atac_qc, min_frags = 1000, min_tss = 10)
#' @export
plot_tss_scatter <- function(atac_qc, min_frags = NULL, min_tss = NULL, bins = 100, apply_styling = TRUE) {
  assert_is(atac_qc, "data.frame")
  if (!is.null(min_frags)) {
    assert_is_numeric(min_frags)
  }
  if (!is.null(min_tss)) {
    assert_is_numeric(min_tss)
  }
  assert_is(apply_styling, "logical")

  assert_has_names(atac_qc, c("nFrags", "TSSEnrichment"))
  plot <- ggplot2::ggplot(atac_qc, ggplot2::aes(log10(nFrags), TSSEnrichment)) +
    ggplot2::stat_bin_hex(bins = bins)

  passing_cells <- rep(TRUE, nrow(atac_qc))
  if (!is.null(min_frags)) {
    plot <- plot + ggplot2::geom_vline(xintercept = log10(min_frags), linetype = "dashed")
    passing_cells <- passing_cells & atac_qc$nFrags > min_frags
  }
  if (!is.null(min_tss)) {
    plot <- plot + ggplot2::geom_hline(yintercept = min_tss, linetype = "dashed")
    passing_cells <- passing_cells & atac_qc$TSSEnrichment > min_tss
  }

  if (!apply_styling) {
    return(plot)
  }

  if (!is.null(min_frags) || !is.null(min_tss)) {
    cell_label <- tibble::tibble(
      label = sprintf("%s cells", scales::label_comma()(sum(passing_cells))),
      x = max(log10(atac_qc$nFrags)),
      y = max(atac_qc$TSSEnrichment)
    )
    if (is.null(min_frags)) {
      min_frags <- 0
    }
    rectangle_highlight <- tibble::tibble(
      xmin = max(-Inf, log10(min_frags)), xmax = Inf,
      ymin = max(-Inf, min_tss), ymax = Inf
    )
    plot <- plot +
      ggplot2::geom_text(data = cell_label, ggplot2::aes(x, y, label = label), hjust = "inward", vjust = "inward") +
      ggplot2::geom_rect(
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, x = NULL, y = NULL),
        data = rectangle_highlight,
        alpha = 0.1
      )
  }

  passing_cells <- sum(atac_qc$nFrags > min_frags)

  plot <- plot +
    ggplot2::scale_x_continuous(labels = scales::label_math(), breaks = scales::breaks_width(1)) +
    ggplot2::scale_fill_distiller(
      palette = "Spectral",
      trans = "log10",
      labels = scales::label_comma()
    ) +
    ggplot2::labs(x = "Fragments", y = "TSS Enrichment", fill = "Cell Density") +
    ggplot2::annotation_logticks(side = "b") +
    ggplot2::theme_classic()
  return(plot)
}

#' Fragment size distribution
#'
#' Plot the distribution of fragment lengths, with length in basepairs on the
#' x-axis, and proportion of fragments on the y-axis. Typical plots will show
#' 10-basepair periodicity, as well as humps spaced at multiples of a nucleosome
#' width (about 150bp).
#'
#' @param fragments Fragments object
#' @param max_length Maximum length to show on the plot
#' @inheritParams plot_embedding
#' @return Numeric vector where index i contans the number of length-i fragments
#' @examples
#' frags <- get_demo_frags(filter_qc = FALSE, subset = FALSE)
#' plot_fragment_length(frags)
#' @export
plot_fragment_length <- function(fragments, max_length = 500, return_data = FALSE, apply_styling = TRUE) {
  assert_is(fragments, "IterableFragments")
  assert_is_wholenumber(max_length)

  iter <- iterate_fragments(fragments)
  res <- fragment_lengths_cpp(iter)

  data <- tibble::tibble(
    length = seq_along(res),
    count = res,
    proportion = res / sum(res)
  )
  if (return_data) {
    return(data)
  }

  plot <- data %>%
    dplyr::filter(length <= max_length) %>%
    ggplot2::ggplot(ggplot2::aes(length, proportion)) +
    ggplot2::geom_line()
  if (!apply_styling) {
    return(plot)
  }
  plot <- plot +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(labels = scales::label_percent()) +
    ggplot2::labs(x = "Fragment Length (bp)", y = "Proportion")
  return(plot)
}

#' Plot TSS profile
#'
#' Plot the enrichmment of insertions relative to transcription start sites (TSS).
#' Typically, this plot shows strong enrichment of insertions near a TSS, and a
#' small bump downstream around 220bp downstream of the TSS for the +1 nucleosome.
#'
#' @inheritParams plot_embedding
#' @inheritParams footprint
#' @param genes Coordinate ranges for genes (must include strand)
#' @param smooth Number of bases to smooth over (rolling average)
#' @seealso `footprint()`, `plot_tf_footprint()`
#' @examples
#' ## Prep data
#' frags <- get_demo_frags(filter_qc = FALSE, subset = FALSE)
#' genes <- read_gencode_transcripts(
#'   file.path(tempdir(), "references"), release = "42",
#'   annotation_set = "basic",
#'   features = "transcript"
#' )
#' 
#' ## Plot tss profile
#' plot_tss_profile(frags, genes)
#' @export
plot_tss_profile <- function(fragments, genes, cell_groups = rlang::rep_along(cellNames(fragments), "all"),
                             flank = 2000L, smooth = 0L, zero_based_coords = !is(genes, "GRanges"),
                             colors = discrete_palette("stallion"),
                             return_data = FALSE, apply_styling = TRUE) {
  assert_is(fragments, "IterableFragments")
  assert_is_wholenumber(smooth)
  ranges <- normalize_ranges(genes, metadata_cols = "strand", zero_based_coords = zero_based_coords)
  data <- footprint(fragments, ranges, cell_groups = cell_groups, flank = flank + smooth)

  # Helper function to compute rolling means using cumsums
  rollmean <- function(x, k) {
    if (k <= 1) {
      return(x)
    }
    idx <- seq_along(x)
    roll_sum <- cumsum(x)[k - 1 + idx] - (cumsum(x) - x)[idx]
    roll_mean <- roll_sum / k
    na_left <- (k - 1) %/% 2
    return(c(rep.int(NA, na_left), roll_mean[seq_len(length(x) - na_left)]))
  }

  data <- data %>%
    dplyr::group_by(group) %>%
    dplyr::arrange(position) %>%
    dplyr::mutate(
      smoothed_enrichment = rollmean(enrichment, smooth)
    ) %>%
    dplyr::filter(abs(position) <= flank)

  if (return_data) {
    return(data)
  }
  plot <- ggplot2::ggplot(data, ggplot2::aes(position, smoothed_enrichment, color = group)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Distance from TSS (bp)", y = "Enrichment")
  if (!apply_styling) {
    return(plot)
  }
  group_name <- argument_name(cell_groups, 2)
  plot <- plot +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(color = group_name) +
    ggplot2::theme_classic()
  if (dplyr::n_distinct(data$group) == 1) {
    plot <- plot + ggplot2::guides(color = "none")
  }
  return(plot)
}

#' Plot TF footprint
#'
#' Plot the footprinting around TF motif sites
#'
#' @inheritParams plot_embedding
#' @inheritParams footprint
#' @param motif_positions Coordinate ranges for motifs (must include strand) and
#'   have constant width
#' @seealso `footprint()`, `plot_tss_profile()`
#' @export
plot_tf_footprint <- function(fragments, motif_positions, cell_groups = rlang::rep_along(cellNames(fragments), "all"),
                              flank = 250L, smooth = 0L, zero_based_coords = !is(genes, "GRanges"),
                              colors = discrete_palette("stallion"),
                              return_data = FALSE, apply_styling = TRUE) {
  assert_is(fragments, "IterableFragments")
  assert_is_wholenumber(smooth)
  ranges <- normalize_ranges(motif_positions, metadata_cols = "strand", zero_based_coords = zero_based_coords)
  if (1 != dplyr::n_distinct(ranges$end - ranges$start)) {
    rlang::abort("motif_positions must all have same width (end - start)")
  }
  width <- ranges$end[1] - ranges$start[1]
  ranges$start <- ranges$start + width %/% 2
  plot <- plot_tss_profile(fragments, ranges,
    cell_groups = cell_groups, flank = flank, smooth = smooth,
    colors = colors, return_data = return_data, apply_styling = apply_styling
  )
  if (return_data) {
    return(plot)
  }
  plot <- plot + ggplot2::labs(x = "Distance from Motif Center (bp)", y = "Enrichment")
  return(plot)
}

#' Collect features for plotting
#'
#' Helper function for data on features to plot from a diverse
#' set of data sources.
#'
#' If `source` is a data.frame, features will be drawn from the columns. If
#' `source` is a matrix object (`IterableMatrix`, `dgCMatrix`, or `matrix`), features
#' will be drawn from rows.
#'
#' @param source Matrix or data frame to pull features from, or a vector of feature values for a single feature. For
#'   a matrix, the features must be rows.
#' @param features Character vector of features names to plot if `source` is not a vector.
#' @param gene_mapping An optional vector for gene name matching with match_gene_symbol().
#'   Ignored if source is a data frame.
#' @param n Internal-use parameter marking the number of nested calls. This is used for
#'   finding the name of the "source" input variable from the caller's perspective
#' @return Data frame with one column for each feature requested
#' @export
collect_features <- function(source, features = NULL, gene_mapping = human_gene_mapping, n = 1) {
  if (!is.null(features)) {
    assert_is(features, "character")
    assert_is(source, c("data.frame", "IterableMatrix", "dgCMatrix", "matrix"))
  } else {
    if (!is.atomic(source) || is.null(source)) {
      rlang::abort("collect_features: if features is NULL, source must be a vector")
    }
  }
  if (is.null(features)) {
    data <- data.frame(source)
    colnames(data) <- argument_name(source, n + 1)
  } else if (is(source, "data.frame")) {
    # Data frame
    assert_has_names(source, features)
    data <- source[, features]
  } else if (is(source, "IterableMatrix") || is(source, "dgCMatrix") || is(source, "matrix")) {
    # Iterable Matrix
    assert_not_null(rownames(source))
    indices <- match(features, rownames(source), incomparables = NA)
    if (!is.null(gene_mapping)) {
      indices <- match_gene_symbol(features, rownames(source), gene_mapping)
    } else if (anyNA(indices)) {
      rlang::warn(sprintf(
        "collect_features could not match %d features: %s",
        sum(is.na(indices)),
        pretty_print_vector(features[is.na(indices)])
      ))
    }
    features <- features[!is.na(indices)]
    indices <- indices[!is.na(indices)]
    data <- source[indices, ]
    if (is(source, "IterableMatrix")) {
      data <- as(data, "dgCMatrix")
    }
    data <- data %>%
      as("matrix") %>%
      t()
    colnames(data) <- features
    rownames(data) <- colnames(source)
  }
  return(tibble::as_tibble(data, rownames = NA))
}

#' Plot UMAP or embeddings
#'
#' Plot one or more features by coloring cells in a UMAP plot.
#'
#' ### Smoothing
#' Smoothing is performed as follows: first, the smoothing matrix is normalized
#' so the sum of incoming weights to every cell is 1. Then, the raw data values
#' are repeatedly multiplied by the smoothing matrix and re-scaled so the average
#' value stays the same.
#'
#' @param source Matrix, or data frame to pull features from, or a vector of feature values for a single feature. For
#'   a matrix, the features must be rows.
#' @param embedding A matrix of dimensions cells x 2 with embedding coordinates
#' @param features Character vector of features to plot if `source` is not a vector.
#' @param quantile_range (optional) Length 2 vector giving the quantiles to clip the minimum and
#'   maximum color scale values, as fractions between 0 and 1. NULL or NA values to skip clipping
#' @param randomize_order If TRUE, shuffle cells to prevent overplotting biases. Can pass
#'   an integer instead to specify a random seed to use.
#' @param smooth (optional) Sparse matrix of dimensions cells x cells with cell-cell distance
#'   weights for smoothing.
#' @param smooth_rounds Number of multiplication rounds to apply when smoothing.
#' @param gene_mapping An optional vector for gene name matching with match_gene_symbol().
#'   Ignored if `source` is a data frame.
#' @param size Point size for plotting
#' @param rasterize Whether to rasterize the point drawing to speed up display in graphics
#' programs.
#' @param raster_pixels Number of pixels to use when rasterizing. Can provide one number for
#' square dimensions, or two numbers for width x height.
#' @param legend_continuous Whether to label continuous features by quantile or value. "auto"
#'   labels by quantile only when all features are continuous and `quantile_range` is not NULL.
#'   Quantile labeling adds text annotation listing the range of displayed values.
#' @param labels_quantile_range Whether to add a text label with the value range of each feature
#'   when the legend is set to quantile
#' @param colors_continuous Vector of colors to use for continuous color palette
#' @param labels_discrete Whether to add text labels at the center of each group for discrete (categorical) features.
#' @param legend_discrete Whether to show the legend for discrete (categorical) features.
#' @param colors_discrete Vector of colors to use for discrete (categorical) features.
#' @param return_plot_list If `TRUE`, return multiple plots as a list, rather than a
#'   single plot combined using patchwork::wrap_plots()
#' @inheritParams plot_read_count_knee
#' @return By default, returns a ggplot2 object with all the requested features plotted
#' in a grid. If `return_data` or `return_plot_list` is called, the return value will
#' match that argument.
#' @examples
## Prep data
#' set.seed(123)
#' mat <- get_demo_mat()
#' ## Normalize matrix
#' mat_norm <- log1p(multiply_cols(mat, 1/colSums(mat)) * 10000) %>% write_matrix_memory(compress = FALSE)
#' ## Get variable genes
#' stats <- matrix_stats(mat, row_stats = "variance")
#' variable_genes <- order(stats$row_stats["variance",], decreasing=TRUE) %>% 
#'   head(1000) %>% 
#'   sort()
#' # Z score normalize genes
#' mat_norm <- mat[variable_genes, ]
#' gene_means <- stats$row_stats['mean', variable_genes]
#' gene_vars <- stats$row_stats['variance', variable_genes]
#' mat_norm <- (mat_norm - gene_means) / gene_vars
#' ## Save matrix to memory
#' mat_norm <- mat_norm %>% write_matrix_memory(compress = FALSE)
#' ## Run SVD
#' svd <- BPCells::svds(mat_norm, k = 10)
#' pca <- multiply_cols(svd$v, svd$d)
#' ## Get UMAP
#' umap <- uwot::umap(pca)
#' ## Get clusters
#' clusts <- knn_hnsw(pca, ef = 500) %>%
#'   knn_to_snn_graph() %>%
#'   cluster_graph_louvain()
#' 
#' 
#' ## Plot embeddings
#' print(length(clusts))
#' 
#' plot_embedding(clusts, umap)
#' 
#' 
#' ### Can also plot by features
#' #plot_embedding(
#' #  source = mat,
#' #  umap,
#' #  features = c("MS4A1", "CD3E"),
#' #)
#' 
#' @export
plot_embedding <- function(source, embedding, features = NULL,
                           quantile_range = c(0.01, 0.99),
                           randomize_order = TRUE,
                           smooth = NULL,
                           smooth_rounds = 3,
                           gene_mapping = human_gene_mapping,
                           size = NULL,
                           rasterize = FALSE, raster_pixels = 512,
                           legend_continuous = c("auto", "quantile", "value"),
                           labels_quantile_range = TRUE,
                           colors_continuous = c("lightgrey", "#4682B4"),
                           legend_discrete = TRUE,
                           labels_discrete = TRUE,
                           colors_discrete = discrete_palette("stallion"),
                           return_data = FALSE, return_plot_list = FALSE, apply_styling = TRUE) {
  assert_is(embedding, "matrix")
  assert_true(ncol(embedding) == 2)
  assert_is(randomize_order, c("logical", "numeric"))
  assert_is_wholenumber(smooth_rounds)
  assert_is(raster_pixels, "numeric")
  legend_continuous <- match.arg(legend_continuous)
  assert_is(colors_continuous, "character")
  assert_is(labels_discrete, "logical")
  assert_is(legend_discrete, "logical")
  assert_is(colors_discrete, "character")
  assert_is(return_data, "logical")
  assert_is(return_plot_list, "logical")
  assert_is(apply_styling, "logical")

  # Check raster + size arguments
  if (is.null(size)) {
    size <- ifelse(rasterize, 1.01, 0.1)
  } else {
    assert_is_numeric(size)
  }

  if (!rasterize) {
    point_geom <- ggplot2::geom_point(size = size)
  } else {
    if (length(raster_pixels) == 1) {
      raster_pixels <- c(raster_pixels, raster_pixels)
    }
    point_geom <- scattermore::geom_scattermore(pointsize = size, pixels = raster_pixels)
  }

  # Fetch the actual data
  data <- list(embedding = embedding)
  data$data <- collect_features(source, features, gene_mapping, n = 2)

  # Check for mismatched cell counts/names
  if (is.null(source) && is.character(features) && length(features) != nrow(data$embedding) && length(features) < 50) {
    rlang::abort("Features length does not match number of cells in embedding.\nDid you forget to pass source?")
  }
  assert_true(nrow(data$data) == nrow(data$embedding))
  if (!is.null(rownames(data$embedding)) && !is.null(rownames(data$data))) {
    mapping <- match(rownames(data$data), rownames(data$embedding), incomparables = NA)
    if (!all(is.na(mapping))) {
      rlang::warn("Cell names are mismatched between source and embedding")
    } else {
      data$data <- data$data[mapping, ]
    }
  }

  # Perform data smoothing if requested
  if (!is.null(smooth)) {
    assert_is(smooth, "dgCMatrix")
    assert_true(nrow(smooth) == ncol(smooth))
    assert_true(nrow(smooth) == nrow(embedding))
    # Normalize smoothing inputs for each cell to sum to 1
    smooth <- t(t(smooth) / rowSums(smooth))

    numeric_cols <- as.logical(lapply(data$data, is.numeric))
    d <- as.matrix(data$data[, numeric_cols]) %>% t()
    target_means <- rowMeans(d)
    for (i in seq_len(smooth_rounds)) {
      d <- d %*% smooth
      d <- d * (target_means / rowMeans(d))
    }
    d <- t(d)
    for (i in seq_along(numeric_cols)) {
      if (numeric_cols[i]) {
        data$data[[i]] <- d[, i]
      }
    }
    rm(d, target_means)
  }

  # Perform quantile clipping if requested
  if (!is.null(quantile_range)) {
    assert_is(quantile_range, "numeric")
    data$quantile_cutoffs <- list()
    for (f in names(data$data)) {
      if (!is.numeric(data$data[[f]])) next()
      cutoffs <- quantile(data$data[[f]], quantile_range)
      data$quantile_cutoffs[[f]] <- cutoffs
      data$data[[f]] <- data$data[[f]] %>%
        pmax(min(cutoffs)) %>%
        pmin(max(cutoffs))
    }
  } else {
    if (legend_continuous == "quantile") {
      rlang::warn("legend_continuous = \"quantile\" not supported when quantile_range is NULL")
      legend_continuous <- "auto"
    }
  }

  # Generate a random order for cells
  if (is.numeric(randomize_order) || randomize_order) {
    if (is.numeric(randomize_order)) {
      # Set seed without permanently changing seed state
      prev_seed <- get_seed()
      set.seed(randomize_order)
    }
    data$random_order <- sample.int(nrow(embedding))
    if (is.numeric(randomize_order)) {
      restore_seed(prev_seed)
    }
  }

  if (return_data) {
    return(data)
  }

  # Apply random order for cells
  if (is.numeric(randomize_order) || randomize_order) {
    data$embedding <- data$embedding[data$random_order, ]
    data$data <- data$data[data$random_order, ]
  }

  # Find names for the embedding dimensions
  embedding_names <- colnames(embedding)
  if (is.null(embedding_names)) {
    embedding_names <- paste0(argument_name(embedding, 2), c(" 1", " 2"))
  }
  label_pos <- c(min(embedding[, 1]), max(embedding[, 2]))

  # Only share color scale if all data is numeric
  all_numeric <- all(as.logical(lapply(data$data, is.numeric)))
  if (legend_continuous == "auto") {
    if (all_numeric && !is.null(quantile_range) && ncol(data$data) > 1) {
      legend_continuous <- "quantile"
    } else {
      legend_continuous <- "value"
    }
  }

  aspect_ratio <- diff(range(data$embedding[, 1])) / diff(range(data$embedding[, 2]))

  # Make plots
  plots <- lapply(colnames(data$data), function(feature) {
    feature_values <- data$data[[feature]]
    # Re-scale if combining color scales
    is_continuous <- is.numeric(data$quantile_cutoffs[[feature]])
    if (is_continuous && legend_continuous == "quantile") {
      cutoffs <- data$quantile_cutoffs[[feature]]
      feature_values <- (feature_values - min(cutoffs)) / max(cutoffs)
    }

    plot <- ggplot2::ggplot(NULL, ggplot2::aes(data$embedding[, 1], data$embedding[, 2], color = feature_values)) +
      point_geom +
      ggplot2::labs(color = "")

    if (!is_continuous && labels_discrete) {
      centers <- tibble::tibble(
        x = data$embedding[, 1],
        y = data$embedding[, 2],
        label = feature_values
      ) %>%
        dplyr::group_by(label) %>%
        dplyr::summarize(x = mean(x), y = mean(y))
      plot <- plot + ggrepel::geom_text_repel(data = centers, mapping = ggplot2::aes(x, y, label = label, color = NULL))
    }

    if (!apply_styling) {
      return(plot)
    }

    # Handle color scale (both color and fill since we're trying to get the "rect" shape in the legend)
    if (is_continuous && legend_continuous == "quantile") {
      cutoffs <- data$quantile_cutoffs[[feature]]
      label <- sprintf("(%.2g-%.2g)", cutoffs[1], cutoffs[2])
      if (labels_quantile_range) {
        plot <- plot + ggplot2::labs(subtitle = label)
      }
      scale_labels <- sprintf("p%.2g", 100 * quantile_range)
      plot <- plot + ggplot2::scale_color_gradientn(
        colors = colors_continuous, breaks = c(0, 1), labels = scale_labels,
        limits = c(0, 1)
      )
    } else if (is_continuous) {
      plot <- plot + ggplot2::scale_color_gradientn(
        colors = colors_continuous,
        limits = c(0, 1)
      )
    } else {
      palette <- colors_discrete
      colors_needed <- dplyr::n_distinct(feature_values)
      if (length(palette) < colors_needed) {
        palette <- grDevices::colorRampPalette(palette, space = "Lab")(colors_needed)
      }
      plot <- plot + ggplot2::scale_color_manual(values = palette)
      if (legend_discrete) {
        # Set color guide to use squares, with a max of 8 per column
        plot <- plot + ggplot2::guides(
          color = ggplot2::guide_legend(nrow = 8, override.aes = list(size = 3, shape = "square"))
        )
      } else {
        plot <- plot + ggplot2::guides(color = "none")
      }
    }
    plot <- plot +
      ggplot2::labs(x = embedding_names[1], y = embedding_names[2], title = feature) +
      ggplot2::coord_fixed(aspect_ratio) +
      ggplot2::theme_classic() #+
    # ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))



    return(plot)
  })

  if (return_plot_list && length(plots) > 1) {
    return(plots)
  } else if (return_plot_list || length(plots) == 1) {
    return(plots[[1]])
  } else {
    guides <- ifelse(all_numeric && legend_continuous == "quantile", "collect", "auto")
    plot <- patchwork::wrap_plots(plots)
    if (apply_styling) {
      plot <- plot + patchwork::plot_layout(guides = guides, byrow = TRUE)
    }
    return(plot)
  }
}

#' Rotate ggplot x axis labels
#' @param degrees Number of degrees to rotate by
#' @export
rotate_x_labels <- function(degrees = 45) {
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = degrees, hjust = 1, vjust = 1))
}

#' Dotplot
#'
#' Plot feature levels per group or cluster as a grid of dots.
#' Dots are colored by z-score normalized average expression, and sized by percent non-zero.
#'
#' @param source Feature x cell matrix or data.frame with features. For best results,
#'   features should be sparse and log-normalized (e.g. run `log1p()` so zero raw counts
#'   map to zero)
#' @param features Character vector of features to plot
#' @inheritParams cluster_membership_matrix
#' @inheritParams plot_embedding
#' @param gene_mapping An optional vector for gene name matching with match_gene_symbol().
#' @param colors Color scale for plot
#' @examples
#' ## Prep data
#' mat <- get_demo_mat()
#' cell_types <- paste("Group", rep(1:3, length.out = length(colnames(mat))))
#' 
#' ## Plot dot
#' plot <- plot_dot(mat, c("MS4A1", "CD3E"), cell_types)
#' 
#' BPCells:::render_plot_from_storage(
#'   plot, width = 4, height = 5
#' )
#' @export
plot_dot <- function(source, features, groups, group_order = NULL, gene_mapping = human_gene_mapping,
                     colors = c("lightgrey", "#4682B4"),
                     return_data = FALSE, apply_styling = TRUE) {
  assert_is(features, "character")
  assert_is(source, c("IterableMatrix", "matrix", "dgCMatrix"))

  feature_data <- collect_features(source, features, gene_mapping, n = 2)

  numeric_features <- as.logical(lapply(feature_data, is.numeric))
  if (!all(numeric_features)) {
    rlang::abort(sprintf(
      "Non-numeric features detected: %s",
      pretty_print_vector(features[!numeric_features])
    ))
  }
  feature_data <- as.matrix(feature_data)
  group_weights <- cluster_membership_matrix(groups, group_order)
  group_weights <- multiply_cols(group_weights, 1 / colSums(group_weights))

  avg <- t(group_weights) %*% scale(feature_data)
  pct <- 100 * (t(group_weights) %*% (feature_data != 0))

  data <- tibble::tibble(
    group = rep(rownames(avg), ncol(avg)),
    feature = rep(colnames(avg), each = nrow(avg)),
    average = as.vector(avg),
    percent = as.vector(pct)
  )
  if (return_data) {
    return(data)
  }

  group_name <- argument_name(groups, 2)
  plot <- ggplot2::ggplot(data, ggplot2::aes(feature, group, color = average, size = percent)) +
    ggplot2::geom_point() +
    ggplot2::scale_size_area()
  if (!apply_styling) {
    return(plot)
  }
  plot <- plot + ggplot2::scale_color_gradientn(colors = colors) +
    ggplot2::labs(x = "Features", y = group_name, size = "% Detected", color = "Mean Z-score") +
    ggplot2::theme_classic() +
    ggplot2::scale_y_discrete(limits = rownames(avg)) +
    ggplot2::scale_x_discrete(limits = colnames(avg)) +
    rotate_x_labels(45)
  return(plot)
}
