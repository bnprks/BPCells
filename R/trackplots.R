# Design: The aim is to make genome-browser style plots,
# where stacked plots are aligned to share a horizontal axis representing
# genome location.

# Each component plot function returns either a plain ggplot or a list of ggplot objects.
# These are then combined using draw_trackplot_grid, which handles alignment of
# all the plots with labels

#' Combine ggplot track plots into an aligned grid.
#' Uses patchwork to perform the alignment
#' @param ... Plots in order from top to bottom, generally plain ggplots.
#'    To better accomodate many bulk tracks, patchwork objects which contain multiple
#'    tracks are also accepted. In this case, plot labels will be drawn from the
#'    attribute `$patchwork$labels` if present, rather than the `labels` argument.
#' @param labels Text labels to display for each track
#' @param title Text for overarching title of the plot
#' @param heights Relative heights for each component plot. It is suggested to use 1 as standard height of a
#'    pseudobulk track.
#' @param label_width Fraction of width that should be used for labels relative to the main track area
#' @param label_style Arguments to pass to geom_text to adjust label text style
#' @return A plot object with aligned genome plots. Each aligned row has
#'    the text label, y-axis, and plot body. The relative height of each row is given
#'    by heights. A shared title and x-axis are put at the top.
#' @export
draw_trackplot_grid <- function(..., labels, title = NULL,
                                heights = rep(1, length(plots)),
                                label_width = 0.2,
                                label_style = list(fontface = "bold", size = 4)) {
  plots <- list(...)
  assert_true(length(plots) == length(labels))
  assert_true(length(plots) == length(heights))
  for (plot in plots) {
    assert_is(plot, "ggplot")
  }

  # Handle lists of plots by flattening appropriately
  # This is a bit awkward, but will hopefully save our function callers from inconvenience
  list_length <- vapply(plots, function(p) {
    if (is(p, "patchwork")) length(p$patches$plots) + 1 else 0
  }, 0)

  # Recalculate the labels and heights needed for the requested layout
  if (any(list_length != 0)) {
    total_elements <- sum(list_length) + sum(list_length == 0)
    next_index <- 1

    new_plots <- list()
    new_labels <- rep("", total_elements)
    new_heights <- rep(0, total_elements)

    for (i in seq_along(plots)) {
      if (list_length[i] == 0) {
        new_plots[[next_index]] <- plots[[i]]
        new_labels[next_index] <- labels[i]
        new_heights[next_index] <- heights[i]
        next_index <- next_index + 1
      } else {
        plot_list <- patchwork:::get_patches(plots[[i]])$plots
        names(plot_list) <- plots[[i]]$patchwork$labels
        if (is.null(names(plot_list))) {
          names(plot_list) <- rep("", length(plot_list))
        }
        for (j in seq_along(plot_list)) {
          new_plots[[next_index]] <- plot_list[[j]]
          new_labels[next_index] <- names(plot_list)[j]
          new_heights[next_index] <- heights[i]
          next_index <- next_index + 1
        }
      }
    }
    plots <- new_plots
    labels <- new_labels
    heights <- new_heights
  }

  x_axis_height <- 0.5

  # Make column of plot labels
  labels_plots <- lapply(labels, function(label_text) {
    ggplot2::ggplot(NULL, ggplot2::aes(0, 0, label = label_text)) +
      do.call(ggplot2::geom_text, label_style) +
      ggplot2::theme_void()
  })

  data_plots <- lapply(plots, function(p) {
    p + ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.ticks.length.x = ggplot2::unit(0, "pt"),
      axis.text.x = ggplot2::element_blank()
    )
  })
  data_plots[[1]] <- plots[[1]]

  patch <- Reduce(`+`, c(labels_plots, data_plots)) +
    patchwork::plot_layout(ncol = 2, byrow = FALSE, widths = c(label_width, 1), heights = heights, guides = "collect")
  if (!is.null(title)) {
    patch <- patch + patchwork::plot_annotation(title = title, theme = theme(plot.title = element_text(hjust = 0.5)))
  }
  return(patch)
}

#' Pseudobulk trackplot
#'
#' Plot a pseudobulk genome track, showing the number of fragment insertions
#' across a region.
#'
#' @inheritParams cluster_membership_matrix
#' @inheritParams plot_embedding
#' @inheritParams convert_to_fragments
#' @param region GRanges of length 1 with region to plot, or list/data.frame with
#'    one entry each for chr, start, end. See `gene_region()` or [genomic-ranges] for details
#' @param fragments Fragments object
#' @param cell_read_counts Numeric vector of read counts for each cell (used for normalization)
#' @param bins Number of bins to plot across the region
#' @param clip_quantile (optional) Quantile of values for clipping y-axis limits. Default of 0.999 will crop out
#'    just the most extreme outliers across the region. NULL to disable clipping
#' @param colors Character vector of color values (optionally named by group)
#' @param legend_label Custom label to put on the legend
#'
#' @return Returns a combined plot of pseudobulk genome tracks. For compatability with
#' `draw_trackplot_grid()`, the extra attribute `$patches$labels` will be added to
#' specify the labels for each track. If `return_data` or `return_plot_list` is
#' `TRUE`, the return value will be modified accordingly.
#' @export
trackplot_bulk <- function(fragments, region, groups,
                           cell_read_counts,
                           group_order = NULL,
                           bins = 200, clip_quantile = 0.999,
                           colors = discrete_palette("stallion"),
                           legend_label = "group",
                           zero_based_coords = !is(region, "GRanges"),
                           return_data = FALSE, return_plot_list = FALSE, apply_styling = TRUE) {
  assert_is(fragments, "IterableFragments")
  assert_not_null(cellNames(fragments))
  region <- normalize_ranges(region)
  assert_true(length(region$chr) == 1)
  assert_is_wholenumber(bins)
  if (!is.null(clip_quantile)) {
    assert_is_numeric(clip_quantile)
    assert_len(clip_quantile, 1)
  }
  assert_is(colors, "character")
  assert_true(length(colors) >= length(unique(groups)))
  assert_is(legend_label, "character")

  groups <- as.factor(groups)
  assert_true(length(cellNames(fragments)) == length(groups))
  assert_true(length(cellNames(fragments)) == length(cell_read_counts))


  region$tile_width <- (region$end - region$start) %/% bins

  membership_matrix <- cluster_membership_matrix(groups, group_order)
  group_read_counts <- multiply_rows(membership_matrix, cell_read_counts) %>%
    colSums()
  group_norm_factors <- 1e9 / (group_read_counts * region$tile_width)

  if (is.null(names(colors))) {
    names(colors) <- colnames(membership_matrix)
  }
  colors <- colors[seq_len(ncol(membership_matrix))]

  bin_centers <- seq(region$start, region$end - 1, region$tile_width) + region$tile_width / 2
  bin_centers <- pmin(bin_centers, region$end - 1)

  mat <- (tile_matrix(fragments, region, explicit_tile_names = TRUE) %*% membership_matrix) %>%
    as("dgCMatrix") %>%
    as("matrix")
  # Discard any partial bins
  mat <- mat[seq_along(bin_centers), ]


  data <- tibble::tibble(
    pos = rep(bin_centers, ncol(mat)),
    group = factor(rep(colnames(mat), each = nrow(mat)), levels=levels(groups)),
    insertions = as.vector(mat),
    # Normalized insertions are "IPKM" - insertions per kilobase per million mapped reads
    normalized_insertions = as.vector(multiply_cols(mat, group_norm_factors))
  )
  if (return_data) {
    return(data)
  }

  ymax <- quantile(data$normalized_insertions, clip_quantile)
  # Set precision of y-axis label to within 1% of the max value
  ymax_accuracy <- 10^as.integer(log10(0.01 * ymax))
  names(ymax) <- scales::label_comma(accuracy = ymax_accuracy)(ymax)

  data$normalized_insertions <- pmin(data$normalized_insertions, ymax)

  trackplots <- list()

  for (i in seq_len(ncol(membership_matrix))) {
    group <- colnames(membership_matrix)[i]

    plot <- dplyr::filter(data, group == .env$group) %>%
      ggplot2::ggplot(ggplot2::aes(pos, pmin(normalized_insertions, ymax), fill = group)) +
      ggplot2::geom_area() +
      ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
      ggplot2::scale_x_continuous(limits = c(region$start, region$end), expand = c(0, 0), position = "top") +
      ggplot2::scale_y_continuous(breaks = ymax, limits = c(0, ymax), expand = c(0, 0))
    if (apply_styling) {
      plot <- plot +
        ggplot2::labs(x = NULL, y = NULL, fill = legend_label) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank(),
          plot.margin = ggplot2::unit(c(0, 0, 0, 0), "pt")
        )
    }
    trackplots[[group]] <- plot
  }
  if (return_plot_list) {
    return(trackplots)
  } else {
    # Clear out axis tickmarks
    for (i in seq_along(trackplots)) {
      if (i != 1) {
        trackplots[[i]] <- trackplots[[i]] + ggplot2::theme(
          axis.title.x = ggplot2::element_blank(),
          axis.ticks.length.x = ggplot2::unit(0, "pt"),
          axis.text.x = ggplot2::element_blank()
        )
      }
    }
    plot <- patchwork::wrap_plots(trackplots, ncol = 1, guides = "collect")
    plot$patchwork$labels <- names(trackplots)
    return(plot)
  }
}


#' Plot transcript models
#' @param transcripts GRanges, list, or data.frame of transcript features to plot.
#' Required attributes are:
#' - `chr`, `start`, `end`: genomic position
#' - `strand`: "+"/"-" or TRUE/FALSE for positive or negative strand
#' - `feature` (only entries marked as `"transcript"` or `"exon"` will be considered)
#' - `transcript_id`
#' - `gene_name`
#' See [genomic-ranges] for more details
#' @inheritParams trackplot_bulk
#' @param labels Character vector with labels for each item in transcripts. NA for items that should not be labeled
#' @param exon_size size for exon lines
#' @param transcript_size size for transcript lines
#' @param label_size size for transcript labels
#' @return Plot of gene locations
#' @export
trackplot_gene <- function(transcripts, region, exon_size = 3, gene_size = 1, label_size = 3, return_data = FALSE) {
  region <- normalize_ranges(region)
  transcripts <- normalize_ranges(transcripts, metadata_cols = c("strand", "feature", "transcript_id", "gene_name"))
  data <- transcripts %>%
    tibble::as_tibble() %>%
    dplyr::filter(as.character(chr) == as.character(region$chr), end > region$start, start < region$end, feature %in% c("transcript", "exon")) %>%
    dplyr::mutate(size = dplyr::if_else(feature == "exon", exon_size, gene_size))

  if (return_data) {
    return(data)
  }

  ############################
  # Calculate y-position layout
  #############################
  # Steps:
  # 1. Calculate the maximum overlap depth of transcripts
  # 2. Iterate through transcript start/end in sorted order
  # 3. Randomly assign each transcript a y-coordinate between 1 and max overlap depth,
  #    with the restriction that a transcript can't have the same y-coordinate
  #    as a transcript it overlaps.
  tx_boundaries <- data %>% # Data frame of pos, tx_id, is_start
    dplyr::filter(feature == "transcript") %>%
    {
      dplyr::bind_rows(
        dplyr::transmute(., pos = start, transcript_id, start = TRUE),
        dplyr::transmute(., pos = end, transcript_id, start = FALSE)
      )
    } %>%
    dplyr::arrange(pos)
  total_positions <- max(cumsum(2 * tx_boundaries$start - 1))
  y_pos <- integer(0) # Names: tx_id, value: y position for transcript
  occupied_y <- rep(FALSE, total_positions) # Boolean vector, marking which y positions are open
  prev_seed <- get_seed()
  set.seed(12057235)
  for (i in seq_len(nrow(tx_boundaries))) {
    row <- tx_boundaries[i, ]
    if (row$start) {
      assigned_pos <- sample(which(!occupied_y), 1)
      y_pos[row$transcript_id] <- assigned_pos
      occupied_y[assigned_pos] <- TRUE
    } else {
      assigned_pos <- y_pos[row$transcript_id]
      occupied_y[assigned_pos] <- FALSE
    }
  }
  restore_seed(prev_seed)

  # Adjust the endpoints of any partially overlapping elements to fit within
  # the plot boundaries
  data <- dplyr::mutate(
    data,
    start = pmax(region$start, pmin(region$end, start)),
    end = pmax(region$start, pmin(region$end, end))
  )
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = dplyr::if_else(strand, start, end), xend = dplyr::if_else(strand, end, start),
      y = y_pos[transcript_id], yend = y_pos[transcript_id],
      size = size
    )
  ) +
    ggplot2::geom_segment(ggplot2::aes(color = dplyr::if_else(strand, "+", "-"))) +
    ggrepel::geom_text_repel(
      data = dplyr::filter(data, feature == "transcript"),
      ggplot2::aes(label = gene_name),
      size = label_size,
      position = ggrepel::position_nudge_repel(y = 0.25)
    ) +
    ggplot2::scale_size(range = c(min(exon_size, gene_size, label_size), max(exon_size, gene_size, label_size))) +
    ggplot2::scale_color_manual(values = c("+" = "black", "-" = "darkgrey")) +
    ggplot2::scale_x_continuous(limits = c(region$start, region$end), expand = c(0, 0)) +
    ggplot2::scale_y_discrete(labels = NULL, breaks = NULL) +
    ggplot2::labs(x = NULL, y = NULL, color = "strand") +
    ggplot2::guides(size = "none") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

  return(plot)
}
