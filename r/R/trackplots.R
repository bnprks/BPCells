# Copyright 2022 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

# Design: The aim is to make genome-browser style plots,
# where stacked plots are aligned to share a horizontal axis representing
# genome location.

# Each component plot function returns a faceted ggplot object, which can then be stacked
# via trackplot-combine

trackplot_theme <- function(base_size=11) {
  ggplot2::theme_bw(base_size=base_size) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    panel.spacing.y = ggplot2::unit(0, "pt"),
    strip.text.y.left = ggplot2::element_text(angle=0, hjust=1, size=ggplot2::rel(1.2)),
    strip.background = ggplot2::element_blank(),
    strip.placement = "outside",
    axis.title.y.left = ggplot2::element_text(size=ggplot2::rel(1))
  )
}

# Attributes on trackplots
# - height (length-1 unit): Height of trackplot. "null" units give relative sizes to other tracks, and other units give absolute sizes
# - takes_sideplot (length-1 bool): Whether this trackplot is an appropriate place to put a sideplot (for now just `trackplot_coverage()`)
# - region (list with chr, start, end): The genome region this plot covers.
# - keep_vertical_margin (length-1 bool): Whether this trackplot should retain vertical margin even when stacking.
#                               `trackplot_scalebar()` is the only FALSE case for now.
wrap_trackplot <- function(plot, height=NULL, takes_sideplot=FALSE, region=NULL, keep_vertical_margin=FALSE) {
  if (!is.null(height)) {
    assert_is(height, "unit")
  }
  if (!is.null(region)) {
    assert_has_names(region, c("chr", "start", "end"))
    region <- sprintf("%s:%d-%d", region[["chr"]], region[["start"]], region[["end"]])
  }
  assert_is(keep_vertical_margin, "logical")
  if (!is(plot, "trackplot")) {
    class(plot) <- c("trackplot", class(plot))
  }
  plot$trackplot <- list(height=height, takes_sideplot=takes_sideplot, region=region, keep_vertical_margin=keep_vertical_margin)
  plot
}

# Internal helper function to return empty track plots if there's no data to be plotted
trackplot_empty <- function(region, label) {
  ggplot2::ggplot(tibble::tibble(label=label)) +
    ggplot2::scale_x_continuous(limits = c(region$start, region$end), expand = c(0, 0), labels=scales::label_number()) +
    ggplot2::labs(x = "Genomic Position (bp)", y = NULL) +
    ggplot2::facet_wrap("label", strip.position="left") +
    trackplot_theme()
}

get_patchwork_plots <- function(patchwork) {
  assert_is(patchwork, "patchwork")
  ret <- plot$patches
  plot$patches <- NULL
  class(plot) <- setdiff(class(plot), "patchwork")
  c(ret, plot)
}

#' Adjust trackplot properties
#'
#' Adjust labels and heights on trackplots. Labels are set as facet labels in ggplot2, and heights
#' are additional properties read by `trackplot_combine()` to determine relative height of input plots.
#'
#' @param plot ggplot object
#' @param labels character vector of labels -- must match existing number of facets in plot
#' @return **set_trackplot_label**: ggplot object with adjusted facet labels
#' @rdname trackplot_utils
#' @export
set_trackplot_label <- function(plot, labels) {
  stopifnot(is(plot, "ggplot"))
  stopifnot(is(labels, "character"))
  
  labeller <- ggplot2::as_labeller(function(x) {
    if (length(x) != length(labels)) stop(sprintf("set_trackplot_label(): Plot has %d facets, but %d labels were provided", length(x), length(labels)))
    labels
  })
  if (is(plot$facet, "FacetNull")) {
    if (length(labels) != 1) stop(sprintf("Can only add single label to non-faceted plot. Got %d labels", length(labels)))
    return(plot + ggplot2::facet_wrap(NULL, strip.position="left", labeller=labeller))
  } else if (is(plot$facet, "FacetWrap")) {
    y <- plot
    y$facet <- ggplot2::ggproto(NULL, ggplot2::FacetWrap, shrink=y$facet$shrink, params=y$facet$params)
    y$facet$params$labeller <- labeller
    return(y)
  } else {
    stop("Unrecognized facet class: only facet_wrap() and facet_null() currently supported")
  }
}

#' @rdname trackplot_utils
#' @param height New height. If numeric, adjusts relative height. If `ggplot2::unit` or `grid::unit` sets absolute height in specified units.
#'     `"null"` units are interpreted as relative height.
#' @return **set_trackplot_height**: ggplot object with adjusted trackplot height
#' @export
set_trackplot_height <- function(plot, height) {
  if (!is(plot, "trackplot")) {
    plot <- wrap_trackplot(plot)
  }
  if (!is(height, "unit")) {
    height <- ggplot2::unit(height, "null")
  }
  plot$trackplot$height <- height
  plot
}

#' @rdname trackplot_utils
#' @return **get_trackplot_height**: `ggplot2::unit` object with height setting
#' @export
get_trackplot_height <- function(plot) {
  if (is(plot, "trackplot")) {
    plot$trackplot$height
  } else {
    ggplot2::unit(1L, "null")
  }
}

#' Calculate y positions for trackplot segments to avoid overlap
#' Steps:
#' 1. Calculate the maximum overlap depth of transcripts
#' 2. Iterate through start/end of segments in sorted order
#' 3. Randomly assign each segment a y-coordinate between 1 and max overlap depth,
#'   with the restriction that a segment can't have the same y-coordinate as an overlapping segment
#' @param data tibble of genome ranges with start and end columns, assumed to be on same chromosome. 
#' @return Vector of y coordinates, one per input row, such that no ranges at the same y coordinate overlap
#' @keywords internal
trackplot_calculate_segment_height <- function(data) {
  data$row_number <- seq_len(nrow(data))
  boundaries <- data %>%
    {
      dplyr::bind_rows(
        dplyr::mutate(., pos = start, row_number, start = TRUE, .keep="none"),
        dplyr::mutate(., pos = end, row_number, start = FALSE, .keep="none")
      )
    } %>%
    dplyr::arrange(pos)
  
  if (nrow(data) > 0) total_positions <- max(cumsum(2 * boundaries$start - 1))
  else total_positions <- 0L
  y_pos <- integer(nrow(data))
  occupied_y <- rep(FALSE, total_positions) # Boolean vector, marking which y positions are open
  prev_seed <- get_seed()
  set.seed(12057235)
  for (i in seq_len(nrow(boundaries))) {
    row <- boundaries[i, ]
    if (row$start) {
      # if else block to account for sample interpreting int of length 1 as a range
      assigned_pos <- which(!occupied_y)
      if (length(assigned_pos) > 1) assigned_pos <- sample(assigned_pos, 1)
      y_pos[row$row_number] <- assigned_pos
      occupied_y[assigned_pos] <- TRUE
    } else {
      assigned_pos <- y_pos[row$row_number]
      occupied_y[assigned_pos] <- FALSE
    }
  }
  restore_seed(prev_seed)
  return(y_pos)
}

#' Break up segments into smaller segments the length of the plot, divided by size
#' @param data Dataframe of full segments to be broken up
#' @param region Region to be plotted with end and start attr
#' @param size int Number of arrows to span the x axis of track
#' @param head_only bool If TRUE, only the head of the segment will be plotted
#' @return Dataframe of segments broken up into smaller segments.  Has columns start, end, and any additional metadata columns in original data
#' @keywords internal
trackplot_create_arrow_segs <- function(data, region, size = 50, head_only = FALSE) {
  # Get region to be plotted
  arrow_spacing <- (region$end - region$start) / size
  # Provide an initial arrow_list element to avoid schema errors if no transcripts are present
  arrow_list <- NULL
  for (i in seq_len(nrow(data))) {
    if (data$strand[i]) {
      endpoints <- seq(data$start[i], data$end[i], arrow_spacing)
    } else {
      endpoints <- rev(seq(data$end[i], data$start[i], -arrow_spacing))
    }
    new_arrow <- tibble::tibble(start = endpoints[-length(endpoints)], end = endpoints[-1])
    if (nrow(new_arrow) > 0) {
      arrow_list <- c(arrow_list, list(dplyr::bind_cols(
        new_arrow,
        data %>% dplyr::select(-c("start", "end")) %>% dplyr::slice(i)
      )))
    } else {
      arrow_list <- c(arrow_list, list(dplyr::slice(data, i)))
    }
  }
  arrows <- dplyr::bind_rows(arrow_list)
  if (head_only) {
    # Set segment size small enough to be invisible if head_only
    # Cannot swap values with if-else mutation block due to sequential evaluation
    arrows <- arrows %>% dplyr::mutate(
        start_tmp = start,
        start = ifelse(strand, end-1e-4, start),
        end = ifelse(strand, end, start_tmp+1e-4),
        start_tmp = NULL
    )
  }
  arrows <- dplyr::filter(arrows, start >= region$start, start < region$end, end >= region$start, end < region$end)
  return(arrows)
}

#' Normalize trackplot ranges data, while handling metadata argument renaming and type conversions
#' Type conversions are as follows:
#'    color -> factor or numeric
#'    label -> string
#' @param data Input ranges-like object
#' @param metadata List of form e.g. list(color=color_by, label=label_by). The values can be either column names or data vectors.
#'    Any NULL values will be skipped
#' @return Tibble with normalized ranges and additional columns populated as requested in `metadata`
#' @keywords internal
trackplot_normalize_ranges_with_metadata <- function(data, metadata) {
  metadata_column_names <- character(0)
  metadata_values <- list()
  # Check which columns need to be fetched in `normalize_ranges`
  for (i in seq_along(metadata)) {
    if (is.character(metadata[[i]]) && length(metadata[[i]]) == 1) {
      metadata_column_names <- union(metadata_column_names, metadata[[i]])
    } 
  }
  data <- normalize_ranges(data, metadata_cols = metadata_column_names)
  # Collect the metadata
  for (i in seq_along(metadata)) {
    if (is.null(metadata[[i]])) next

    col_name <- names(metadata)[i]
    if (is.character(metadata[[i]]) && length(metadata[[i]]) == 1) {
      metadata_values[[col_name]] <- data[[metadata[[i]]]]
    } else {
      if (length(metadata[[i]]) != nrow(data)) {
        rlang::abort(sprintf("Metadata for '%s' must match length of input data frame if given as a data vector", col_name))
      }
      metadata_values[[col_name]] <- metadata[[i]]
    }

    if (col_name == "label") metadata_values[[col_name]] <- as.character(metadata_values[[col_name]])
    if (col_name == "color") {
      if (!is.numeric(metadata_values[[col_name]])) metadata_values[[col_name]] <- as.factor(metadata_values[[col_name]]) 
    }
  }

  if (length(metadata_column_names) > 0) data <- dplyr::select(data, !dplyr::any_of(metadata_column_names))
  if (length(metadata_values) > 0) data <- dplyr::bind_cols(data, tibble::as_tibble(metadata_values))
  return(data)
}

#' Render a plot with intermediate disk storage step
#' 
#' Take a plotting object and save in temp storage, so it can be outputted with exact dimensions.
#' Primarily used to allow for adjusting plot dimensions within function reference examples.
#' @param plot (ggplot) ggplot output from a plotting function
#' @param width (numeric) width of rendered plot
#' @param height (numeric) height of rendered plot
#' @keywords internal
render_plot_from_storage <- function(plot, width, height) {
  assert_is(plot, "ggplot")
  image_path <- tempfile(fileext = ".png")
  ggplot2::ggsave(image_path, plot, width = width, height = height)
  img <- png::readPNG(image_path)
  grid::grid.raster(img)
}

#' Combine track plots
#' 
#' Combines multiple track plots of the same region into a single grid.
#' Uses the `patchwork` package to perform the alignment.
#'
#' @param tracks List of tracks in order from top to bottom, generally ggplots as output from
#'    the other `trackplot_*()` functions.
#' @param side_plot Optional plot to align to the right (e.g. RNA expression per cluster). Will be aligned to the first
#'    `trackplot_coverage()` output if present, or else the first generic ggplot in the alignment. Should be in horizontal orientation and 
#'    in the same cluster ordering as the coverage plots.
#' @param title Text for overarching title of the plot
#' @param side_plot_width Fraction of width that should be used for the side plot relative to the main track area
#' @return A plot object with aligned genome plots. Each aligned row has
#'    the text label, y-axis, and plot body. The relative height of each row is given
#'    by heights. A shared title and x-axis are put at the top.
#' @seealso `trackplot_coverage()`, `trackplot_gene()`, `trackplot_loop()`, `trackplot_scalebar()`
#' @examples
#' ## Prep data
#' frags <- get_demo_frags()
#' 
#' ## Use genes and blacklist to determine proper number of reads per cell
#' genes <- read_gencode_transcripts(
#'   file.path(tempdir(), "references"), release = "42",
#'   annotation_set = "basic",
#'   features = "transcript"
#' )
#' blacklist <- read_encode_blacklist(file.path(tempdir(), "references"), genome="hg38")
#' read_counts <- qc_scATAC(frags, genes, blacklist)$nFrags
#' region <- "chr4:3034877-4034877"
#' cell_types <- paste("Group", rep(1:3, length.out = length(cellNames(frags))))
#' transcripts <- read_gencode_transcripts(
#'   file.path(tempdir(), "references"), release = "42",
#'   annotation_set = "basic"
#' )
#' region <- "chr4:3034877-4034877"
#' 
#' 
#' ## Get all trackplots and scalebars to combine
#' plot_scalebar <- trackplot_scalebar(region)
#' plot_gene <- trackplot_gene(transcripts, region)
#' plot_coverage <- trackplot_coverage(frags, region, groups = cell_types, cell_read_counts = read_counts)
#' 
#' 
#' ## Combine trackplots and render
#' ## Also remove colors from gene track
#' plot <- trackplot_combine(
#'     list(plot_scalebar, plot_coverage, plot_gene + ggplot2::guides(color = "none"))
#' )
#' BPCells:::render_plot_from_storage(plot, width = 6, height = 4)
#' @export
trackplot_combine <- function(tracks, side_plot = NULL, title = NULL, side_plot_width = 0.3) {
  for (plot in tracks) {
    assert_is(plot, "ggplot")
  }
  if (!is.null(side_plot)) {
    assert_is(side_plot, "ggplot")
  }

  # Calculate layout information on the plots
  heights <- list()
  collapse_upper_margin <- rep.int(TRUE, length(tracks))
  side_plot_row <- NULL
  areas <- NULL
  last_region <- NULL
  for (i in seq_along(tracks)) {
    if (is(tracks[[i]], "trackplot")) {
      if (tracks[[i]]$trackplot$takes_sideplot && is.null(side_plot_row)) {
        side_plot_row <- i
      }
      # If we switch regions, don't collapse margins into the track above
      if (!is.null(tracks[[i]]$trackplot$region)) {
        if (!is.null(last_region) && last_region != tracks[[i]]$trackplot$region && i > 1) collapse_upper_margin[i] <- FALSE 
        last_region <- tracks[[i]]$trackplot$region
      }

      # Preserve top and bottom margins if `keep_vertical_margin`
      if (tracks[[i]]$trackplot$keep_vertical_margin) {
        collapse_upper_margin[i] <- FALSE
        if (i < length(tracks)) collapse_upper_margin[i+1] <- FALSE
      }
    } else {
      if (is.null(side_plot_row)) side_plot_row <- i
    }
    heights <- c(heights, list(get_trackplot_height(tracks[[i]])))
    areas <- c(areas, list(patchwork::area(i, 1)))
  }
  heights <- do.call(grid::unit.c, heights)
  if (!is.null(side_plot) && is.null(side_plot_row)) {
    rlang::warn("Did not find a row to place the side_plot: no trackplot_coverage() or base ggplot tracks found. Defaulting to first row")
    side_plot_row <- 1L
  }

  # Collapse margins as needed among plots
  for (i in seq_along(tracks)) {
    plot.margin <- c(TRUE, TRUE, TRUE, TRUE) # Top, right, bottom, left
    if (!is.null(side_plot)) plot.margin[2] <- FALSE # Side plot should be flush on right side
    if (i < length(tracks) && collapse_upper_margin[i+1]) plot.margin[3] <- FALSE # Plot below should be flush
    if (collapse_upper_margin[i]) plot.margin[1] <- FALSE

    if (!plot.margin[3]) {
      tracks[[i]] <- tracks[[i]] +
        ggplot2::guides(x="none") +
        ggplot2::labs(x=NULL)
    }

    # Independent of showing the axis, we'll remove the bottom margin if the next row has the side_plot, since the
    # axis tick labels will add in some natural margin already
    if (i+1 == side_plot_row) plot.margin[3] <- FALSE

    tracks[[i]] <- tracks[[i]] + ggplot2::theme(plot.margin=ggplot2::unit(5.5*plot.margin, "pt"))
  }

  # Reduce cut-off y-axis labels. Put plots with y axis labels later in the plot list, as they will take layer priority with patchwork
  has_y_axis <- vapply(tracks, function(t) is(t, "ggplot") && !is.null(t$labels$y), logical(1))
  tracks <- c(tracks[!has_y_axis], tracks[has_y_axis])
  areas <- c(areas[!has_y_axis], areas[has_y_axis])

  if (is.null(side_plot)) {
    widths <- c(1)
  } else {
    # Decide whether to put legends below/above side plot by adding up the height of all relatively-sized tracks
    height_above <- sum(as.vector(heights)[seq_along(heights) < side_plot_row & grid::unitType(heights) == "null"])
    height_below <- sum(as.vector(heights)[seq_along(heights) > side_plot_row & grid::unitType(heights) == "null"])
    if (height_above < height_below) {
      guide_position <- patchwork::area(side_plot_row+1L, 2, length(tracks))
    } else {
      guide_position <- patchwork::area(1L, 2, side_plot_row-1L)
    }
    
    widths <- c(1, side_plot_width)
    areas <- c(areas, list(patchwork::area(side_plot_row, 2), guide_position))
    # Make adjustments to the side plot style to fit in with tracks
    side_plot <- side_plot + 
      ggplot2::scale_x_discrete(limits=rev, position="top") +
      ggplot2::scale_y_continuous(position="right") +
      ggplot2::coord_flip() +
      ggplot2::labs(x=side_plot$labels$y, y=NULL) +
      ggplot2::theme(
        plot.margin=ggplot2::unit(c(0,0,0,0), "pt"),
        axis.ticks.length.x.top=grid::unit(-2.75, "pt")
      ) 
    guide_area <- patchwork::guide_area() + ggplot2::theme(plot.margin=ggplot2::unit(c(0,0,0,0), "pt"))
    tracks <- c(tracks, list(side_plot, guide_area))
  }

  patch <- patchwork::wrap_plots(tracks) +
    patchwork::plot_layout(ncol = 1, byrow = FALSE, heights = heights, widths=widths, guides = "collect", design=do.call(c, areas))
  
  if (!is.null(side_plot)) {
    # If a side plot is present, switch legend layout to use the horizontal space better
    patch <- patch * ggplot2::theme(
      legend.direction="horizontal", 
      legend.title.position = "top"
    )
  }
  if (!is.null(title)) {
    patch <- patch + patchwork::plot_annotation(title = title, theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  }
  return(patch)
}




#' Pseudobulk coverage trackplot
#'
#' Plot a pseudobulk genome track, showing the number of fragment insertions
#' across a region for each cell type or group. 
#'
#' @inheritParams cluster_membership_matrix
#' @inheritParams plot_embedding
#' @inheritParams convert_to_fragments
#' @param region Region to plot, e.g. output from `gene_region()`. String of format "chr1:100-200",
#'    or list/data.frame/GRanges of length 1 specifying chr, start, end. See `help("genomic-ranges-like")` for details
#' @param fragments Fragments object
#' @param cell_read_counts Numeric vector of read counts for each cell (used for normalization)
#' @param scale_bar Whether to include a scale bar in the top track (`TRUE` or `FALSE`)
#' @param bins Number of bins to plot across the region
#' @param clip_quantile (optional) Quantile of values for clipping y-axis limits. Default of 0.999 will crop out
#'    just the most extreme outliers across the region. NULL to disable clipping
#' @param colors Character vector of color values (optionally named by group)
#' @param legend_label `r lifecycle::badge("deprecated")` Custom label to put on the legend (no longer used as color legend is not shown anymore)
#'
#' @return Returns a combined plot of pseudobulk genome tracks. For compatability with
#' `draw_trackplot_grid()`, the extra attribute `$patches$labels` will be added to
#' specify the labels for each track. If `return_data` or `return_plot_list` is
#' `TRUE`, the return value will be modified accordingly.
#' @seealso `trackplot_combine()`, `trackplot_gene()`, `trackplot_loop()`, `trackplot_scalebar()`
#' @examples
## Prep data
#' frags <- get_demo_frags()
#' 
#' ## Use genes and blacklist to determine proper number of reads per cell
#' genes <- read_gencode_transcripts(
#'   file.path(tempdir(), "references"), release = "42",
#'   annotation_set = "basic",
#'   features = "transcript"
#' )
#' blacklist <- read_encode_blacklist(file.path(tempdir(), "references"), genome="hg38")
#' read_counts <- qc_scATAC(frags, genes, blacklist)$nFrags
#' region <- "chr4:3034877-4034877"
#' cell_types <- paste("Group", rep(1:3, length.out = length(cellNames(frags))))
#' 
#' 
#' BPCells:::render_plot_from_storage(
#'   trackplot_coverage(frags, region, groups = cell_types, cell_read_counts = read_counts),
#'   width = 6, height = 3
#' )
#' @export
trackplot_coverage <- function(fragments, region, groups,
                           cell_read_counts,
                           group_order = NULL,
                           bins = 500, clip_quantile = 0.999,
                           colors = discrete_palette("stallion"),
                           legend_label = NULL,
                           zero_based_coords = !is(region, "GRanges"),
                           return_data = FALSE) {
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
  if (!is.null(legend_label)) {
    lifecycle::deprecate_warn("0.3.0", "trackplot_coverage(legend_label)", details="Argument value is no longer used since color legend is not shown.")
  }

  groups <- as.factor(groups)
  assert_true(length(cellNames(fragments)) == length(groups))
  assert_true(length(cellNames(fragments)) == length(cell_read_counts))


  region$tile_width <- max(((region$end - region$start) %/% bins), 1)

  membership_matrix <- cluster_membership_matrix(groups, group_order)
  group_read_counts <- multiply_rows(membership_matrix, cell_read_counts) %>%
    colSums()
  group_norm_factors <- 1e9 / (group_read_counts * region$tile_width)

  if (is.null(names(colors))) {
    names(colors) <- colnames(membership_matrix)
  }
  colors <- colors[seq_len(ncol(membership_matrix))]

  bin_centers <- seq(region$start, region$end - 1, region$tile_width) + ((region$tile_width-1) / 2)
  bin_centers <- pmin(bin_centers, region$end - 1)

  mat <- (tile_matrix(fragments, region, explicit_tile_names = TRUE) %*% membership_matrix) %>%
    as("dgCMatrix") %>%
    as("matrix")
  # Discard any partial bins
  mat <- mat[seq_along(bin_centers), , drop = FALSE]

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
  # Set precision of y-axis range label to within 1% of the max value
  ymax_accuracy <- 10^floor(log10(0.01 * ymax))
  range_label <- sprintf("[0-%s]", scales::label_comma(accuracy = ymax_accuracy, big.mark=" ")(ymax))

  data$normalized_insertions <- pmin(data$normalized_insertions, ymax)

  plot <- ggplot2::ggplot(data) +
    ggplot2::geom_area(ggplot2::aes(pos, pmin(normalized_insertions, ymax), fill = group)) +
    ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
    ggplot2::scale_x_continuous(limits = c(region$start, region$end), expand = c(0, 0), labels=scales::label_number()) +
    ggplot2::scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
    ggplot2::annotate("text", x=region$start, y=ymax, label=range_label, vjust=1.5, hjust=-0.1, size=11*.8/ggplot2::.pt) +
    ggplot2::labs(x = "Genomic Position (bp)", y = "Insertions (RPKM)") +
    ggplot2::guides(y="none", fill="none") +
    ggplot2::facet_wrap("group", ncol=1, strip.position="left") +
    trackplot_theme() 

  wrap_trackplot(plot, ggplot2::unit(ncol(mat), "null"), takes_sideplot=TRUE, region=region)
}

#    GRanges, list, or data.frame of transcript features to plot.
# Required attributes are:
# - `chr`, `start`, `end`: genomic position
# - `strand`: "+"/"-" or TRUE/FALSE for positive or negative strand
# - `feature` (only entries marked as `"transcript"` or `"exon"` will be considered)
# - `transcript_id`
# - `gene_name`

#' Plot transcript models
#' @param transcripts `r document_granges("Transcipt features", strand="default", extras=c("feature"="Only entries marked as \U60\"transcript\"\U60 or \U60\"exon\"\U60 will be considered", "gene_name"="Symbol or gene ID to display", "transcript_id"="Transcritp identifier to link transcripts and exons"))`
#'    
#'    Usually given as the output from `read_gencode_transcripts()`
#' @inheritParams trackplot_coverage
#' @param labels Character vector with labels for each item in transcripts. NA for items that should not be labeled
#' @param exon_size size for exon lines in units of mm
#' @param gene_size size for intron/gene lines in units of mm
#' @param transcript_size size for transcript lines in units of mm
#' @param label_size size for transcript labels in units of mm
#' @return Plot of gene locations
#' @seealso `trackplot_combine()`, `trackplot_coverage()`, `trackplot_loop()`, `trackplot_scalebar()`
#' @examples
#' ## Prep data
#' transcripts <- read_gencode_transcripts(
#'   file.path(tempdir(), "references"), release = "42",
#'   annotation_set = "basic", features = "transcript"
#' )
#' region <- "chr4:3034877-4034877"
#' 
#' 
#' ## Plot gene trackplot
#' plot <- trackplot_gene(transcripts, region)
#' BPCells:::render_plot_from_storage(plot, width = 6, height = 1)
#' @export
trackplot_gene <- function(transcripts, region, exon_size = 2.5, gene_size = 0.5, label_size = 11*.8/ggplot2::.pt, track_label="Genes", return_data = FALSE) {
  region <- normalize_ranges(region)
  transcripts <- normalize_ranges(transcripts, metadata_cols = c("strand", "feature", "transcript_id", "gene_name"))

  # Adjust for the fact that exon_size and gene_size are now in units of linewidth = 0.75mm 
  # whereas size is given in units of 1mm (https://ggplot2.tidyverse.org/articles/ggplot2-specs.html#linewidth)
  size_range <- c(min(exon_size, gene_size, label_size), max(exon_size, gene_size, label_size))
  linewidth_range <- size_range/.75
  exon_size <- exon_size/.75
  gene_size <- gene_size/.75
  data <- transcripts %>%
    tibble::as_tibble() %>%
    dplyr::filter(as.character(chr) == as.character(region$chr), end > region$start, start < region$end, feature %in% c("transcript", "exon")) %>%
    dplyr::mutate(size = dplyr::if_else(feature == "exon", exon_size, gene_size))

  if (nrow(data) == 0) {
    if (return_data) {
      data$y <- numeric(0)
      data$facet_label <- character(0)
      return(list(data=data, arrows=data))
    } else {
      return(trackplot_empty(region, label=track_label))
    }
  }
  # Calculate y positions of each transcript to avoid overlap
  transcript_data <- dplyr::filter(data, feature=="transcript") %>%
    dplyr::mutate(y = trackplot_calculate_segment_height(.))
  if(anyDuplicated(transcript_data$transcript_id)) {
    rlang::abort(sprintf("Found multiple feature == \"transcript\" rows with same ID (e.g. %s)", transcript_data$transcript_id[anyDuplicated(transcript_data$transcript_id)]))
  }
  data <- dplyr::left_join(
    data,
    dplyr::select(transcript_data, transcript_id, y),
    by = "transcript_id"
  )
  if(anyNA(data$y)) {
    rlang::abort(sprintf("Found transcript ID with no feature == \"transcript\" rows: %s", data$transcript_id[which(is.na(data$y))[1]]))
  }
  # Set up data for arrows
  transcript_coords <- dplyr::filter(data, feature == "transcript")
  arrows <- trackplot_create_arrow_segs(transcript_coords, region, 50)
  # Adjust the endpoints of any partially overlapping elements to fit within
  # the plot boundaries
  data <- dplyr::mutate(
    data,
    start = pmax(region$start, pmin(region$end, start)),
    end = pmax(region$start, pmin(region$end, end)),
    facet_label = track_label
  )

  if (return_data) {
    return(list(data=data, arrows=arrows))
  }

  plot <- ggplot2::ggplot(
      data,
      ggplot2::aes(
        x = dplyr::if_else(strand, start, end), xend = dplyr::if_else(strand, end, start),
        y = y, yend = y,
        linewidth = size
      )
    ) +
    # Keep both levels in the color legend even if only one direction is in the viewport by using factors and show.legend=TRUE
    ggplot2::geom_segment(ggplot2::aes(color = factor(strand, levels=c(TRUE,FALSE), labels=c("+","-"))), show.legend=TRUE) +
    ggplot2::geom_segment(ggplot2::aes(color = factor(strand, levels=c(TRUE,FALSE), labels=c("+","-"))), data=arrows, arrow=grid::arrow(length=grid::unit(.4*exon_size, "mm")), show.legend=TRUE) +
    ggrepel::geom_text_repel(
      data = dplyr::filter(data, feature == "transcript"),
      ggplot2::aes(label = gene_name),
      size = label_size,
      position = ggrepel::position_nudge_repel(y = 0.25)
    ) +
    ggplot2::scale_size(range = size_range, limits = size_range) +
    ggplot2::scale_linewidth(range = linewidth_range, limits = linewidth_range) +
    ggplot2::scale_color_manual(values = c("+" = "black", "-" = "darkgrey"), drop=FALSE) +
    ggplot2::scale_x_continuous(limits = c(region$start, region$end), expand = c(0, 0), labels=scales::label_number()) +
    ggplot2::scale_y_discrete(labels = NULL, breaks = NULL) +
    ggplot2::labs(x = "Genomic Position (bp)", y = NULL, color = "strand") +
    ggplot2::guides(size="none", linewidth="none") +
    ggplot2::facet_wrap("facet_label", strip.position="left") +
    trackplot_theme()

  wrap_trackplot(plot, height=ggplot2::unit(1, "null"), region=region)
}

#' Plot range-based annotation tracks (e.g. peaks)
#' @param loci `r document_granges("Genomic loci")`
#' @inheritParams trackplot_coverage
#' @param annotation_size size for annotation lines in mm
#' @param label_by Name of a metadata column in `loci` to use for labeling, or a data vector with same length as loci. Column must hold string data.
#' @param label_size size for labels in units of mm
#' @param color_by Name of a metadata column in `loci` to use for coloring, or a data vector with same length as loci. Column must be numeric or convertible to a factor.
#' @param colors Vector of hex color codes to use for the color scale. For numeric `color_by` data, this is passed to `ggplot2::scale_color_gradientn()`,
#'               otherwise it is interpreted as a discrete color palette in `ggplot2::scale_color_manual()`
#' @param show_strand If TRUE, show strand direction as arrows
#' @return Plot of genomic loci if return_data is FALSE, otherwise returns the data frame used to generate the plot
#' @seealso `trackplot_combine()`, `trackplot_coverage()`, `trackplot_loop()`, `trackplot_scalebar()`, `trackplot_gene()`
#' @examples
#' ## Prep data
#' ## Peaks generated from demo frags, as input into `call_peaks_tile()`
#' peaks <- tibble::tibble(
#'   chr = factor(rep("chr4", 16)),
#'   start = c(3041400, 3041733, 3037400, 3041933, 3040466, 3041200, 
#'             3038200, 3038000, 3040266, 3037733, 3040800, 3042133, 
#'             3038466, 3037200, 3043333, 3040066),
#'   end = c(3041600, 3041933, 3037600, 3042133, 3040666, 3041400, 
#'           3038400, 3038200, 3040466, 3037933, 3041000, 3042333, 
#'           3038666, 3037400, 3043533, 3040266),
#'   enrichment = c(46.4, 43.5, 28.4, 27.3, 17.3, 11.7, 
#'                  10.5, 7.95, 7.22, 6.86, 6.32, 6.14, 
#'                  5.96, 5.06, 4.51, 3.43)
#' )
#' region <- "chr4:3034877-3044877"
#' 
#' ## Plot peaks
#' BPCells:::render_plot_from_storage(
#'   trackplot_genome_annotation(peaks, region, color_by = "enrichment"),
#'   width = 6, height = 1
#' )
#' @export
trackplot_genome_annotation <- function(loci, region, color_by = NULL, colors = NULL, label_by = NULL, label_size = 11*.8/ggplot2::.pt, show_strand=FALSE,
                                        annotation_size = 2.5, track_label="Peaks", return_data = FALSE) {
  region <- normalize_ranges(region)
  assert_true(is.null(label_by) || is.character(label_by))

  data <- trackplot_normalize_ranges_with_metadata(
    loci,
    list("color" = color_by, "label" = label_by, "strand" = if (show_strand) "strand" else NULL)
  )

  if (!is.null(colors) && is.factor(data$color)) {
    if (length(colors) < length(levels(data[["color"]]))) {
      rlang::abort(sprintf("Insufficient values in manual color scale. %d needed but only %d provided.", length(levels(data[["color"]])), length(colors)))
    }
  }

  # Auto-detect color_by label for data vector argument
  if (!is.null(color_by) && !(is.character(color_by) && length(color_by) == 1)) {
    color_by <- argument_name(color_by, 2)
  }

  annotation_size <- annotation_size/.75
  data <- dplyr::filter(data, as.character(chr) == as.character(region$chr), end > region$start, start < region$end)

  # if data is empty, return data or empty plot
  if (nrow(data) == 0) {
    if (return_data) {
      data$y <- numeric(0)
      data$facet_label <- character(0)
      if (show_strand) return(list(data=data, arrows=data))
      else return(list(data=data, arrows=NULL))
    } else {
      return(trackplot_empty(region, track_label))
    }
  }
  data$y <- trackplot_calculate_segment_height(data)
  # Adjust the endpoints of any partially overlapping elements to fit within
  # the plot boundaries
  data <- dplyr::mutate(
    data,
    start = pmax(region$start, pmin(region$end, start)),
    end = pmax(region$start, pmin(region$end, end)),
    facet_label = track_label
  )
  if (show_strand) arrows <- trackplot_create_arrow_segs(data, region, 50, head_only=TRUE)
  else arrows <- NULL
  if (return_data) {
    return(list(data=data, arrows=arrows))
  }

  plot <- ggplot2::ggplot(data, ggplot2::aes(x = start, xend = end, y =y, yend = y)) + 
    ggplot2::geom_segment(linewidth=annotation_size)  + 
    ggplot2::scale_x_continuous(limits = c(region$start, region$end), expand = c(0, 0), labels=scales::label_number()) +
    ggplot2::scale_y_discrete(labels = NULL, breaks = NULL) +
    ggplot2::labs(x = "Genomic Position (bp)", y = NULL) +
    ggplot2::guides(size="none", linewidth="none") +
    ggplot2::facet_wrap("facet_label", strip.position="left") +
    trackplot_theme()
  # add in arrows
  if (show_strand) {
    arrow_angle <- 45
    default_linewidth <- 0.5
    # use trig identity to find arrow length based on annotation size, and account for the linewidth of the ends of arrow
    arrow_length <- ((annotation_size * 0.75/2) + default_linewidth*0.5*0.75*cos(arrow_angle*pi/180))/ sin(arrow_angle*pi/180)
    
    plot <- plot +
      ggplot2::geom_segment(arrows, mapping=ggplot2::aes(
        x = ifelse(strand, start, end),
        xend = ifelse(strand, end, start),
        y = y, yend = y,
      ),
      arrow=grid::arrow(length=grid::unit(arrow_length, "mm"), angle = arrow_angle),
      color="white")
  }

  # add labels if given
  if (!is.null(label_by)) {
    plot <- plot +
      ggrepel::geom_text_repel(
        data = data,
        ggplot2::aes(label = label),
        size = label_size,
        position = ggrepel::position_nudge_repel(y = 0.25),
        color="black"
      )
  }

  # Set up continuous/discrete color scale based on data$color type
  if (!is.null(color_by)) {
    if (is.factor(data$color)) {
      if (is.null(colors)) colors <- discrete_palette("tableau", length(levels(data$color)))
      color_scale <- ggplot2::scale_color_manual(values=colors)
    } else {
      # default continuous palette `c(munsell::mnsl("5P 3/12"), munsell::mnsl("5P 7/12"))`
      if (is.null(colors)) colors <- c("#D099F2", "#6D2884")
      color_scale <- ggplot2::scale_color_gradientn(colors=colors, limits=c(min(data$color), max(data$color)))
    }
    plot <- plot + 
      ggplot2::aes(color = color) +
      color_scale +
      ggplot2::labs(color = color_by)
  }

  wrap_trackplot(plot, height=ggplot2::unit(1, "null"), region=region)
}




#' Plot loops
#'
#' @param loops `r document_granges()`
#' @param color_by Name of a metadata column in `loops` to use for coloring, or a data vector with same length as loci. Column must be numeric or convertible to a factor.
#' @param colors Vector of hex color codes to use for the color scale. For numeric `color_by` data, this is passed to `ggplot2::scale_color_gradientn()`,
#'               otherwise it is interpreted as a discrete color palette in `ggplot2::scale_color_manual()`
#' @param allow_truncated If FALSE, remove any loops that are not fully contained within `region`
#' @param curvature Curvature value between 0 and 1. 1 is a 180-degree arc, and 0 is flat lines.
#' @inheritParams trackplot_coverage
#' 
#' @return Plot of loops connecting genomic coordinates
#' @seealso `trackplot_combine()`, `trackplot_coverage()`, `trackplot_gene()`, `trackplot_scalebar()`, `trackplot_genome_annotation()`
#' @examples
#' peaks <- c(3054877, 3334877, 3534877, 3634877, 3734877)
#' loops <- tibble::tibble(
#'   chr = "chr4",
#'   start = peaks[c(1,1,2,3)],
#'   end = peaks[c(2,3,4,5)],
#'   score = c(4,1,3,2)
#' )
#' region <- "chr4:3034877-4034877"
#' 
#' ## Plot loops
#' plot <- trackplot_loop(loops, region, color_by = "score")
#' BPCells:::render_plot_from_storage(plot, width = 6, height = 1.5)
#' @export
trackplot_loop <- function(loops, region, color_by=NULL, colors=NULL, allow_truncated=TRUE, curvature=0.75, track_label="Links", return_data = FALSE) {
  region <- normalize_ranges(region)
  assert_true(is.null(color_by) || is.numeric(color_by) || is.character(color_by))
  assert_is_numeric(curvature)
  loops <- trackplot_normalize_ranges_with_metadata(loops, metadata=list("color"=color_by))
  
  if (!is.null(colors) && is.factor(loops$color)) {
    if (length(colors) < length(levels(loops[["color"]]))) {
      rlang::abort(sprintf("Insufficient values in manual color scale. %d needed but only %d provided.", length(levels(loops[["color"]])), length(colors)))
    }
  }

  # Auto-detect color_by label for data vector argument
  if (!is.null(color_by) && !(is.character(color_by) && length(color_by) == 1)) {
    color_by <- argument_name(color_by, 2)
  }

  # Calculate curve points (segment of a circle) without scaling
  # Curve is scaled to start at x=0 and end at x=1
  resolution <- 100
  step <- (seq_len(resolution+1)-1)/resolution
  total_angle <- pi*curvature
  min_angle <- 1.5*pi - 0.5*total_angle
  max_angle <- min_angle + total_angle
  curve_x <- (cos(min_angle + step*total_angle) - cos(min_angle))/(cos(max_angle)-cos(min_angle))
  curve_y <- (sin(min_angle + step*total_angle) - sin(min_angle))/(cos(max_angle)-cos(min_angle))

  data <- loops %>%
    tibble::as_tibble() %>%
    dplyr::filter(as.character(chr) == as.character(region$chr),
      ((end > region$start & start < region$end) & !(end > region$end & start < region$start))
    ) %>%
    dplyr::mutate(loop_id = dplyr::row_number()) %>%
    dplyr::cross_join(tibble::tibble(x=curve_x, y=curve_y)) %>%
    dplyr::mutate(
      x = x * (end - start) + start,
      y = y * (end - start),
      facet_label = track_label
    )
  
  if (return_data) {
    return(data)
  }

  if (nrow(data) == 0) {
    return(trackplot_empty(region, track_label))
  }
  
  ymin <- data %>%
    dplyr::filter(end <= region$end, start >= region$start) %>%
    dplyr::summarize(ymin=min(y)) %>%
    dplyr::pull(ymin)

  data <- data %>% 
    dplyr::mutate(
      y = pmax(y, 1.05*ymin),
      x = pmax(region$start, pmin(region$end, x))
    )
  
  plot <- ggplot2::ggplot(data, ggplot2::aes(x, y, group=loop_id)) + 
    ggplot2::geom_line() +
    ggplot2::scale_x_continuous(limits = c(region$start, region$end), expand = c(0, 0), labels=scales::label_number()) +
    ggplot2::scale_y_continuous(labels=NULL, breaks=NULL, expand = c(0.05, 0, 0, 0)) +
    ggplot2::guides(size = "none") +
    ggplot2::labs(x = "Genomic Position (bp)", y = NULL) +
    ggplot2::facet_wrap("facet_label", strip.position="left") + 
    trackplot_theme()

  # Set up continuous/discrete color scale based on data$color type
  if (!is.null(color_by)) {
    if (is.factor(data$color)) {
      if (is.null(colors)) colors <- discrete_palette("tableau", length(levels(data$color)))
      color_scale <- ggplot2::scale_color_manual(values=colors)
    } else {
      # default continuous palette subset from colorbrewer BuPu scheme
      if (is.null(colors)) colors <- c("#bfd3e6","#8c96c6","#88419d","#4d004b")
      color_scale <- ggplot2::scale_color_gradientn(colors=colors, limits=c(min(data$color), max(data$color)))
    }
    plot <- plot + 
      ggplot2::aes(color = color) +
      color_scale +
      ggplot2::labs(color = color_by)
  } 

  wrap_trackplot(plot, ggplot2::unit(1, "null"), region=region)
}

#' Plot scale bar
#'
#' Plots a human-readable scale bar and coordinates of the region being plotted
#' 
#' @param font_pt Font size for scale bar labels in units of pt.
#' @inheritParams trackplot_coverage
#' 
#' @return Plot with coordinates and scalebar for plotted genomic region
#' @seealso `trackplot_combine()`, `trackplot_coverage()`, `trackplot_gene()`, `trackplot_loop()`
#' @examples
#' region <- "chr4:3034877-3044877"
#' BPCells:::render_plot_from_storage(
#'   trackplot_scalebar(region), width = 6, height = 1
#' )
#' @export
trackplot_scalebar <- function(region, font_pt=11) {
  region <- normalize_ranges(region)
  breaks <- pretty(c(region$start, region$end))
  width <- diff(breaks)[1]

  width_text <- scales::label_number(scale_cut = scales::cut_si("b"))(width)
  bar_data <- tibble::tibble(
    right = region$end - 0.05*(region$end-region$start),
    left = right - width,
  )
  scale_label_data <- tibble::tibble(
    right = bar_data$left,
    text = sprintf("%s", width_text)
  )
  number_format <- scales::label_number()
  region_data <- tibble::tibble(
    left = region$start,
    text = sprintf("%s: %s - %s", region$chr, number_format(region$start), number_format(region$end-1L))
  )

  plot <- ggplot2::ggplot() +
    ggplot2::geom_text(data=region_data, ggplot2::aes(x=left, y=0, label=text), size=font_pt/ggplot2::.pt, hjust="left") +
    ggplot2::geom_text(data=scale_label_data, ggplot2::aes(x=right, y=0, label=text), size=font_pt/ggplot2::.pt, hjust=1.1) +
    ggplot2::geom_errorbar(data=bar_data, ggplot2::aes(xmin=left, xmax=right, y=0), width=1) +
    ggplot2::scale_x_continuous(limits = c(region$start, region$end), expand = c(0, 0), labels=scales::label_number()) +
    ggplot2::theme_void() 
  
  wrap_trackplot(plot, height=ggplot2::unit(font_pt*1.1, "pt"), region=region, keep_vertical_margin=TRUE)
}

# This is a still-private all-in-one helper function that is subject to change.
# The arguments will change prior to stabilization and becoming part of the public API
trackplot_helper <- function(gene, clusters, fragments, cell_read_counts, transcripts, loops, loop_color, expression_matrix, flank=1e5) {
  region <- gene_region(transcripts, gene, extend_bp = flank)
  pal <- discrete_palette("stallion")

  base_size <- 11
  scale_plot <- trackplot_scalebar(region, font_pt=base_size)
  bulk_plot <- trackplot_coverage(
    fragments,
    region=region, 
    groups=clusters,
    cell_read_counts=cell_read_counts,
    colors=pal,
    bins=500
  )
  
  gene_plot <- trackplot_gene(transcripts, region) + ggplot2::guides(color="none")

  expression <- collect_features(expression_matrix, gene)
  names(expression) <- "gene"

  expression_plot <- ggplot2::ggplot(expression, ggplot2::aes(clusters, gene, fill=clusters)) +
    ggplot2::geom_boxplot() + 
    ggplot2::guides(y="none", fill="none") + 
    ggplot2::labs(x=NULL, y="RNA") +
    ggplot2::scale_fill_manual(values=pal, drop=FALSE) +
    trackplot_theme()

  loop_plot <- trackplot_loop(loops, region, color_by=loop_color, track_label="Co-Accessibility")

  trackplot_combine(list(scale_plot, bulk_plot, gene_plot, loop_plot), side_plot=expression_plot, title=gene)
}




### Deprecated code:

#' Combine ggplot track plots into an aligned grid.
#' 
#' @description
#' `r lifecycle::badge("deprecated")` 
#' 
#' This function has been renamed to `trackplot_combine()`.
#'
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
#' @keywords internal
draw_trackplot_grid <- function(..., labels, title = NULL,
                                heights = rep(1, length(plots)),
                                label_width = 0.2,
                                label_style = list(fontface = "bold", size = 4)) {
  lifecycle::deprecate_warn("0.2.0", "draw_trackplot_grid()", "trackplot_combine()")
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
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function has been renamed to `trackplot_coverage()`
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
#' @keywords internal
trackplot_bulk <- function(fragments, region, groups,
                           cell_read_counts,
                           group_order = NULL,
                           bins = 200, clip_quantile = 0.999,
                           colors = discrete_palette("stallion"),
                           legend_label = "group",
                           zero_based_coords = !is(region, "GRanges"),
                           return_data = FALSE, return_plot_list = FALSE, apply_styling = TRUE) {
  lifecycle::deprecate_warn("0.2.0", "trackplot_bulk()", "trackplot_coverage()")

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
      ggplot2::geom_area(show.legend = TRUE) +
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
