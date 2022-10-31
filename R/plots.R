

#' Make a knee plot of single cell read counts
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
#' @export
plot_read_count_knee <- function(read_counts, cutoff = NULL, return_data = FALSE, apply_styling = TRUE) {
    # Keep ~1k cells per order of magnitude to speed up plotting
    ranks <- unique(floor(10^(1:(1000 * ceiling(log10(length(read_counts)))) / 1000)))

    data <- list(
        table = tibble::tibble(
            ranks = ranks[ranks < length(read_counts)],
            reads = sort(read_counts, decreasing = TRUE)[ranks]
        ),
        cells_passing = NA
    )

    if (!is.null(cutoff)) data$cells_passing <- sum(read_counts > cutoff)

    if (return_data) {
        return(data)
    }


    plot <- ggplot2::ggplot(data$table, ggplot2::aes(log10(ranks), log10(reads))) +
        ggplot2::geom_point() +
        ggplot2::scale_x_continuous(labels = scales::label_math(), breaks = scales::breaks_width(1)) +
        ggplot2::scale_y_continuous(labels = scales::label_math(), breaks = scales::breaks_width(1))

    subtitle <- NULL
    if (!is.null(cutoff)) {
        plot <- plot +
            ggplot2::geom_vline(xintercept = log10(data$cells_passing), linetype = "dashed")
        subtitle <- sprintf("%s cells with > %s reads", scales::comma(data$cells_passing), scales::comma(floor(cutoff)))
    }

    plot <- plot + ggplot2::labs(x = "Barcode Rank", y = "Reads", subtitle = subtitle)

    if (apply_styling) {
        plot <- plot +
            ggplot2::annotation_logticks() +
            ggplot2::theme_classic()
    }
    return(plot)
}


#' Notes:
#' - scattermore seems like a better choice for rasterized plots
#' Important functionality:
#' - Can pass
#'

# ' Plot a UMAP, t-SNE, or PCA with color highlighting
# '
# '
# '
# ' @param feature Character vector of features to plot if source is not NULL,
# '  or a vector of data to plot if source is NULL.
# ' @param embedding A matrix of dimensions cells x 2 with embedding coordinates
# ' @param source (optional) Matrix, list, or data frame to pull features from. For
# '   a matrix, the features must be rows.
# ' @param ncol Number of columns to use for multi-feature plots
# ' @param quantile_range (optional) Length 2 vector giving the quantiles to clip the minimum and
# '   maximum color scale values, as fractions between 0 and 1. NULL or NA values to skip clipping
# ' @param randomize_order (optional) If TRUE, randomize the order of cells for plotting
# ' @param point_geom ggplot2 point geom for plotting. Can be used to adjust point size,
# '   or use an alternative geom like ggrastr::geom_point_rast
# ' @param
# ' @param match_gene_symbol (optional) A function to match gene symbols/ids from different
# '   systems. Used to match entries in `feature` if source represents a gene matrix.
# '   Function must work like the built-in `match` function, where first argument is the
# '   requested symbols/ids; second argument is the available symbols/ids; return value
# '   is a integer vector giving the correct index for each requested entry in the available
# '   vector (or NA for unmatched IDs)
# ' @inheritParams plot_read_count_knee
# ' @return ggplot2 object
# plot_embedding <- function(feature, embedding, source = NULL, ncol = 3) {}
#                            quantile_range = c(0.01, 0.99),
#                            randomize_order = TRUE
#                            match_gene_symbol = match_gene_symbol_human,
#                            return_data = FALSE, apply_styling = TRUE) {


# }