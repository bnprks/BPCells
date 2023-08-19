
#' Download a file with a custom timeout
#'
#' @param path Output path to write file
#' @param url to download from
#' @param timeout timeout in seconds
#' @keywords internal
ensure_downloaded <- function(path, backup_url, timeout) {
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout))
  options(timeout = timeout)
  if (!is.null(backup_url) && !file.exists(path)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    download.file(backup_url, path, mode="wb")
  }
}

#' Read GTF gene annotations
#'
#' Read gene annotations from gtf format into a data frame. The source can be
#' a URL, a gtf file on disk, or a gencode release version.
#'
#' @rdname read_gtf
#' @details **read_gtf**
#'
#' Read gtf from a file or URL
#' @param path Path to file (or desired save location if backup_url is used)
#' @param attributes Vector of GTF attribute names to parse out as columns
#' @param tags Vector of tags to parse out as boolean presence/absence
#' @param features List of features types to keep from the GTF (e.g. gene, transcript, exon, intron)
#' @param keep_attribute_column Boolean for whether to preserve the raw attribute text column
#' @param backup_url If path does not exist, provides a URL to download the gtf from
#' @param timeout Maximum time in seconds to wait for download from backup_url
#' @return Data frame with coordinates using the 0-based convention. Columns are:
#'
#' - chr
#' - source
#' - feature
#' - start
#' - end
#' - score
#' - strand
#' - frame
#' - attributes (optional; named according to listed attributes)
#' - tags (named according to listed tags)
#' @seealso [read_bed()], [read_encode_blacklist()]
#' @export
read_gtf <- function(path, attributes = c("gene_id"), tags = character(0), features = c("gene"), keep_attribute_column = FALSE, backup_url = NULL, timeout = 300) {
  assert_has_package("readr")
  assert_is_file(path, must_exist = is.null(backup_url), multiple_ok = FALSE)
  ensure_downloaded(path, backup_url, timeout)
  gtf_colnames <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
  gtf <- readr::read_tsv(path, comment = "#", col_names = gtf_colnames, col_types = "ccciicccc")
  gtf <- gtf[gtf$feature %in% features, ]
  gtf$start <- gtf$start - 1
  for (a in attributes) {
    gtf[[a]] <- stringr::str_match(gtf[["attributes"]], sprintf('%s "([^"]*)"', a))[, 2]
  }
  for (t in tags) {
    gtf[[t]] <- stringr::str_detect(gtf[["attributes"]], stringr::fixed(sprintf('tag "%s"', t)))
  }
  if (!keep_attribute_column) {
    gtf[["attributes"]] <- NULL
  }
  return(gtf)
}

#' @rdname read_gtf
#' @details **read_gencode_genes**
#'
#' Read gene annotations directly from GENCODE. The file name will vary depending
#' on the release and annotation set requested, but will be of the format
#' `gencode.v42.annotation.gtf.gz`. GENCODE currently recommends the basic set:
#' <https://www.gencodegenes.org/human/>. In release 42, both the comprehensive and
#' basic sets had identical gene-level annotations, but the comprehensive set had
#' additional transcript variants annotated.
#'
#' @param dir Output directory to cache the downloaded gtf file
#' @param release release version (prefix with M for mouse versions). For
#'   most recent version, use "latest" or "latest_mouse"
#' @param annotation_set Either "basic" or "comprehensive" annotation sets (see details section).
#' @param gene_type Regular expression with which gene types to keep.
#' Defaults to protein_coding, lncRNA, and IG/TR genes
#' @export
read_gencode_genes <- function(dir, release = "latest",
                               annotation_set = c("basic", "comprehensive"),
                               gene_type = "lncRNA|protein_coding|IG_.*_gene|TR_.*_gene",
                               attributes = c("gene_id", "gene_type", "gene_name"), tags = character(0),
                               features = c("gene"), timeout = 300) {
  assert_true("gene_type" %in% attributes)
  annotation_set <- match.arg(annotation_set)
  # Detect latest GENCODE release if needed
  if (release %in% c("latest", "latest_mouse")) {
    species <- if (release == "latest") "human" else "mouse"
    md5_url <- sprintf("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_%s/latest_release/MD5SUMS", species)
    conn <- url(md5_url)
    lines <- readLines(conn)
    close(conn)
    matches <- stringr::str_match(lines, "gencode\\.v([^.]+)\\.annotation.gtf.gz")[, 2] #nolint
    matches <- matches[!is.na(matches)]
    if (length(matches) != 1) {
      stop("Error fetching latest release number from gencode. Please specify explicit release")
    }
    release <- matches[1]
  }

  species <- if (stringr::str_detect(release, "^M")) "mouse" else "human"
  set <- if (annotation_set == "basic") "basic." else ""
  url <- sprintf(
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_%s/release_%s/gencode.v%2$s.%sannotation.gtf.gz",
    species, release, set
  )
  path <- file.path(dir, basename(url))
  gtf <- read_gtf(path, attributes = attributes, tags = tags, features = features, backup_url = url, timeout = timeout)
  gtf <- gtf[stringr::str_detect(gtf$gene_type, paste0("^(", gene_type, ")$")), ]
  return(gtf)
}

#' @rdname read_gtf
#' @param transcript_choice Method for selecting representative transcripts. Choices are:
#'    - MANE_Select: human-only, most conservative
#'    - Ensembl_Canonical: human+mouse, superset of MANE_Select for human
#'    - all: Preserve all transcript models (not recommended for plotting)
#' @details **read_gencode_transcripts**
#'
#' Read transcript models from GENCODE, for use with trackplot_gene()
#' @export 
read_gencode_transcripts <- function(dir, release = "latest", transcript_choice = c("MANE_Select", "Ensembl_Canonical", "all"),
                                     annotation_set = c("basic", "comprehensive"),
                                     gene_type = "lncRNA|protein_coding|IG_.*_gene|TR_.*_gene",
                                     attributes = c("gene_id", "gene_type", "gene_name", "transcript_id"),
                                     features = c("transcript", "exon"),
                                     timeout = 300) {
  transcript_choice <- match.arg(transcript_choice)
  species <- if (stringr::str_detect(release, "(^M)|(latest_mouse)")) "mouse" else "human"
  if (species == "mouse" && transcript_choice == "MANE_Select") {
    rlang::warn("MANE_Select transcripts not available for mouse. Defaulting to Ensembl_Canonical")
    transcript_choice <- "Ensembl_Canonical"
  }
  tags <- if (transcript_choice != "all") transcript_choice else character(0)
  gtf <- read_gencode_genes(
    dir,
    release = release,
    annotation_set = annotation_set,
    gene_type = gene_type,
    attributes = attributes,
    tags = tags,
    features = features,
    timeout = timeout
  )
  if (transcript_choice == "MANE_Select") {
    return(gtf[gtf$MANE_Select, ])
  } else if (transcript_choice == "Ensembl_Canonical") {
    return(gtf[gtf$Ensembl_Canonical, ])
  } else {
    return(gtf)
  }
}




#' Read a bed file into a data frame
#'
#' Bed files can contain peak or blacklist annotations. These utilities help
#' read thos annotations
#' @rdname read_bed
#' @details **read_bed**
#'
#' Read a bed file from disk or a url.
#' @inheritParams read_gtf
#' @param additional_columns Names for additional columns in the bed file
#' @return Data frame with coordinates using the 0-based convention.
#' @seealso [read_gtf()], [read_gencode_genes()]
#' @export
read_bed <- function(path, additional_columns = character(0), backup_url = NULL, timeout = 300) {
  assert_has_package("readr")
  assert_is_file(path, must_exist = is.null(backup_url), multiple_ok = FALSE)
  ensure_downloaded(path, backup_url, timeout)
  bed <- readr::read_tsv(path, col_names = c("chr", "start", "end", additional_columns), col_types = "") %>%
    dplyr::select(c(chr, start, end, additional_columns))
  return(bed)
}


#' @rdname read_bed
#' @details **read_encode_blacklist**
#'
#' Downloads the Boyle Lab blacklist, as described in <https://doi.org/10.1038/s41598-019-45839-z>
#' @param genome genome name
#' @export
read_encode_blacklist <- function(dir, genome = c("hg38", "mm10", "hg19", "dm6", "dm3", "ce11", "ce10"), timeout = 300) {
  genome <- match.arg(genome)
  url <- sprintf("https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/%s-blacklist.v2.bed.gz", genome)
  path <- file.path(dir, basename(url))
  return(read_bed(path, backup_url = url, additional_columns = "reason"))
}

#' Read UCSC chromosome sizes
#'
#' Read chromosome sizes from UCSC and return as a tibble with one row per
#' chromosome.
#' The underlying data is pulled from here: <https://hgdownload.soe.ucsc.edu/downloads.html>
#'
#' @export
read_ucsc_chrom_sizes <- function(dir, genome = c("hg38", "mm39", "mm10", "mm9", "hg19"),
                                  keep_chromosomes = "chr[0-9]+|chrX|chrY", timeout = 300) {
  genome <- match.arg(genome)
  url <- sprintf("https://hgdownload.soe.ucsc.edu/goldenPath/%s/bigZips/%1$s.chrom.sizes", genome)
  path <- file.path(dir, basename(url))
  ensure_downloaded(path, url, timeout)
  sizes <- readr::read_tsv(path, col_names = c("chr", "end"), col_types = "ci") %>%
    dplyr::transmute(chr, start = 0, end)
  if (!is.null(keep_chromosomes)) {
    sizes <- sizes %>%
      dplyr::filter(stringr::str_detect(chr, paste0("^(", keep_chromosomes, ")$")))
  }
  return(sizes)
}


#' Get end-sorted ordering for genome ranges
#'
#' Use this function to order regioins prior to calling `peak_matrix()` or `tile_matrix()`.
#'
#' @inheritParams normalize_ranges
#' @param chr_levels Ordering of chromosome names
#' @param sort_by_end If TRUE (defualt), sort by (chr, end, start). Else sort by (chr, start, end)
#' @return Numeric vector analagous to the `order` function. Provides an index
#' selection that will reorder the input ranges to be sorted by chr, end, start
#' @export
order_ranges <- function(ranges, chr_levels, sort_by_end = TRUE) {
  ranges <- normalize_ranges(ranges)
  if (sort_by_end) {
    return(order(match(as.character(ranges$chr), chr_levels), ranges$end, ranges$start))
  } else {
    return(order(match(as.character(ranges$chr), chr_levels), ranges$start, ranges$end))
  }
}

# Helper function to remove ensembl version suffixes from gene names
remove_ensembl_version <- function(vec) {
  is_versioned_ensembl <- stringr::str_detect(vec, "^ENS(MUS)?G[0-9]{11}\\.[0-9]+")
  vec[is_versioned_ensembl] <- stringr::str_remove(vec[is_versioned_ensembl], "\\.[0-9]+")
  vec
}

#' Gene symbol matching
#'
#' Correct alias gene symbols, Ensembl IDs, and Entrez IDs to canonical gene
#' symbols. This is useful for matching gene names between different datasets
#' which might not always use the same gene naming conventions.
#'
#' @rdname gene_matching
#' @param query Character vector of gene symbols or IDs
#' @param subject Vector of gene symbols or IDs to index into
#' @param gene_mapping Named vector where names are gene symbols or IDs and values are
#'   canonical gene symbols
#' @return **match_gene_symbol**
#'
#' Integer vector of indices `v` such that `subject[v]` corresponds to the gene
#' symbols in `query`
#' @export
match_gene_symbol <- function(query, subject, gene_mapping = human_gene_mapping) {
  assert_is(query, "character")
  assert_is(subject, "character")
  assert_is(gene_mapping, "character")
  # Any exact match counts
  res <- match(query, subject, incomparables = NA)

  initial_missing <- sum(is.na(res))
  orig_query <- query
  query <- human_gene_mapping[remove_ensembl_version(query)]
  subject <- human_gene_mapping[remove_ensembl_version(subject)]

  corrected_match <- match(query, subject, incomparables = NA)

  res[is.na(res)] <- corrected_match[is.na(res)]

  corrected_missing <- sum(is.na(res))
  if (corrected_missing > 0) {
    rlang::warn(sprintf(
      "Could not match %d symbols: %s",
      corrected_missing,
      pretty_print_vector(orig_query[is.na(res)])
    ))
  }
  return(res)
}

#' @rdname gene_matching
#' @return **canonical_gene_symbol**
#'
#' Character vector of canonical gene symbols for each symbol in `query`
#' @export
canonical_gene_symbol <- function(query, gene_mapping = human_gene_mapping) {
  assert_is(query, "character")
  assert_is(gene_mapping, "character")
  res <- remove_ensembl_version(query)
  res <- human_gene_mapping[res]
  names(res) <- query
  if (anyNA(res)) {
    rlang::warn(sprintf(
      "Could not match %d symbols: %s",
      sum(is.na(res)),
      pretty_print_vector(query[is.na(res)])
    ))
  }
  return(res)
}

#' Find gene region
#'
#' Conveniently look up the region of a gene by gene symbol. The value returned by this function
#' can be used as the `region` argument for trackplot functions such as
#' `trackplot_bulk()` or `trackplot_gene()`
#'
#' @inheritParams match_gene_symbol
#' @param genes GRanges, list, or data.frame of transcript features to plot.
#' Required attributes are:
#' - `chr`, `start`, `end`: genomic position
#' - `gene_name`: Symbol or gene ID
#' @param gene_symbol Name of gene symbol or ID
#' @param extend_bp Bases to extend region upstream and downstream of gene
#' @return list of chr, start, end positions for use with , etc.
#' @export
gene_region <- function(genes, gene_symbol, extend_bp = 1e4, gene_mapping = human_gene_mapping) {
  genes <- normalize_ranges(genes, metadata_cols = c("gene_name"))
  idx <- match_gene_symbol(gene_symbol, genes$gene_name)
  if (is.na(idx)) {
    rlang::stop("Could not locate gene")
  }
  return(list(
    chr = as.character(genes$chr[idx]),
    start = genes$start[idx] - extend_bp,
    end = genes$end[idx] + extend_bp
  ))
}
