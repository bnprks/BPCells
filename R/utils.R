
#' Download a file with a custom timeout
#' 
#' @param path Output path to write file
#' @param url to download from
#' @param timeout timeout in seconds
ensure_downloaded <- function(path, backup_url, timeout) {
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout))
  options(timeout = timeout)
  if (!is.null(backup_url) && !file.exists(path)) {
    dir.create(dirname(path), recursive=TRUE, showWarnings = FALSE)
    download.file(backup_url, path)
  }
}

#' Read a gtf file into a data frame
#' @param path Path to file (or desired save location if backup_url is used)
#' @param attributes Vector of GTF attribute names to parse out as columns
#' @param features List of features types to keep from the GTF (e.g. gene, transcript, exon, intron)
#' @param keep_attribute_column Boolean for whether to preserve the raw attribute text column
#' @param backup_url If path does not exist, provides a URL to download the gtf from
#' @param timeout Maximum time in seconds to wait for download from backup_url
#' @return data frame with coordinates using the 0-based convention
#' @export
read_gtf <- function(path, attributes=c("gene_id"), features=c("gene"), keep_attribute_column=FALSE, backup_url=NULL, timeout=300) {
  assert_has_package("readr")
  assert_is_file(path, must_exist = is.null(backup_url), multiple_ok = FALSE)
  ensure_downloaded(path, backup_url, timeout)
  gtf_colnames <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
  gtf <- readr::read_tsv(path, comment = "#", col_names = gtf_colnames, col_types = "ccciicccc")
  gtf <- gtf[gtf$feature %in% features,]
  gtf$start <- gtf$start - 1
  for (a in attributes) {
    gtf[[a]] <- stringr::str_match(gtf[["attributes"]], sprintf('%s "([^"]*)"', a))[,2]
  }
  if(!keep_attribute_column)
    gtf[["attributes"]] <- NULL
  return(gtf)
}

#' Read a bed file file into a data frame
#' @inheritParams read_gtf
#' @param additional_columns Names for additional columns in the bed file
#' @return data frame with coordinates using the 0-based convention
#' @export
read_bed <- function(path, additional_columns=character(0), backup_url=NULL, timeout=300) {
  assert_has_package("readr")
  assert_is_file(path, must_exist = is.null(backup_url), multiple_ok = FALSE)
  ensure_downloaded(path, backup_url, timeout)
  bed <- readr::read_tsv(path, col_names=c("chr", "start", "end", additional_columns), col_types="") %>%
    dplyr::select(c(col_names, additional_columns))
  return(bed)
}

#' Download gene annotations from gencode
#' @param dir output directory to cache the downloaded file
#' @param release release version (prefix with M for mouse versions). For
#'   most recent version, use "latest" or "latest_mouse"
#' @param gene_type Regular expression with which gene types to keep.  
#' Defaults to protein_coding, lncRNA, and IG/TR genes
#' @inheritParams read_gtf
#' @return Data frame with coordinates using the 0-based convention
#' @export
read_gencode_genes <- function(dir, release="41", gene_type="lncRNA|protein_coding|IG_.*_gene|TR_.*_gene", timeout=300) {
  # Detect latest GENCODE release if needed
  if (release %in% c("latest", "latest_mouse")) {
    species <- if (release == "latest") "human" else "mouse"
    md5_url <- sprintf("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_%s/latest_release/MD5SUMS", species)
    conn <- url(md5_url)
    lines <- readLines(conn)
    close(conn)
    matches <- stringr::str_match(lines, "gencode\\.v([^.]+)\\.annotation.gtf.gz")[,2]
    matches <- matches[!is.na(matches)]
    if (length(matches) != 1) {
      stop("Error fetching latest release number from gencode. Please specify explicit release")
    }
    release <- matches[1]
  }

  species <- if (stringr::str_detect(release, "^M")) "mouse" else "human"
  url <- sprintf(
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_%s/release_%s/gencode.v%2$s.annotation.gtf.gz",
     species, release
  )
  path <- file.path(dir, basename(url))
  gtf <- read_gtf(path, attributes=c("gene_id", "gene_type", "gene_name"), backup_url=url)
  gtf <- gtf[stringr::str_detect(gtf$gene_type, gene_type),]
  return(gtf)
}

#' Download ENCODE blacklist annotations
#' 
#' Downloads the Boyle Lab blacklist, as described in https://doi.org/10.1038/s41598-019-45839-z
#' @inheritParams read_gencode_genes
#' @param genome genome name
#' @return Data frame with coordinates using the 0-based convention
#' @export
read_encode_blacklist <- function(dir, genome=c("hg38", "mm10", "hg19", "dm6", "dm3", "ce11", "ce10"), timeout=300) {
  genome <- match.arg(genome)
  url <- sprintf("https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/%s-blacklist.v2.bed.gz", genome)
  path <- file.path(dir, basename(url))
  return(read_bed(path, backup_url=url, additional_columns = "reason"))
}


#' Get end-sorted ordering for genome ranges
#' 
#' Use this function to order regioins prior to calling peakMatrix or tileMatrix.
#' 
#' @inheritParams normalize_ranges
#' @param chr_levels Ordering of chromosome names
#' @return Numeric vector analagous to the `order` function. Provides an index
#' selection that will reorder the input ranges to be sorted by chr, end, start
#' @export
order_ranges <- function(ranges, chr_levels) {
  ranges <- normalize_ranges(ranges)
  order(match(as.character(ranges$chr), chr_levels), ranges$end, ranges$start)
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
  # Any exact match counts
  res <- match(query, subject, incomparables = NA)

  initial_missing <- sum(is.na(res))

  query <- human_gene_mapping[remove_ensembl_version(query)]
  subject <- human_gene_mapping[remove_ensembl_version(subject)]

  corrected_match <- match(query, subject, incomparables = NA)

  res[is.na(res)] <- corrected_match[is.na(res)]

  corrected_missing <- sum(is.na(res))
  if (corrected_missing > 0) {
    rlang::warn(sprintf(
      "Could not match %d symbols: %s",
      corrected_missing,
      pretty_print_vector(query[is.na(res)])
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
  res <- remove_ensembl_version(query)
  res <- human_gene_mapping[res]
  names(res) <- query
  if (any(is.na(res))) {
    rlang::warn(sprintf(
      "Could not match %d symbols: %s",
      sum(is.na(res)),
      pretty_print_vector(query[is.na(res)])
    ))
  }
  return(res)
}