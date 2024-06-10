# Copyright 2021 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

# Error checking utilities
# Design: Provide assert_* functions that will raise an error in the event of issues,
#   with enough context to help the user find the source of the problem.
#   normalize_* functions return a normalized verison of the first argument, while
#   raising an error via assert if the argument cannot be normalized.
#   To avoid having the assert functions themselves in the call stack, all of
#   them take an argument "n", which specifies how many calls on the stack to
#   skip over during error printing.

argument_name <- function(arg, n) {
  arg_name <- substitute(arg)
  for (f in seq_len(n)) {
    arg_name <- do.call(substitute, list(arg_name, parent.frame(f)))
  }
  arg_name <- deparse(arg_name)
  return(arg_name)
}

# Pretty print an error on behalf of a calling function
# Goes up n levels on the call stack, and prints the offending
# function name, argument name, and stack trace along with the provided
# error message.
# Calls stop, causing an error.
# @param arg The argument which caused the error
# @param msg The message to print
# @param n Number of call levels to go up
# @keywords internal
pretty_error <- function(arg, msg, n) {
  arg_name <- argument_name(arg, n + 1)

  fn_name <- rlang::call_name(sys.call(-(1 + n)))

  trace <- paste0(format(rlang::trace_back(bottom = parent.frame(n))), collapse = "\n")

  message <- sprintf(
    "In function %s: \"%s\" %s\n%s",
    fn_name, arg_name, msg, trace
  )

  stop(message, call. = FALSE)
}

assert_is_wholenumber <- function(x, n = 1) {
  assert_is_numeric(x, n + 1)
  # Taken from documentation for is.integer
  if (any(abs(x - round(x)) >= .Machine$double.eps^0.5)) {
    pretty_error(x, "must be a whole number", n)
  }
}

assert_is_numeric <- function(x, n = 1) {
  if (!is.numeric(x)) pretty_error(x, "must be numeric", n)
}

assert_is_character <- function(x, n = 1) {
  if (!is.character(x)) pretty_error(x, "must be a character vector", n)
}

assert_len <- function(x, len, n = 1) {
  if (length(x) != len) pretty_error(x, sprintf("must have length %s", len), n)
}

# Assertions for file paths, with adjustable criteria.
# Always requires that the path is a character vector
# @param path File path to test
# @param must_exist Whether the path must exist on disk
# @param multiple_ok Whether the length of path can be > 1
# @param extension If not null, list of valid extensions (including the ".", e.g. ".tsv")
# @param n How many caller frames above to use as the environment for error printing
# @keywords internal
assert_is_file <- function(path, must_exist = TRUE, multiple_ok = FALSE, extension = NULL, n = 1) {
  assert_is_character(path, n + 1)
  if (!multiple_ok) assert_len(path, 1, n + 1)

  if (!is.null(extension)) {
    for (p in path) {
      if (!any(endsWith(p, extension))) {
        pretty_error(
          path,
          sprintf("must have file extension %s", paste0(extension, collapse = ", or ")),
          n
        )
      }
    }
  }

  if (must_exist && !all(file.exists(path))) pretty_error(path, "is not a valid file path (file not found)", n)
}

assert_distinct <- function(vector, n = 1) {
  if (anyDuplicated(vector)) pretty_error(vector, "must not have duplicated values", n)
}
assert_greater_than_zero <- function(vector, n = 1) {
  assert_is_numeric(vector, n + 1)
  if (any(vector <= 0)) pretty_error(vector, "must be greater than zero", n)
}
assert_is <- function(object, class, n = 1) {
  if (length(class) == 1) {
    if (!is(object, class)) pretty_error(object, sprintf("must have class \"%s\"", class), n)
  } else {
    match <- FALSE
    for (c in class) {
      match <- is(object, c) || match
    }
    if (!match) pretty_error(object, sprintf("must have class %s", paste0(class, collapse = ", or ")), n)
  }
}

assert_true <- function(expr, n = 1) {
  if (!expr) pretty_error(expr, "is not true", n)
}

assert_has_names <- function(x, names, n = 1) {
  if (!all(names %in% names(x))) {
    pretty_error(x, sprintf("is missing names: %s", toString(setdiff(names, names(x)))), n)
  }
}

assert_not_na <- function(x, n = 1) {
  if (anyNA(x)) pretty_error(x, "has NA entries", n)
}

assert_not_null <- function(x, n = 1) {
  if (is.null(x)) pretty_error(x, "is NULL", n)
}


normalize_length <- function(x, len, n = 1) {
  if (length(x) == 1) {
    return(rep(x, len))
  }
  assert_len(x, len, n + 1)
  return(x)
}

assert_has_package <- function(packages, n = 1) {
  missing <- c()
  for (p in packages) {
    if (!requireNamespace(p, quietly = TRUE)) {
      missing <- c(missing, p)
    }
  }
  if (length(missing) > 0) {
    stop(sprintf("Missing required package(s): %s", toString(missing)), .call = FALSE)
  }
}


#' Normalize an object representing genomic ranges
#'
#' @param ranges `r document_granges()`
#' @param metadata_cols Optional list of metadata columns to require & extract
#' @param zero_based_coords If true, coordinates start and 0 and the end coordinate is not included in the range.
#'  If false, coordinates start at 1 and the end coordinate is included in the range
#' @return data frame with zero-based coordinates, and elements chr (factor), start (int), and end (int).
#' If `ranges` does not have chr level information, chr levels are the sorted unique values of chr.
#'
#' If strand is in metadata_cols, then the output strand element will be TRUE for positive strand,
#' and FALSE for negative strand. (Converted from a character vector of "+"/"-" if necessary)
normalize_ranges <- function(ranges, metadata_cols = character(0), zero_based_coords = !is(ranges, "GRanges"), n = 1) {
  assert_is(ranges, c("GRanges", "list", "data.frame", "character"), n = n + 1)
  assert_is(metadata_cols, "character", n = n + 1)
  assert_is(zero_based_coords, "logical", n = n + 1)

  ret <- list()
  if (is(ranges, "character")) {
    assert_len(metadata_cols, 0, n = n + 1)
    ret_parse <- stringr::str_match(ranges, "^([^:]+):([0-9,]+)-([0-9,]+)$")
    if (anyNA(ret_parse)) {
      rlang::abort("Could not parse ranges from character. Must have form chr1:2501-2900")
    }
    chr_levels <- sort(unique(ret_parse[, 2]))
    ret$chr <- factor(ret_parse[, 2], levels = chr_levels)
    ret$start <- as.integer(stringr::str_remove_all(ret_parse[, 3], ",")) - !zero_based_coords
    ret$end <- as.integer(stringr::str_remove_all(ret_parse[, 4], ","))
  } else if (is(ranges, "GRanges")) {
    if ("strand" %in% metadata_cols) ret$strand <- as.logical(GenomicRanges::strand(ranges) != "-")
    metadata_cols <- metadata_cols[metadata_cols != "strand"]

    assert_has_names(GenomicRanges::mcols(ranges), metadata_cols, n = n + 1)
    chr_levels <- GenomeInfoDb::seqlevels(ranges)
    ret$chr <- factor(as.character(GenomicRanges::seqnames(ranges)), levels = chr_levels)
    ret$start <- GenomicRanges::start(ranges) - !zero_based_coords
    ret$end <- GenomicRanges::end(ranges)

    for (c in metadata_cols) {
      ret[[c]] <- GenomicRanges::mcols(ranges)[[c]]
    }
  } else {
    assert_has_names(ranges, c("chr", "start", "end", metadata_cols), n = n + 1)
    chr_levels <- levels(ranges$chr)
    if (is.null(chr_levels)) chr_levels <- sort(unique(as.character(ranges$chr)))
    ret$chr <- factor(as.character(ranges$chr), levels = chr_levels)
    ret$start <- as.integer(ranges$start) - !zero_based_coords
    ret$end <- as.integer(ranges$end)
    if ("strand" %in% metadata_cols) {
      assert_is(ranges$strand, c("character", "logical"), n = n + 1)
      if (is(ranges$strand, "character")) {
        ret$strand <- ranges$strand != "-"
      } else {
        ret$strand <- ranges$strand
      } # If strand is already logical
    }
    metadata_cols <- metadata_cols[metadata_cols != "strand"]
    for (c in metadata_cols) {
      ret[[c]] <- ranges[[c]]
    }
  }
  return(tibble::as_tibble(ret))
}
