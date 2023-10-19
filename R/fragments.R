#' IterableFragments methods
#'
#' Methods for IterableFragments objects
#'
#' @name IterableFragments-methods
#' @rdname IterableFragments-methods
NULL

setClass("IterableFragments")


# Get an external pointer to a C++ FragmentsLoader object. This is only called
# when an operation is actively ready to run
setGeneric("iterate_fragments", function(x) standardGeneric("iterate_fragments"))

# Provide a short description of the IterableFragments, used for
# pretty-printing pipelines.
# Returns a vector of strings, each of which describes one step of the
# transformation pipeline
setGeneric("short_description", function(x) standardGeneric("short_description"))
setMethod("short_description", "IterableFragments", function(x) {
  character(0)
})

#' @describeIn IterableFragments-methods Print IterableFragments
#' @param object IterableFragments object
setMethod("show", "IterableFragments", function(object) {
  cat(sprintf("IterableFragments object of class \"%s\"\n", class(object)))

  cat("\n")

  if (is.null(cellNames(object))) {
    cat("Cells: count unknown\n")
  } else {
    cat(sprintf(
      "Cells: %d cells with %s\n", length(cellNames(object)),
      pretty_print_vector(cellNames(object), prefix = "names ", empty = "unknown names")
    ))
  }

  if (is.null(chrNames(object))) {
    cat("Chromosomes: count unknown\n")
  } else {
    cat(sprintf(
      "Chromosomes: %d chromosomes with %s\n", length(chrNames(object)),
      pretty_print_vector(chrNames(object), prefix = "names ", empty = "unknown names")
    ))
  }

  cat("\n")
  description <- short_description(object)
  if (length(description) > 0) cat("Queued Operations:\n")
  for (i in seq_along(description)) {
    cat(sprintf("%d. %s\n", i, description[i]))
  }
})

# Generic setters/getters for cell and chromosome names

#' Get cell names
#' @param x an IterableFragments object
#' @return * `cellNames()` Character vector of cell names, or NULL if none are known
#' @describeIn IterableFragments-methods Get cell names
#' @export
setGeneric("cellNames", function(x) standardGeneric("cellNames"))
setMethod("cellNames", "IterableFragments", function(x) {
  if (.hasSlot(x, "fragments")) {
    return(cellNames(x@fragments))
  }
  stop(sprintf("Error: cellNames not implemented on object of class %s", class(x)))
})

#' Set cell names
#' @param x an IterableFragments object
#' @param value Character vector of new names
#' @details * `cellNames<-` It is only possible to replace names, not add new names.
#' @describeIn IterableFragments-methods Set cell names
#' @export
setGeneric("cellNames<-", function(x, ..., value) standardGeneric("cellNames<-"))
setMethod("cellNames<-", "IterableFragments", function(x, ..., value) {
  if (is.null(cellNames(x))) {
    stop("Assigning new cellNames is not allowed, only renaming")
  }
  assert_is_character(value)
  assert_len(value, length(cellNames(x)))
  new("CellRename", fragments = x, cell_names = value)
})

#' Get chromosome names
#' @param x an IterableFragments object
#' @return * `chrNames()`: Character vector of chromosome names, or NULL if none are known
#' @describeIn IterableFragments-methods Set chromosome names
#' @export
setGeneric("chrNames", function(x) standardGeneric("chrNames"))
setMethod("chrNames", "IterableFragments", function(x) {
  if (.hasSlot(x, "fragments")) {
    return(chrNames(x@fragments))
  }
  stop(sprintf("Error: chrNames not implemented on object of class %s", class(x)))
})

#' Set chromosome names
#' @param x an IterableFragments object
#' @param value Character vector of new names
#' @details * `chrNames<-` It is only possible to replace names, not add new names.
#' @describeIn IterableFragments-methods Set chromosome names
#' @export
setGeneric("chrNames<-", function(x, ..., value) standardGeneric("chrNames<-"))
setMethod("chrNames<-", "IterableFragments", function(x, ..., value) {
  if (is.null(chrNames(x))) {
    stop("Assigning new chrNames is not allowed, only renaming")
  }
  assert_is_character(value)
  assert_len(value, length(chrNames(x)))
  new("ChrRename", fragments = x, chr_names = value)
})

# Read/write from 10x fragment files
setClass("FragmentsTsv",
  contains = "IterableFragments",
  slots = c(
    path = "character",
    comment = "character"
  ),
  prototype = list(
    path = NA_character_,
    comment = ""
  )
)
setMethod("chrNames", "FragmentsTsv", function(x) NULL)
setMethod("cellNames", "FragmentsTsv", function(x) NULL)

setMethod("iterate_fragments", "FragmentsTsv", function(x) iterate_10x_fragments_cpp(normalizePath(x@path), x@comment))
setMethod("short_description", "FragmentsTsv", function(x) {
  sprintf("Load 10x fragments file from %s", x@path)
})

#' Read/write a 10x fragments file
#'
#' 10x fragment files come in a bed-like format, with columns chr, start, end,
#' cell_id, and pcr_duplicates. Unlike a standard bed format, the format from
#' cellranger has an *inclusive* end-coordinate, meaning the end coordinate itself
#' is what should be counted as the tagmentation site, rather than offset by 1.
#'
#' @details **open_fragments_10x**
#'
#' No disk operations will take place until the fragments are used in a function
#' @param path File path (e.g. fragments.tsv or fragments.tsv.gz)
#' @param comment Skip lines at beginning of file which start with comment string
#' @param end_inclusive Whether the end coordinate of the bed is inclusive -- i.e. there was an
#'     insertion at the end coordinate rather than the base before the end coordinate. This is the
#'     10x default, though it's not quite standard for the bed file format.
#' @return 10x fragments file object
#' @export
open_fragments_10x <- function(path, comment = "#", end_inclusive = TRUE) {
  assert_is_file(path, extension = c(".tsv", ".tsv.gz"))
  assert_is_character(comment)
  assert_len(comment, 1)
  path <- normalizePath(path)
  res <- new("FragmentsTsv", path = path, comment = comment)
  if (end_inclusive) {
    res <- shift_fragments(res, shift_end = 1)
  }
  res
}


#' @rdname open_fragments_10x
#' @param fragments Input fragments object
#' @param append_5th_column Whether to include 5th column of all 0 for compatibility
#'        with 10x fragment file outputs (defaults to 4 columns chr,start,end,cell)
#' @details **write_fragments_10x**
#'
#' Fragments will be written to disk immediately, then returned in a readable object.
#' @export
write_fragments_10x <- function(fragments, path, end_inclusive = TRUE, append_5th_column = FALSE) {
  assert_is_file(path, must_exist = FALSE, extension = c(".tsv", ".tsv.gz"))
  if (end_inclusive) {
    fragments <- shift_fragments(fragments, shift_end = -1)
  }
  write_10x_fragments_cpp(
    normalizePath(path, mustWork = FALSE),
    iterate_fragments(fragments),
    append_5th_column
  )

  open_fragments_10x(path, comment = "", end_inclusive = end_inclusive)
}


setClass("UnpackedMemFragments",
  contains = "IterableFragments",
  slots = c(
    cell = "integer",
    start = "integer",
    end = "integer",
    end_max = "integer",
    chr_ptr = "numeric",
    chr_names = "character",
    cell_names = "character",
    version = "character"
  ),
  prototype = list(
    cell = integer(0),
    start = integer(0),
    end = integer(0),
    end_max = integer(0),
    chr_ptr = numeric(0),
    chr_names = character(0),
    cell_names = character(0),
    version = character(0)
  )
)
setMethod("chrNames", "UnpackedMemFragments", function(x) x@chr_names)
setMethod("cellNames", "UnpackedMemFragments", function(x) x@cell_names)
setMethod("iterate_fragments", "UnpackedMemFragments", function(x) {
  iterate_unpacked_fragments_cpp(x)
})
setMethod("short_description", "UnpackedMemFragments", function(x) {
  "Read uncompressed fragments from memory"
})

setClass("PackedMemFragments",
  contains = "IterableFragments",
  slots = c(
    cell_data = "integer",
    cell_idx = "integer",
    cell_idx_offsets = "numeric",
    start_data = "integer",
    start_idx = "integer",
    start_idx_offsets = "numeric",
    start_starts = "integer",
    end_data = "integer",
    end_idx = "integer",
    end_idx_offsets = "numeric",
    end_max = "integer",
    chr_ptr = "numeric",
    chr_names = "character",
    cell_names = "character",
    version = "character"
  ),
  prototype = list(
    cell_data = integer(0),
    cell_idx = integer(0),
    cell_idx_offsets = numeric(0),
    start_data = integer(0),
    start_idx = integer(0),
    start_idx_offsets = numeric(0),
    start_starts = integer(0),
    end_data = integer(0),
    end_idx = integer(0),
    end_idx_offsets = numeric(0),
    end_max = integer(0),
    chr_ptr = numeric(0),
    chr_names = character(0),
    cell_names = character(0),
    version = character(0)
  )
)
setMethod("chrNames", "PackedMemFragments", function(x) x@chr_names)
setMethod("cellNames", "PackedMemFragments", function(x) x@cell_names)
setMethod("iterate_fragments", "PackedMemFragments", function(x) {
  iterate_packed_fragments_cpp(x)
})
setMethod("short_description", "PackedMemFragments", function(x) {
  "Read compressed fragments from memory"
})

#' Read/write BPCells fragment objects
#'
#' BPCells fragments can be read/written in compressed (bitpacked) or
#' uncompressed form in a variety of storage locations: in memory (as an R
#' object), in an hdf5 file, or in a directory on disk (containing binary
#' files).
#'
#' Saving in a directory on disk is a good default for local analysis, as it
#' provides the best I/O performance and lowest memory usage. The HDF5 format
#' allows saving within existing hdf5 files to group data together, and the in
#' memory format provides the fastest performance in the event memory usage is
#' unimportant.
#'
#' @param fragments Input fragments object
#' @param compress Whether or not to compress the data. With compression, storage size is
#' be about half the size of a gzip-compressed 10x fragments file.
#' @rdname fragment_io
#' @export
write_fragments_memory <- function(fragments, compress = TRUE) {
  assert_is(fragments, "IterableFragments")
  assert_is(compress, "logical")
  it <- iterate_fragments(fragments)
  if (compress) {
    res <- write_packed_fragments_cpp(it)
    do.call(new, c("PackedMemFragments", res))
  } else {
    res <- write_unpacked_fragments_cpp(it)
    do.call(new, c("UnpackedMemFragments", res))
  }
}

setClass("FragmentsDir",
  contains = "IterableFragments",
  slots = c(
    dir = "character",
    compressed = "logical",
    buffer_size = "integer",
    chr_names = "character",
    cell_names = "character"
  ),
  prototype = list(
    dir = character(0),
    compressed = TRUE,
    buffer_size = 1024L,
    chr_names = character(0),
    cell_names = character(0)
  )
)
setMethod("chrNames", "FragmentsDir", function(x) x@chr_names)
setMethod("cellNames", "FragmentsDir", function(x) x@cell_names)
setMethod("iterate_fragments", "FragmentsDir", function(x) {
  if (x@compressed) {
    iterate_packed_fragments_file_cpp(x@dir, x@buffer_size, x@chr_names, x@cell_names)
  } else {
    iterate_unpacked_fragments_file_cpp(x@dir, x@buffer_size, x@chr_names, x@cell_names)
  }
})
setMethod("short_description", "FragmentsDir", function(x) {
  sprintf(
    "Read %s fragments from directory %s",
    if (x@compressed) "compressed" else "uncompressed",
    x@dir
  )
})


#' @param dir Directory to read/write the data from
#' @param buffer_size For performance tuning only. The number of items to be bufferred
#' in memory before calling writes to disk.
#' @param overwrite If `TRUE`, write to a temp dir then overwrite existing data. Alternatively,
#'   pass a temp path as a string to customize the temp dir location.
#' @return Fragment object
#' @rdname fragment_io
#' @export
write_fragments_dir <- function(fragments, dir, compress = TRUE, buffer_size = 1024L, overwrite = FALSE) {
  assert_is(fragments, "IterableFragments")
  assert_is(dir, "character")
  assert_is(compress, "logical")
  assert_is(buffer_size, "integer")
  assert_is(overwrite, c("logical", "character"))
  if (is(overwrite, "character")) {
    assert_true(dir.exists(overwrite))
    overwrite_path <- tempfile("overwrite", tmpdir=overwrite)
    overwrite <- TRUE
  } else if (overwrite) {
    overwrite_path <- tempfile("overwrite")
  }

  dir <- path.expand(dir)
  did_tmp_copy <- FALSE
  if (overwrite && dir.exists(dir)) {
    fragments <- write_fragments_dir(fragments, overwrite_path, compress, buffer_size)
    did_tmp_copy <- TRUE
  }

  it <- iterate_fragments(fragments)
  if (compress) {
    write_packed_fragments_file_cpp(it, dir, buffer_size, overwrite)
  } else {
    write_unpacked_fragments_file_cpp(it, dir, buffer_size, overwrite)
  }

  if (did_tmp_copy) {
    unlink(overwrite_path, recursive=TRUE)
  }
  open_fragments_dir(dir, buffer_size)
}

#' @rdname fragment_io
#' @export
open_fragments_dir <- function(dir, buffer_size = 1024L) {
  assert_is_file(dir)
  assert_is(buffer_size, "integer")

  dir <- path.expand(dir)
  info <- info_fragments_file_cpp(dir, buffer_size)
  new("FragmentsDir", dir = path.expand(dir), compressed = info$compressed, buffer_size = buffer_size, cell_names = info$cell_names, chr_names = info$chr_names)
}

setClass("FragmentsHDF5",
  contains = "IterableFragments",
  slots = c(
    path = "character",
    group = "character",
    compressed = "logical",
    buffer_size = "integer",
    chr_names = "character",
    cell_names = "character"
  ),
  prototype = list(
    path = character(0),
    group = "",
    compressed = TRUE,
    buffer_size = 8192L,
    chr_names = character(0),
    cell_names = character(0)
  )
)
setMethod("chrNames", "FragmentsHDF5", function(x) x@chr_names)
setMethod("cellNames", "FragmentsHDF5", function(x) x@cell_names)
setMethod("iterate_fragments", "FragmentsHDF5", function(x) {
  if (x@compressed) {
    iterate_packed_fragments_hdf5_cpp(x@path, x@group, x@buffer_size, x@chr_names, x@cell_names)
  } else {
    iterate_unpacked_fragments_hdf5_cpp(x@path, x@group, x@buffer_size, x@chr_names, x@cell_names)
  }
})
setMethod("short_description", "FragmentsHDF5", function(x) {
  sprintf(
    "Read %s fragments from %s, group %s",
    if (x@compressed) "compressed" else "uncompressed",
    x@path,
    x@group
  )
})

#' @rdname fragment_io
#' @param path Path to the hdf5 file on disk
#' @param group The group within the hdf5 file to write the data to. If writing
#' to an existing hdf5 file this group must not already be in use
#' @param chunk_size For performance tuning only. The chunk size used for the HDF5 array storage.
#' @param gzip_level Gzip compression level. Default is 0 (no compression). This is recommended when both compression and
#'     compatibility with outside programs is required. Otherwise, using compress=TRUE is recommended
#'     as it is >10x faster with often similar compression levels. 
#' @export
write_fragments_hdf5 <- function(
    fragments, 
    path, 
    group = "fragments", 
    compress = TRUE, 
    buffer_size = 8192L, 
    chunk_size = 1024L, 
    overwrite = FALSE,
    gzip_level = 0L
) {
  assert_is(fragments, "IterableFragments")
  assert_is(path, "character")
  assert_is(group, "character")
  assert_is(compress, "logical")
  assert_is(buffer_size, "integer")
  assert_is(chunk_size, "integer")
  assert_is(gzip_level, "integer")
  assert_is(overwrite, c("logical", "character"))
  if (is(overwrite, "character")) {
    assert_true(dir.exists(overwrite))
    overwrite_path <- tempfile("overwrite", tmpdir=overwrite)
    overwrite <- TRUE
  } else if (overwrite) {
    overwrite_path <- tempfile("overwrite")
  }

  if (gzip_level != 0L && compress) {
     rlang::inform(c(
        "Warning: Mixing gzip compression (gzip_level > 0) with bitpacking compression (compress=TRUE) may be slower than bitpacking compression alone, with little space savings"
     ))
  }

  path <- path.expand(path)
  did_tmp_copy <- FALSE
  if (overwrite && hdf5_group_exists_cpp(path, group)) {
    rlang::inform(c(
      "Warning: Overwriting an hdf5 dataset does not free old storage"
    ), .frequency = "regularly", .frequency_id = "hdf5_overwrite")
    did_tmp_copy <- TRUE
    fragments <- write_fragments_dir(fragments, overwrite_path, compress, buffer_size)
  }
  it <- iterate_fragments(fragments)
  if (compress) {
    write_packed_fragments_hdf5_cpp(it, path, group, buffer_size, chunk_size, overwrite, gzip_level)
  } else {
    write_unpacked_fragments_hdf5_cpp(it, path, group, buffer_size, chunk_size, overwrite, gzip_level)
  }

  if (did_tmp_copy) {
    unlink(overwrite_path, recursive=TRUE)
  }
  open_fragments_hdf5(path, group, buffer_size)
}

#' @rdname fragment_io
#' @export
open_fragments_hdf5 <- function(path, group = "fragments", buffer_size = 16384L) {
  assert_is_file(path)
  assert_is(group, "character")
  assert_is(buffer_size, "integer")

  path <- path.expand(path)
  info <- info_fragments_hdf5_cpp(path, group, buffer_size)
  new("FragmentsHDF5", path = path, group = group, compressed = info$compressed, buffer_size = buffer_size, cell_names = info$cell_names, chr_names = info$chr_names)
}




#' Convert between BPCells fragments and R objects.
#'
#' BPCells fragments can be interconverted with GRanges and data.frame R objects.
#' The main conversion method is R's builtin `as()` function, though the 
#' `convert_to_fragments()` helper is also available. For all R objects except 
#' GRanges, BPCells assumes a 0-based, end-exclusive coordinate system. (See
#' [genomic-ranges-like] reference for details)
#'
#' @usage 
#' # Convert from R to BPCells
#' convert_to_fragments(x, zero_based_coords = !is(x, "GRanges"))
#' as(x, "IterableFragments")
#' 
#' # Convert from BPCells to R
#' as.data.frame(bpcells_fragments)
#' as(bpcells_fragments, "data.frame")
#' as(bpcells_fragments, "GRanges")
#'
#' @param x `r document_granges("Fragment coordinates", extras=c("cell_id" = "cell barcodes or unique identifiers as string or factor"))`
#' @param zero_based_coords Whether to convert the ranges from a 1-based end-inclusive
#'    coordinate system to a 0-based end-exclusive coordinate system. Defaults to true
#'    for GRanges and false for other formats
#'    (see this [archived UCSC blogpost](https://web.archive.org/web/20210920203703/http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/))
#' @return **convert_to_fragments()**: IterableFragments object
#' @rdname fragment_R_conversion
#' @export
convert_to_fragments <- function(x, zero_based_coords = !is(x, "GRanges")) {
  assert_is(zero_based_coords, "logical")
  x <- normalize_ranges(x, metadata_cols = "cell_id", zero_based_coords = zero_based_coords)
  x <- data.frame(x)
  x <- x[order(x$chr, x$start), ]

  x$cell_id <- as.factor(x$cell_id)
  x$chr <- as.factor(x$chr)

  chr_ptr <- rep(cumsum(table(x$chr)), each = 2)
  chr_ptr <- c(0, chr_ptr[-length(chr_ptr)])
  names(chr_ptr) <- NULL

  end_max <- calculate_end_max_cpp(as.integer(x$end), chr_ptr)
  new("UnpackedMemFragments",
    cell = as.integer(x$cell_id) - 1L,
    start = as.integer(x$start),
    end = as.integer(x$end),
    end_max = end_max,
    chr_ptr = chr_ptr,
    cell_names = levels(x$cell_id),
    chr_names = levels(x$chr),
    version = "unpacked-fragments-v2"
  )
}


# Define converters with GRanges if we have the GenomicRanges package available
if (requireNamespace("GenomicRanges", quietly = TRUE)) {
  setAs("IterableFragments", "GRanges", function(from) {
    if (is(from, "UnpackedMemFragments")) {
      frags <- from
    } else {
      frags <- write_fragments_memory(from, compress = FALSE)
    }
    # Get the differences between adjacent elements on the chr_ptr list
    # to calculate how many values are in each chromosome
    chr_counts <- frags@chr_ptr[c(FALSE, TRUE)] - frags@chr_ptr[c(TRUE, FALSE)]
    # order by where this chromosome range starts
    chr_order <- order(frags@chr_ptr[c(TRUE, FALSE)])

    GenomicRanges::GRanges(
      seqnames = S4Vectors::Rle(values = factor(chrNames(frags)[chr_order], chrNames(frags)), lengths = chr_counts[chr_order]),
      ranges = IRanges::IRanges(
        start = frags@start + 1,
        end = frags@end
      ),
      cell_id = factor(frags@cell_names[frags@cell + 1], frags@cell_names)
    )
  })
  setAs("GRanges", "IterableFragments", function(from) {
    convert_to_fragments(from)
  })
}

setAs("IterableFragments", "data.frame", function(from) {
  if (is(from, "UnpackedMemFragments")) {
    frags <- from
  } else {
    frags <- write_fragments_memory(from, compress = FALSE)
  }
  # Get the differences between adjacent elements on the chr_ptr list
  # to calculate how many values are in each chromosome
  chr_counts <- frags@chr_ptr[c(FALSE, TRUE)] - frags@chr_ptr[c(TRUE, FALSE)]
  # order by where this chromosome range starts
  chr_order <- order(frags@chr_ptr[c(TRUE, FALSE)])
  data.frame(
    chr = factor(rep.int(chrNames(frags)[chr_order], chr_counts[chr_order]), levels=chrNames(frags)),
    start = frags@start,
    end = frags@end,
    cell_id = factor(cellNames(frags)[frags@cell + 1], cellNames(frags))
  )
})

#' @exportS3Method base::as.data.frame
as.data.frame.IterableFragments <- function(x, ...) as(x, "data.frame")
#' @export
setMethod("as.data.frame", signature(x = "IterableFragments"), function(x, ...) as(x, "data.frame"))
setAs("data.frame", "IterableFragments", function(from) {
  convert_to_fragments(from)
}) 

# Shift start/end coordinates of fragments
setClass("ShiftFragments",
  contains = "IterableFragments",
  slots = c(
    fragments = "IterableFragments",
    shift_start = "integer",
    shift_end = "integer"
  ),
  prototype = list(
    fragments = NULL,
    shift_start = 0L,
    shift_end = 0L
  )
)
setMethod("iterate_fragments", "ShiftFragments", function(x) {
  iterate_shift_cpp(iterate_fragments(x@fragments), x@shift_start, x@shift_end)
})
setMethod("short_description", "ShiftFragments", function(x) {
  c(
    short_description(x@fragments),
    sprintf("Shift start %+dbp, end %+dbp", x@shift_start, x@shift_end)
  )
})
#' Shift start or end coordinates
#'
#' Shifts start or end of fragments by a fixed amount, which can be useful to
#' correct the Tn5 offset.
#'
#' The correct Tn5 offset is +/- 4bp since the Tn5 cut sites on opposite strands
#' are offset by 9bp. However, +4/-5 bp is often applied to bed-format files,
#' since the end coordinate in bed files is 1 past the last basepair of the sequenced
#' DNA fragment. This results in a bed-like format except with inclusive end coordinates.
#' @param fragments Input fragments object
#' @param shift_start How many basepairs to shift the start coords
#' @param shift_end How many basepairs to shift the end coords
#' @return Shifted fragments object
#' @export
shift_fragments <- function(fragments, shift_start = 0L, shift_end = 0L) {
  assert_is_wholenumber(shift_start)
  assert_is_wholenumber(shift_end)
  assert_len(shift_start, 1)
  assert_len(shift_end, 1)
  new("ShiftFragments", fragments = fragments, shift_start = as.integer(shift_start), shift_end = as.integer(shift_end))
}

# Select fragments by length
setClass("SelectLength",
  contains = "IterableFragments",
  slots = c(
    fragments = "IterableFragments",
    min_len = "integer",
    max_len = "integer"
  ),
  prototype = list(
    fragments = NULL,
    min_len = NA_integer_,
    max_len = NA_integer_
  )
)
setMethod("iterate_fragments", "SelectLength", function(x) {
  max_len <- if (is.na(x@max_len)) .Machine$integer.max else x@max_len
  iterate_length_select_cpp(iterate_fragments(x@fragments), x@min_len, max_len)
})
setMethod("short_description", "SelectLength", function(x) {
  text_min <- if (!is.na(x@min_len)) sprintf("%d", x@min_len) else "0"
  text_max <- if (!is.na(x@max_len)) sprintf("%d", x@max_len) else "Inf"
  c(
    short_description(x@fragments),
    sprintf("Filter to fragment sizes %s bp-%s bp", x@min_len, x@max_len)
  )
})
#' Subset fragments by length
#' @param fragments Input fragments object
#' @param min_len Minimum bases in fragment (inclusive)
#' @param max_len Maximum bases in fragment (inclusive)
#' @return Fragments object
#' @details Fragment length is calculated as `end-start`
#' @export
subset_lengths <- function(fragments, min_len = 0L, max_len = NA_integer_) {
  assert_is_wholenumber(min_len)
  if (!is.na(max_len)) assert_is_wholenumber(max_len)
  assert_len(min_len, 1)
  assert_len(max_len, 1)
  new("SelectLength", fragments = fragments, min_len = as.integer(min_len), max_len = as.integer(max_len))
}


# Select fragments by chromosome
setClass("ChrSelectName",
  contains = "IterableFragments",
  slots = c(
    fragments = "IterableFragments",
    chr_names = "character"
  ),
  prototype = list(
    fragments = NULL,
    chr_names = character(0)
  )
)
setMethod("chrNames", "ChrSelectName", function(x) x@chr_names)
setMethod("iterate_fragments", "ChrSelectName", function(x) {
  iterate_chr_name_select_cpp(iterate_fragments(x@fragments), x@chr_names)
})
setMethod("short_description", "ChrSelectName", function(x) {
  c(
    short_description(x@fragments),
    sprintf(
      "Select %d chromosomes by name%s",
      length(x@chr_names),
      pretty_print_vector(x@chr_names, max_len = 3, prefix = ": ")
    )
  )
})
setClass("ChrSelectIndex",
  contains = "IterableFragments",
  slots = c(
    fragments = "IterableFragments",
    chr_index_selection = "integer"
  ),
  prototype = list(
    fragments = NULL,
    chr_index_selection = NA_integer_
  )
)
setMethod("chrNames", "ChrSelectIndex", function(x) chrNames(x@fragments)[x@chr_index_selection])
setMethod("iterate_fragments", "ChrSelectIndex", function(x) {
  iterate_chr_index_select_cpp(iterate_fragments(x@fragments), x@chr_index_selection - 1)
})
setMethod("short_description", "ChrSelectIndex", function(x) {
  c(
    short_description(x@fragments),
    sprintf(
      "Select %d chromosomes by index%s",
      length(x@chr_index_selection),
      pretty_print_vector(x@chr_index_selection, max_len = 3, prefix = ": ")
    )
  )
})

#' Subset, translate, or reorder chromosome IDs
#' 
#' @param fragments Input fragments object
#' @param chromosome_selection List of chromosme IDs (numeric), or names (character).
#'    The output chromosome ID `n` will be taken
#'    from the input fragments chromosome with ID/name `chromosome_selection[n]`.
#' @export
select_chromosomes <- function(fragments, chromosome_selection) {
  assert_is(fragments, "IterableFragments")
  assert_distinct(chromosome_selection)
  assert_is(chromosome_selection, c("character", "numeric"))
  if (is.numeric(chromosome_selection)) {
    assert_greater_than_zero(chromosome_selection)
    new_names <- NA_character_
    if (!is.null(chrNames(fragments))) {
      assert_true(all(chromosome_selection <= length(chrNames(fragments))))
    }
    return(
      new("ChrSelectIndex", fragments = fragments, chr_index_selection = as.integer(chromosome_selection))
    )
  } else {
    return(
      new("ChrSelectName", fragments = fragments, chr_names = chromosome_selection)
    )
  }
}



# Select fragments by cell
setClass("CellSelectName",
  contains = "IterableFragments",
  slots = c(
    fragments = "IterableFragments",
    cell_names = "character"
  ),
  prototype = list(
    fragments = NULL,
    cell_names = character(0)
  )
)
setMethod("cellNames", "CellSelectName", function(x) x@cell_names)
setMethod("iterate_fragments", "CellSelectName", function(x) {
  iterate_cell_name_select_cpp(iterate_fragments(x@fragments), x@cell_names)
})
setMethod("short_description", "CellSelectName", function(x) {
  c(
    short_description(x@fragments),
    sprintf(
      "Select %d cells by name%s",
      length(x@cell_names),
      pretty_print_vector(x@cell_names, max_len = 3, prefix = ": ")
    )
  )
})
setClass("CellSelectIndex",
  contains = "IterableFragments",
  slots = c(
    fragments = "IterableFragments",
    cell_index_selection = "integer"
  ),
  prototype = list(
    fragments = NULL,
    cell_index_selection = NA_integer_
  )
)
setMethod("cellNames", "CellSelectIndex", function(x) cellNames(x@fragments)[x@cell_index_selection])
setMethod("iterate_fragments", "CellSelectIndex", function(x) {
  iterate_cell_index_select_cpp(iterate_fragments(x@fragments), x@cell_index_selection - 1)
})
setMethod("short_description", "CellSelectIndex", function(x) {
  c(
    short_description(x@fragments),
    sprintf(
      "Select %d cells by index%s",
      length(x@cell_index_selection),
      pretty_print_vector(x@cell_index_selection, max_len = 3, prefix = ": ")
    )
  )
})
#' Subset, translate, or reorder cell IDs
#' @param fragments Input fragments object
#' @param cell_selection List of chromosme IDs (numeric), or names (character).
#'    The output cell ID `n` will be taken
#'    from the input cell with ID/name `cell_selection[n]`.
#' @export
select_cells <- function(fragments, cell_selection) {
  assert_is(fragments, "IterableFragments")
  assert_distinct(cell_selection)
  assert_is(cell_selection, c("character", "numeric"))
  if (is.numeric(cell_selection)) {
    assert_greater_than_zero(cell_selection)
    new_names <- NA_character_
    if (!is.null(cellNames(fragments))) {
      assert_true(all(cell_selection <= length(cellNames(fragments))))
    }
    return(
      new("CellSelectIndex", fragments = fragments, cell_index_selection = as.integer(cell_selection))
    )
  } else {
    return(
      new("CellSelectName", fragments = fragments, cell_names = cell_selection)
    )
  }
}

setClass("CellMerge",
  contains = "IterableFragments",
  slots = c(
    fragments = "IterableFragments",
    group_ids = "integer",
    group_names = "character"
  ),
  prototype = list(
    fragments = NULL,
    group_ids = integer(0),
    group_names = character(0)
  )
)
setMethod("cellNames", "CellMerge", function(x) x@group_names)
setMethod("iterate_fragments", "CellMerge", function(x) {
  iterate_cell_merge_cpp(iterate_fragments(x@fragments), x@group_ids, x@group_names)
})
setMethod("short_description", "CellMerge", function(x) {
  c(
    short_description(x@fragments),
    sprintf(
      "Merge %d cells into %d groups",
      length(x@group_ids),
      length(x@group_names)
    )
  )
})
#' Merge cells into pseudobulks
#'
#' Peak and tile matrix calculations can be sped up by reducing the number of
#' cells. For cases where the outputs are going to be added together afterwards, this can
#' provide a performance improvement
#'
#' @param fragments Input fragments object
#' @param cell_groups Character or factor vector providing a group for each cell.
#' Ordering is the same as `cellNames(fragments)`
#' @export
merge_cells <- function(fragments, cell_groups) {
  assert_is(fragments, "IterableFragments")
  assert_not_null(cellNames(fragments))
  assert_is(cell_groups, c("character", "factor"))
  cell_groups <- as.factor(cell_groups)
  new("CellMerge", fragments = fragments, group_ids = as.integer(cell_groups) - 1L, group_names = levels(cell_groups))
}

# Rename cells or chromosomes
# Select fragments by chromosome
setClass("ChrRename",
  contains = "IterableFragments",
  slots = c(
    fragments = "IterableFragments",
    chr_names = "character"
  ),
  prototype = list(
    fragments = NULL,
    chr_names = character(0)
  )
)
setMethod("chrNames", "ChrRename", function(x) x@chr_names)
setMethod("iterate_fragments", "ChrRename", function(x) {
  iterate_chr_rename_cpp(iterate_fragments(x@fragments), x@chr_names)
})
setMethod("short_description", "ChrRename", function(x) {
  c(
    short_description(x@fragments),
    sprintf("Rename chromosomes to: %s", pretty_print_vector(x@chr_names, max_len = 3))
  )
})

# Rename chr/cell names for cases where names are not all known ahead-of-time
setClass("CellRename",
  contains = "IterableFragments",
  slots = c(
    fragments = "IterableFragments",
    cell_names = "character"
  ),
  prototype = list(
    fragments = NULL,
    cell_names = character(0)
  )
)
setMethod("cellNames", "CellRename", function(x) x@cell_names)
setMethod("iterate_fragments", "CellRename", function(x) {
  iterate_cell_rename_cpp(iterate_fragments(x@fragments), x@cell_names)
})
setMethod("short_description", "CellRename", function(x) {
  c(
    short_description(x@fragments),
    sprintf("Rename cells to: %s", pretty_print_vector(x@cell_names, max_len = 3))
  )
})

setClass("CellPrefix",
  contains = "IterableFragments",
  slots = c(
    fragments = "IterableFragments",
    prefix = "character"
  ),
  prototype = list(
    fragments = NULL,
    prefix = ""
  )
)
setMethod("cellNames", "CellPrefix", function(x) {
  if (is.null(cellNames(x@fragments))) return(NULL)
  paste0(x@prefix, cellNames(x@fragments))
})
setMethod("iterate_fragments", "CellPrefix", function(x) {
  iterate_cell_prefix_cpp(iterate_fragments(x@fragments), x@prefix)
})
setMethod("short_description", "CellPrefix", function(x) {
  c(
    short_description(x@fragments),
    sprintf("Prefix cells names with: %s", x@prefix)
  )
})

#' Add sample prefix to cell names
#' 
#' Rename cells by adding a prefix to the names. Most commonly this will be a sample name.
#' All cells will recieve the exact text of `prefix` added to the beginning, so any separator
#' characters like "_" must be included in the given `prefix`. Use this prior to merging
#' fragments from different experiments with `c()` in order to help prevent cell name
#' clashes.
#' @param fragments Input fragments object.
#' @param prefix String to add as the prefix
#' @return Fragments object with prefixed names
#' @export
prefix_cell_names <- function(fragments, prefix) {
  assert_is(fragments, "IterableFragments")
  assert_is(prefix, "character")
  assert_len(prefix, 1)
  new("CellPrefix", fragments = fragments, prefix = prefix)
}

# Select fragments by chromosome
setClass("RegionSelect",
  contains = "IterableFragments",
  slots = c(
    fragments = "IterableFragments",
    chr_id = "integer",
    start = "integer",
    end = "integer",
    chr_levels = "character",
    invert_selection = "logical"
  ),
  prototype = list(
    fragments = NULL,
    chr_id = integer(0),
    start = integer(0),
    end = integer(0),
    chr_levels = character(0),
    invert_selection = FALSE
  )
)
setMethod("iterate_fragments", "RegionSelect", function(x) {
  iterate_region_select_cpp(iterate_fragments(x@fragments), x@chr_id, x@start, x@end, x@chr_levels, x@invert_selection)
})
setMethod("short_description", "RegionSelect", function(x) {
  # Subset strings first to avoid a very slow string concatenation process
  indices <- c(head(seq_along(x@chr_id), 3), tail(seq_along(x@chr_id), 1))
  labels <- paste0(x@chr_levels[1 + x@chr_id[indices]], ":", x@start[indices] + 1, "-", x@end[indices])

  c(
    short_description(x@fragments),
    sprintf(
      "Subset to fragments %soverlapping %d ranges%s",
      ifelse(x@invert_selection, "not ", ""),
      length(x@chr_id),
      pretty_print_vector(labels, prefix = ": ", max_len = 2)
    )
  )
})

#' Subset fragments by genomic region
#'
#' Fragments can be subset based on overlapping (or not overlapping) a set of regions
#' @param fragments Input fragments object.
#' @param invert_selection If TRUE, select fragments *not* overlapping selected regions
#'  instead of only fragments overlapping the selected regions.
#' @inheritParams peak_matrix
#' @return Fragments object filtered according to the selected regions
#' @export
select_regions <- function(fragments, ranges, invert_selection = FALSE, zero_based_coords = !is(ranges, "GRanges")) {
  assert_is(fragments, "IterableFragments")
  ranges <- normalize_ranges(ranges, zero_based_coords = zero_based_coords)

  assert_is(zero_based_coords, "logical")

  chr_levels <- levels(ranges$chr)
  chr_id <- as.integer(ranges$chr) - 1L

  new("RegionSelect", fragments = fragments, chr_id = chr_id, start = ranges$start, end = ranges$end, chr_levels = chr_levels, invert_selection = invert_selection)
}

setClass("MergeFragments",
  contains = "IterableFragments",
  slots = c(
    fragments_list = "list"
  ),
  prototype = list(
    fragments_list = list()
  )
)
setMethod("chrNames", "MergeFragments", function(x) {
  Reduce(union, lapply(x@fragments_list, chrNames))
})
setMethod("cellNames", "MergeFragments", function(x) {
  do.call(c, lapply(x@fragments_list, cellNames))
})

setMethod("iterate_fragments", "MergeFragments", function(x) {
  iterate_merge_fragments_cpp(lapply(x@fragments_list, iterate_fragments), chrNames(x))
})
setMethod("short_description", "MergeFragments", function(x) {
  # Subset strings first to avoid a very slow string concatenation process
  sprintf(
    "Merge %d fragments objects of classes%s",
    length(x@fragments_list),
    pretty_print_vector(vapply(x@fragments_list, class, character(1)), prefix = ": ", max_len = 3)
  )
})

# Allow merging fragments using standard concatenation method
setMethod("c", "IterableFragments", function(x, ...) {
  tail <- list(...)
  if (length(tail) == 1 && is.null(tail[[1]])) {
    args <- list(x)
  } else {
    args <- c(list(x), list(...))
  }
  chr_names <- chrNames(x)
  assert_not_null(chr_names)
  fragments_list <- list()
  for (i in seq_along(args)) {
    assert_is(args[[i]], "IterableFragments")
    assert_not_null(chrNames(args[[i]]))
    if (is(args[[i]], "MergeFragments")) {
      fragments_list <- c(fragments_list, args[[i]]@fragments_list)
    } else {
      fragments_list <- c(fragments_list, args[[i]])
    }
  }
  if (length(fragments_list) == 1) {
    return(fragments_list)[[1]]
  }
  res <- new("MergeFragments", fragments_list = fragments_list)
  if (anyDuplicated(cellNames(res))) {
    rlang::inform(c("Warning: duplicicate cell names detected when merging fragments.", "Try using prefix_cell_names() to disambiguate"))
  }
  return(res)
})


#' Check if two fragments objects are identical
#' @param fragments1 First IterableFragments to compare
#' @param fragments2 Second IterableFragments to compare
#' @return boolean for whether the fragments objects are identical
fragments_identical <- function(fragments1, fragments2) {
  assert_is(fragments1, "IterableFragments")
  assert_is(fragments2, "IterableFragments")
  i1 <- iterate_fragments(fragments1)
  i2 <- iterate_fragments(fragments2)
  fragments_identical_cpp(i1, i2)
}

#' Scan through fragments without performing any operations (used for benchmarking)
#' @param fragments Fragments object to scan
#' @return Length 4 vector with fragment count, then sums of chr, starts, and ends
#' @keywords internal
scan_fragments <- function(fragments) {
  assert_is(fragments, "IterableFragments")
  i <- iterate_fragments(fragments)
  scan_fragments_cpp(i)
}


pretty_print_vector <- function(x, max_len = 3, sep = ", ", prefix = "", empty = "") {
  # Print whole vector if possible, or else
  # print first few, ellipsis, then the last one
  # if x is length 0, just return the empty value
  if (length(x) == 0) {
    return(empty)
  }

  if (length(x) <= max_len) {
    res <- paste0(x, collapse = sep)
  } else {
    res <- paste0(x[seq_len(max_len - 1)], collapse = sep)
    res <- paste0(res, " ... ", x[length(x)], collapse = sep)
  }

  res <- paste0(prefix, res)
  res
}
