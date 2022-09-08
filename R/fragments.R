#' IterableFragments methods
#'
#' Methods for IterableFragments objects
#'
#' @name IterableFragments-methods
#' @rdname IterableFragments-methods
NULL

setClass("IterableFragments")

#' This class is used to hold a list of external pointers.
setClass("XPtrList", slots=c("pointers"="list", "type"="character"))
ptr <- function(x) {x@pointers[[1]]}
#' Wrap an inner XPtrList with an outer externalptr object (aka Rcpp::XPtr)
#' We also track the type of the current pointer to provide a little safety when converting
#' to C++.
#' WARNING: Always make sure to pass the correct inner object, because if it is not passed then
#' the chain of pointers will be broken and C++ objects could be freed at random by the GC while we're
#' still using them. ONLY pass inner=new("XPtrList") if the C++ object is not wrapping any other fragments/matrices
wrapFragments <- function(outer, inner) {
    inner@pointers <- c(outer, inner@pointers)
    inner@type <- "fragments"
    inner
}
wrapMat_uint32_t <- function(outer, inner) {
    inner@pointers <- c(outer, inner@pointers)
    inner@type <- "mat_uint32_t"
    inner
}
wrapMat_double <- function(outer, inner) {
    inner@pointers <- c(outer, inner@pointers)
    inner@type <- "mat_double"
    inner
}
wrapMat_float <- function(outer, inner) {
    inner@pointers <- c(outer, inner@pointers)
    inner@type <- "mat_float"
    inner
}
setMethod("show", "XPtrList", function(object) {
    cat(sprintf("C++ XPtrList with %d pointers\n", length(object@pointers)))
})

# Get an external pointer to a C++ FragmentsLoader object. This is only called
# when an operation is actively ready to run
setGeneric("iterate_fragments", function(x) standardGeneric("iterate_fragments"))
setMethod("iterate_fragments", "XPtrList", function(x) {
    stopifnot(x@type == "fragments")
    x
})

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

    if(is.null(cellNames(object)))
        cat("Cells: count unknown\n")
    else
        cat(sprintf("Cells: %d cells with %s\n", length(cellNames(object)),
            pretty_print_vector(cellNames(object), prefix="names ", empty="unknown names")
        ))
    
    if(is.null(chrNames(object)))
        cat("Chromosomes: count unknown\n")
    else
        cat(sprintf("Chromosomes: %d chromosomes with %s\n", length(chrNames(object)),
            pretty_print_vector(chrNames(object), prefix="names ", empty="unknown names")
        ))

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
    if (.hasSlot(x, "fragments"))
        return(cellNames(x@fragments))
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
    if (.hasSlot(x, "fragments")) {
        cellNames(x@fragments) <- value
        return(x)
    }
    stop(sprintf("Error: cellNames<- not implemented on object of class %s", class(x)))
})

#' Get chromosome names
#' @param x an IterableFragments object
#' @return * `chrNames()`: Character vector of chromosome names, or NULL if none are known
#' @describeIn IterableFragments-methods Set chromosome names
#' @export
setGeneric("chrNames", function(x) standardGeneric("chrNames"))
setMethod("chrNames", "IterableFragments", function(x) {
    if (.hasSlot(x, "fragments"))
        return(chrNames(x@fragments))
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
    assert_is_character(value)
    assert_len(value, length(chrNames(x)))
    if (.hasSlot(x, "fragments")) {
        chrNames(x@fragments) <- value
        return(x)
    }
    stop(sprintf("Error: chrNames<- not implemented on object of class %s", class(x)))
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

setMethod("chrNames<-", "FragmentsTsv", function(x, ..., value) {
    new("ChrRename", x, chr_names=value)
})
setMethod("cellNames<-", "FragmentsTsv", function(x, ..., value) {
    new("CellRename", x, cell_names=value)
})

setMethod("iterate_fragments", "FragmentsTsv", function(x) wrapFragments(iterate_10x_fragments_cpp(normalizePath(x@path), x@comment), new("XPtrList")))
setMethod("short_description", "FragmentsTsv", function(x) {
    sprintf("Load 10x fragments file from %s", x@path)
})

#' Read a 10x fragments file
#' @details Note: no disk operations will take place until the fragments are used in a function
#' @param path Path of 10x fragments file
#' @param comment Skip lines at beginning of file which start with comment
#' @param end_inclusive Whether the end coordinate of the bed is inclusive -- i.e. there was an
#'     insertion at the end coordinate rather than the base before the end coordinate. This is the
#'     10x default, though it's not quite standard for the bed file format.
#' @return 10x fragments file object
#' @export
open_fragments_10x <- function(path, comment="#", end_inclusive=TRUE) {
    assert_is_file(path, extension=c(".tsv", ".tsv.gz"))
    assert_is_character(comment)
    assert_len(comment, 1)
    path <- normalizePath(path)
    res <- new("FragmentsTsv", path=path, comment=comment)
    if (end_inclusive)
        res <- shift_fragments(res, shift_end=1)
    res
}


#' Write to a 10x fragments file
#' @param fragments Input fragments object
#' @inheritParams open_fragments_10x
#' @param append_5th_column Whether to include 5th column of all 0 for compatibility
#'        with 10x fragment file outputs (defaults to 4 columns chr,start,end,cell)
#' @export
write_fragments_10x <- function(fragments, path, end_inclusive=TRUE, append_5th_column=FALSE) {
    assert_is_file(path, must_exist=FALSE, extension=c(".tsv", ".tsv.gz"))
    if (end_inclusive) 
        fragments <- shift_fragments(fragments, shift_end=-1)
    p <- ptr(iterate_fragments(fragments))
    write_10x_fragments_cpp(
        normalizePath(path, mustWork=FALSE), 
        p,
        append_5th_column
    )

    open_fragments_10x(path, comment="", end_inclusive=end_inclusive)
}


setClass("UnpackedMemFragments",
    contains = "IterableFragments",
    slots = c(
        cell = "integer",
        start = "integer", 
        end = "integer", 
        end_max = "integer", 
        chr_ptr = "integer",
        chr_names = "character",
        cell_names = "character",
        version = "character"
    ),
    prototype = list(
        cell = integer(0),
        start = integer(0),
        end = integer(0),
        end_max = integer(0),
        chr_ptr = integer(0),
        chr_names = character(0),
        cell_names = character(0),
        version = character(0)
    )
)
setMethod("chrNames", "UnpackedMemFragments", function(x) x@chr_names)
setMethod("cellNames", "UnpackedMemFragments", function(x) x@cell_names)
setMethod("chrNames<-", "UnpackedMemFragments", function(x, ..., value) {
    assert_is_character(value)
    assert_len(value, length(x@chr_names))
    x@chr_names <- value
    x
})
setMethod("cellNames<-", "UnpackedMemFragments", function(x, ..., value) {
    assert_is_character(value)
    assert_len(value, length(x@cell_names))
    x@cell_names <- value
    x
})

setMethod("iterate_fragments", "UnpackedMemFragments", function(x) {
    wrapFragments(iterate_unpacked_fragments_cpp(x), new("XPtrList"))
})
setMethod("short_description", "UnpackedMemFragments", function(x) {
    "Read uncompressed fragments from memory"
})

setClass("PackedMemFragments",
    contains = "IterableFragments",
    slots = c(
        cell_data = "integer",
        cell_idx = "integer",
        start_data = "integer",
        start_idx = "integer",
        start_starts = "integer",
        end_data = "integer",
        end_idx = "integer",
        end_max = "integer",
        chr_ptr = "integer",
        
        chr_names = "character",
        cell_names = "character",
        version = "character"
    ),
    prototype = list(
        cell_data = integer(0),
        cell_idx = integer(0),
        start_data = integer(0),
        start_idx = integer(0),
        start_starts = integer(0),
        end_data = integer(0),
        end_idx = integer(0),
        end_max = integer(0),
        chr_ptr = integer(0),

        chr_names = character(0),
        cell_names = character(0),
        version = character(0)
    )
)
setMethod("chrNames", "PackedMemFragments", function(x) x@chr_names)
setMethod("cellNames", "PackedMemFragments", function(x) x@cell_names)
setMethod("chrNames<-", "PackedMemFragments", function(x, ..., value) {
    assert_is_character(value)
    assert_len(value, length(x@chr_names))
    x@chr_names <- value
    x
})
setMethod("cellNames<-", "PackedMemFragments", function(x, ..., value) {
    assert_is_character(value)
    assert_len(value, length(x@cell_names))
    x@cell_names <- value
    x
})

setMethod("iterate_fragments", "PackedMemFragments", function(x) {
    wrapFragments(iterate_packed_fragments_cpp(x), new("XPtrList"))
})
setMethod("short_description", "PackedMemFragments", function(x) {
    "Read compressed fragments from memory"
})

#' Make a fragments object in memory. 
#' @param fragments Input fragments object
#' @param compress Whether or not to compress the data. With compression, memory usage should
#' be about half the size of a gzip-compressed 10x fragments file holding the same data. 
#' Without compression memory usage can be 10-20% larger.
#' @return MemFragments object
#' @export
write_fragments_memory <- function(fragments, compress=TRUE) {
    assert_is(fragments, "IterableFragments")
    assert_is(compress, "logical")
    p <- ptr(iterate_fragments(fragments))
    if (compress) {
        res <- write_packed_fragments_cpp(p)
        do.call(new, c("PackedMemFragments", res))
    } else {
        res <- write_unpacked_fragments_cpp(p)
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
setMethod("chrNames<-", "FragmentsDir", function(x, ..., value) {
    assert_is_character(value)
    assert_len(value, length(x@chr_names))
    x@chr_names <- value
    x
})
setMethod("cellNames<-", "FragmentsDir", function(x, ..., value) {
    assert_is_character(value)
    assert_len(value, length(x@cell_names))
    x@cell_names <- value
    x
})
setMethod("iterate_fragments", "FragmentsDir", function(x) {
    if (x@compressed)
        wrapFragments(iterate_packed_fragments_file_cpp(x@dir, x@buffer_size, x@chr_names, x@cell_names), new("XPtrList"))
    else
        wrapFragments(iterate_unpacked_fragments_file_cpp(x@dir, x@buffer_size, x@chr_names, x@cell_names), new("XPtrList"))
})
setMethod("short_description", "FragmentsDir", function(x) {
    sprintf("Read %s fragments from directory %s", 
        if(x@compressed) "compressed" else "uncompressed",
        x@dir
    )
})
#' Make a fragments object in binary files within a directory. 
#' @inheritParams write_fragments_memory
#' @param dir Directory to save the data into
#' @param buffer_size For performance tuning only. The number of items to be bufferred
#' in memory before calling writes to disk.
#' @return FragmentsDir object
#' @export
write_fragments_dir <- function(fragments, dir, compress=TRUE, buffer_size=1024L) {
    assert_is(fragments, "IterableFragments")
    assert_is(dir, "character")
    assert_is(compress, "logical")
    assert_is(buffer_size, "integer")
    dir <- path.expand(dir)
    p <- ptr(iterate_fragments(fragments))
    if (compress)
        write_packed_fragments_file_cpp(p, dir, buffer_size)
    else
        write_unpacked_fragments_file_cpp(p, dir, buffer_size)
    open_fragments_dir(dir, buffer_size)
}

#' Open an existing fragments object written with [write_fragments_dir]
#' @param dir Directory to read the data from
#' @param buffer_size For performance tuning only. The number of fragments to be buffered
#' in memory for each read
#' @return FragmentsDir object
#' @export
open_fragments_dir <- function(dir, buffer_size=1024L) {
    assert_is_file(dir)
    assert_is(buffer_size, "integer")
    
    dir <- path.expand(dir)
    info <- info_fragments_file_cpp(dir, buffer_size)
    new("FragmentsDir", dir=path.expand(dir), compressed=info$compressed, buffer_size=buffer_size, cell_names=info$cell_names, chr_names=info$chr_names)
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
        compressed=TRUE,
        buffer_size = 8192L,

        chr_names = character(0),
        cell_names = character(0)
    )
)
setMethod("chrNames", "FragmentsHDF5", function(x) x@chr_names)
setMethod("cellNames", "FragmentsHDF5", function(x) x@cell_names)
setMethod("chrNames<-", "FragmentsHDF5", function(x, ..., value) {
    assert_is_character(value)
    assert_len(value, length(x@chr_names))
    x@chr_names <- value
    x
})
setMethod("cellNames<-", "FragmentsHDF5", function(x, ..., value) {
    assert_is_character(value)
    assert_len(value, length(x@cell_names))
    x@cell_names <- value
    x
})
setMethod("iterate_fragments", "FragmentsHDF5", function(x) {
    if (x@compressed)
        wrapFragments(iterate_packed_fragments_hdf5_cpp(x@path, x@group, x@buffer_size, x@chr_names, x@cell_names), new("XPtrList"))
    else
        wrapFragments(iterate_unpacked_fragments_hdf5_cpp(x@path, x@group, x@buffer_size, x@chr_names, x@cell_names), new("XPtrList"))
})
setMethod("short_description", "FragmentsHDF5", function(x) {
    sprintf("Read %s fragments from %s, group %s", 
        if(x@compressed) "compressed" else "uncompressed",
        x@path, 
        x@group
    )
})
#' Write a fragments object in an hdf5 file
#' @inheritParams write_fragments_dir
#' @param path Path to the hdf5 file on disk
#' @param group The group within the hdf5 file to write the data to. If writing
#' to an existing hdf5 file this group must not already be in use
#' @param chunk_size For performance tuning only. The chunk size used for the HDF5 array storage.
#' @return FragmentsHDF5 object
#' @export
write_fragments_hdf5 <- function(fragments, path, group="fragments", compress=TRUE, buffer_size=8192L, chunk_size=1024L) {
    assert_is(fragments, "IterableFragments")
    assert_is(path, "character")
    assert_is(group, "character")
    assert_is(compress, "logical")
    assert_is(buffer_size, "integer")
    assert_is(chunk_size, "integer")
    path <- path.expand(path)
    p <- ptr(iterate_fragments(fragments))
    if (compress)
        write_packed_fragments_hdf5_cpp(p, path, group, buffer_size, chunk_size)
    else
        write_unpacked_fragments_hdf5_cpp(p, path, group, buffer_size, chunk_size)
    open_fragments_hdf5(path, group, buffer_size)
}

#' Read a fragments object from an hdf5 file
#' @inheritParams write_fragments_hdf5
#' @param buffer_size For performance tuning only. The number of items to be bufferred
#' in memory before calling reads from disk.
#' @return FragmentsHDF5 object
#' @export
open_fragments_hdf5 <- function(path, group="fragments", buffer_size=16384L) {
    assert_is_file(path)
    assert_is(group, "character")
    assert_is(buffer_size, "integer")
    
    path <- path.expand(path)
    info <- info_fragments_hdf5_cpp(path, group, buffer_size)
    new("FragmentsHDF5", path=path, group=group, compressed=info$compressed, buffer_size=buffer_size, cell_names=info$cell_names, chr_names=info$chr_names)
}




#' Build a Fragments object from an R data frame or GRanges object
#' @param x An input GRanges, list or data frame. Lists and dataframes 
#'    must have chr, start, end, cell_id. GRanges must have a metadata column
#'    for cell_id
#' @param zero_based_coords Whether to convert the ranges from a 1-based end-inclusive
#'    coordinate system to a 0-based end-exclusive coordinate system. Defaults to true
#'    for GRanges and false for other formats
#'    (see http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/)
#' @return UnpackedMemFragments object representing the given fragments
#' @export
convert_to_fragments <- function(x, zero_based_coords=!is(x, "GRanges")) {
    assert_is(x, c("list", "data.frame", "GRanges"))
    assert_is(zero_based_coords, "logical")
    assert_not_null(x$cell_id)
    if(is(x, "GRanges")) {
        assert_has_names(GenomicRanges::mcols(x), c("cell_id"))
        x <- as.data.frame(x)
        x$chr <- x$seqnames
    } else {
        assert_has_names(x, c("chr", "start", "end", "cell_id"))
    }
    x <- x[order(x$chr, x$start),]
    if(!zero_based_coords) {
        x$start <- x$start - 1
    }
    x$cell_id <- as.factor(x$cell_id)
    x$chr <- as.factor(x$chr)

    chr_ptr <- rep(cumsum(table(x$chr)), each=2)
    chr_ptr <- c(0, chr_ptr[-length(chr_ptr)])

    end_max <- calculate_end_max_cpp(as.integer(x$end), chr_ptr)
    new("UnpackedMemFragments", 
        cell=as.integer(x$cell_id) - 1L, 
        start=as.integer(x$start),
        end=as.integer(x$end),
        end_max=end_max,
        chr_ptr=as.integer(chr_ptr),
        cell_names = levels(x$cell_id),
        chr_names = levels(x$chr), 
        version = "unpacked-fragments-v1"
    )
}


# Define converters with GRanges if we have the GenomicRanges package available
if (requireNamespace("GenomicRanges", quietly=TRUE)) {
    setAs("IterableFragments", "GRanges", function(from) {
        if (is(from, "UnpackedMemFragments")) {
            frags <- from
        } else {
            frags <- write_fragments_memory(from, compress=FALSE)
        }
        # Get the differences between adjacent elements on the chr_ptr list
        # to calculate how many values are in each chromosome
        chr_counts <- frags@chr_ptr[c(FALSE,TRUE)] - frags@chr_ptr[c(TRUE,FALSE)]
        # order by where this chromosome range starts
        chr_order <- order(frags@chr_ptr[c(TRUE,FALSE)])
        
        GenomicRanges::GRanges(
            seqnames=S4Vectors::Rle(values=factor(chrNames(frags)[order(chr_order)], chrNames(frags)), lengths=chr_counts),
            ranges = IRanges::IRanges(
                start = frags@start + 1,
                end = frags@end
            ),
            cell_id = factor(frags@cell_names[frags@cell + 1], frags@cell_names)
        )
    })
}

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
    inner <- iterate_fragments(x@fragments)
    wrapFragments(
        iterate_shift_cpp(ptr(inner), x@shift_start, x@shift_end), 
        inner
    )
})
setMethod("short_description", "ShiftFragments", function(x) {
    c(
        short_description(x@fragments),
        sprintf("Shift start %+dbp, end %+dbp", x@shift_start, x@shift_end)
    )
})
#' Shift start or end coordinates of fragments by a fixed amount
#' @param fragments Input fragments object
#' @param shift_start How many basepairs to shift the start coords
#' @param shift_end How many basepairs to shift the end coords
#' @return Shifted fragments object
#' @export
shift_fragments <- function(fragments, shift_start=0L, shift_end=0L) {
    assert_wholenumber(shift_start)
    assert_wholenumber(shift_end)
    assert_len(shift_start, 1)
    assert_len(shift_end, 1)
    new("ShiftFragments", fragments=fragments, shift_start=as.integer(shift_start), shift_end=as.integer(shift_end))
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
    inner <- iterate_fragments(x@fragments)
    max_len <- if(is.na(x@max_len)) .Machine$integer.max else x@max_len
    wrapFragments(
        iterate_length_select_cpp(ptr(inner), x@min_len, max_len), 
        inner
    )
})
setMethod("short_description", "SelectLength", function(x) {
    text_min <- if(!is.na(x@min_len)) sprintf("%d", x@min_len) else "0"
    text_max <- if(!is.na(x@max_len)) sprintf("%d", x@max_len) else "Inf"
    c(
        short_description(x@fragments),
        sprintf("Filter to fragment sizes %s bp-%s bp", x@min_len, x@max_len)
    )
})
#' Subset fragments to only include those in a given size range
#' @param fragments Input fragments object
#' @param min_len Minimum bases in fragment (inclusive)
#' @param max_len Maximum bases in fragment (inclusive)
#' @return Fragments object
#' @details Fragment length is calculated as end-start
#' @export
subset_lengths <- function(fragments, min_len=0L, max_len=NA_integer_) {
    assert_wholenumber(min_len)
    if (!is.na(max_len)) assert_wholenumber(max_len)
    assert_len(min_len, 1)
    assert_len(max_len, 1)
    new("SelectLength", fragments=fragments, min_len=as.integer(min_len), max_len=as.integer(max_len))
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
setMethod("chrNames<-", "ChrSelectName", function(x, ..., value) {
    assert_is_character(value)
    assert_len(value, length(x@chr_names))
    x@chr_names <- value
    x
})
setMethod("iterate_fragments", "ChrSelectName", function(x) {
    inner <- iterate_fragments(x@fragments)
    wrapFragments(
        iterate_chr_name_select_cpp(ptr(inner), x@chr_names),
        inner
    )
})
setMethod("short_description", "ChrSelectName", function(x) {
    c(
        short_description(x@fragments),
        sprintf(
            "Select %d chromosomes by name%s", 
            length(x@chr_names),
            pretty_print_vector(x@chr_names, max_len=3, prefix=": ")
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
setMethod("chrNames<-", "ChrSelectIndex", function(x, ..., value) {
    assert_is_character(value)
    assert_len(value, length(x@chr_index_selection))
    new("ChrRename", fragments=x, chr_names=value)
})

setMethod("iterate_fragments", "ChrSelectIndex", function(x) {
    inner <- iterate_fragments(x@fragments)
    wrapFragments(
        iterate_chr_index_select_cpp(ptr(inner), x@chr_index_selection - 1),
        inner
    )
})
setMethod("short_description", "ChrSelectIndex", function(x) {
    c(
        short_description(x@fragments),
        sprintf(
            "Select %d chromosomes by index%s", 
            length(x@chr_index_selection),
            pretty_print_vector(x@chr_index_selection, max_len=3, prefix=": ")
        )
    )
})

#' Select chromosmes for subsetting or translating chromosome IDs in a fragments object
#' @param fragments Input fragments object
#' @param chromosome_selection List of chromosme IDs (numeric), or names (character). 
#'    The output chromosome ID n will be taken
#'    from the input fragments chromosome with ID/name chromosome_selection[n]. 
#' @export
select_chromosomes <- function(fragments, chromosome_selection) {
    assert_is(fragments, "IterableFragments")
    assert_distinct(chromosome_selection)
    assert_is(chromosome_selection, c("character", "numeric"))
    if(is.numeric(chromosome_selection)) {
        assert_greater_than_zero(chromosome_selection)
        new_names <- NA_character_
        if (!is.null(chrNames(fragments))) {
            assert_true(all(chromosome_selection <= length(chrNames(fragments))))
        }
        return(
            new("ChrSelectIndex", fragments=fragments, chr_index_selection=as.integer(chromosome_selection))
        )
    } else {
        return(
            new("ChrSelectName", fragments=fragments, chr_names=chromosome_selection)
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
setMethod("cellNames<-", "CellSelectName", function(x, ..., value) {
    assert_is_character(value)
    assert_len(value, length(x@cell_names))
    x@cell_names <- value
    x
})
setMethod("iterate_fragments", "CellSelectName", function(x) {
    inner <- iterate_fragments(x@fragments)
    wrapFragments(
        iterate_cell_name_select_cpp(ptr(inner), x@cell_names),
        inner
    )
})
setMethod("short_description", "CellSelectName", function(x) {
    c(
        short_description(x@fragments),
        sprintf(
            "Select %d cells by name%s", 
            length(x@cell_names),
            pretty_print_vector(x@cell_names, max_len=3, prefix=": ")
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
setMethod("cellNames<-", "CellSelectIndex", function(x, ..., value) {
    assert_is_character(value)
    assert_len(value, length(x@cell_index_selection))
    new("CellRename", fragments=x, cell_names=value)
})
setMethod("iterate_fragments", "CellSelectIndex", function(x) {
    inner <- iterate_fragments(x@fragments)
    wrapFragments(
        iterate_cell_index_select_cpp(ptr(inner), x@cell_index_selection - 1),
        inner
    )
})
setMethod("short_description", "CellSelectIndex", function(x) {
    c(
        short_description(x@fragments),
        sprintf(
            "Select %d cells by index%s", 
            length(x@cell_index_selection),
            pretty_print_vector(x@cell_index_selection, max_len=3, prefix=": ")
        )
    )
})
#' Select cells for subsetting or translating cells IDs in a fragments object
#' @param fragments Input fragments object
#' @param cell_selection List of chromosme IDs (numeric), or names (character). 
#'    The output chromosome ID n will be taken
#'    from the input fragments chromosome with ID/name cell_selection[n]. 
#' @export
select_cells <- function(fragments, cell_selection) {
    assert_is(fragments, "IterableFragments")
    assert_distinct(cell_selection)
    assert_is(cell_selection, c("character", "numeric"))
    if(is.numeric(cell_selection)) {
        assert_greater_than_zero(cell_selection)
        new_names <- NA_character_
        if (!is.null(cellNames(fragments))) {
            assert_true(all(cell_selection <= length(cellNames(fragments))))
        }
        return(
            new("CellSelectIndex", fragments=fragments, cell_index_selection=as.integer(cell_selection))
        )
    } else {
        return(
            new("CellSelectName", fragments=fragments, cell_names=cell_selection)
        )
    }
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
setMethod("chrNames<-", "ChrRename", function(x, ..., value) {
    assert_is_character(value)
    assert_len(value, length(x@chr_names))
    x@chr_names <- value
    x
})
setMethod("iterate_fragments", "ChrRename", function(x) {
    inner <- iterate_fragments(x@fragments)
    wrapFragments(
        iterate_chr_rename_cpp(ptr(inner), x@chr_names),
        inner
    )
})
setMethod("short_description", "ChrRename", function(x) {
    c(
        short_description(x@fragments),
        sprintf("Rename chromosomes to: %s", pretty_print_vector(x@chr_names, max_len=3))
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
setMethod("cellNames<-", "CellRename", function(x, ..., value) {
    assert_is_character(value)
    assert_len(value, length(x@cell_names))
    x@cell_names <- value
    x
})
setMethod("iterate_fragments", "CellRename", function(x) {
    inner <- iterate_fragments(x@fragments)
    wrapFragments(
        iterate_cell_rename_cpp(ptr(inner), x@cell_names),
        inner
    )
})
setMethod("short_description", "CellRename", function(x) {
    c(
        short_description(x@fragments),
        sprintf("Rename cells to: %s", pretty_print_vector(x@cell_names, max_len=3))
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
setMethod("iterate_fragments", "CellPrefix", function(x) {
    inner <- iterate_fragments(x@fragments)
    wrapFragments(
        iterate_cell_prefix_cpp(ptr(inner), x@prefix),
        inner
    )
})
setMethod("short_description", "CellPrefix", function(x) {
    c(
        short_description(x@fragments),
        sprintf("Prefix cells names with: %s", x@prefix)
    )
})

#' Rename cells by adding a prefix to the names (e.g. sample name)
#' @param fragments Input fragments object.
#' @param prefix String to add as the prefix
#' @return Fragments object with prefixed names
#' @export
prefix_cell_names <- function(fragments, prefix, invert_selection=FALSE, zero_based_coords=TRUE) {
    assert_is(fragments, "IterableFragments")
    assert_is(prefix, "character")
    assert_len(prefix, 1)
    new("CellPrefix", fragments=fragments, prefix=prefix)
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
    inner <- iterate_fragments(x@fragments)
    wrapFragments(
        iterate_region_select_cpp(ptr(inner), x@chr_id, x@start, x@end, x@chr_levels, x@invert_selection),
        inner
    )
})
setMethod("short_description", "RegionSelect", function(x) {
    # Subset strings first to avoid a very slow string concatenation process
    indices <- c(head(seq_along(x@chr_id), 3), tail(seq_along(x@chr_id), 1))
    labels <- paste0(x@chr_levels[1+x@chr_id[indices]], ":", x@start[indices]+1, "-", x@end[indices])

    c(
        short_description(x@fragments),
        sprintf("Subset to fragments %soverlapping %d ranges%s", 
            ifelse(x@invert_selection, "not ", ""),
            length(x@chr_id), 
            pretty_print_vector(labels, prefix=": ", max_len=2)
        )
    )
})

#' Select fragments overlapping (or not overlapping) selected regions
#' @param fragments Input fragments object.
#' @param invert_selection If TRUE, select fragments *not* overlapping selected regions
#'  instead of only fragments overlapping the selected regions.
#' @inheritParams peakMatrix
#' @return Fragments object filtered according to the selected regions
#' @export
selectRegions <- function(fragments, ranges, invert_selection=FALSE, zero_based_coords=TRUE) {
    assert_is(fragments, "IterableFragments")
    assert_is(ranges, c("GRanges", "list", "data.frame"))
    
    assert_is(zero_based_coords, "logical")
    
    if(is.list(ranges)) {
        assert_has_names(ranges, c("chr", "start", "end"))
        chr_levels <- as.character(unique(ranges$chr))
        chr_id <- as.integer(factor(as.character(ranges$chr), chr_levels)) - 1L
        start <- as.integer(ranges$start) - !zero_based_coords
        end <- as.integer(ranges$end)
    } else {
        chr_levels <- as.character(unique(GenomicRanges::seqnames(ranges)))
        chr_id <- as.integer(factor(as.character(GenomicRanges::seqnames(ranges)), chr_levels)) - 1L
        start <- GenomicRanges::start(ranges) - !zero_based_coords
        end <- GenomicRanges::end(ranges)
    }
    # Check to make sure regions are end-sorted
    pair_1 <- seq_len(length(chr_id)-1)
    pair_2 <- pair_1 + 1
    if (any(chr_id[pair_1] == chr_id[pair_2] & end[pair_1] > end[pair_2])) {
        stop("Tile regions must be non-overlapping")
    }

    new("RegionSelect", fragments=fragments, chr_id=chr_id, start=start, end=end, chr_levels=chr_levels, invert_selection=invert_selection)
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
    chrNames(x@fragments_list[[1]])
})
setMethod("cellNames", "MergeFragments", function(x) {
    do.call(c, lapply(x@fragments_list, cellNames))
})

setMethod("chrNames<-", "MergeFragments", function(x, ..., value) {
    for (i in seq_along(x@fragments_list)) {
        chrNames(x@fragments_list[[i]]) <- value
    }
    x
})
setMethod("cellNames<-", "MergeFragments", function(x, ..., value) {
    # TODO: Add
    start_idx <- 1
    for (i in seq_along(x@fragments_list)) {
        end_idx <- start_idx - 1 + length(cellNames(x@fragments_list[[i]])) 
        cellNames(x@fragments_list[[i]]) <- value[start_idx:end_idx]
        start_idx <- end_idx + 1
    }
    x
})

setMethod("iterate_fragments", "MergeFragments", function(x) {
    inner <- lapply(x@fragments_list, iterate_fragments)
    wrapFragments(
        iterate_merge_fragments_cpp(lapply(inner, ptr)),
        new("XPtrList", pointers=inner)
    )
})
setMethod("short_description", "MergeFragments", function(x) {
    # Subset strings first to avoid a very slow string concatenation process
    sprintf("Merge %d fragments objects of classes%s",
        length(x@fragments_list),
        pretty_print_vector(vapply(x@fragments_list, class, character(1)), prefix=": ", max_len=3)
    )
})

# Allow merging fragments using standard concatenation method
setMethod("c", "IterableFragments", function(x, ...) {
    tail <- list(...)
    if (length(tail) == 1 && is.null(tail[[1]]))
        args <- list(x)
    else
        args <- c(list(x), list(...))
    chr_names <- chrNames(x)
    assert_not_null(chr_names)
    fragments_list <- list()
    for (i in seq_along(args)) {
        assert_is(args[[i]], "IterableFragments")
        assert_not_null(chrNames(args[[i]]))
        if(!all(chr_names == chrNames(args[[i]])))
            stop("To concatenate fragments, all chrNames must be identical")
        if (is(args[[i]], "MergeFragments"))
            fragments_list <- c(fragments_list, args[[i]]@fragments_list)
        else
            fragments_list <- c(fragments_list, args[[i]])
    }
    new("MergeFragments", fragments_list=fragments_list)
    # TODO: check that cellNames are unique
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
    fragments_identical_cpp(ptr(i1), ptr(i2))
}

#' Scan through fragments without performing any operations (used for benchmarking)
#' @param fragments Fragments object to scan
#' @return Length 4 vector with fragment count, then sums of chr, starts, and ends
scan_fragments <- function(fragments) {
    assert_is(fragments, "IterableFragments")
    i <- iterate_fragments(fragments)
    scan_fragments_cpp(ptr(i))
}


pretty_print_vector <- function(x, max_len=3, sep=", ", prefix="", empty="") {
    # Print whole vector if possible, or else
    # print first few, ellipsis, then the last one
    # if x is length 0, just return the empty value
    if (length(x) == 0) return(empty)

    if (length(x) <= max_len) {
        res <- paste0(x, collapse=sep)
    } else {
        res <- paste0(x[seq_len(max_len-1)], collapse=sep)
        res <- paste0(res, " ... ", x[length(x)], collapse=sep)
    }
    
    res <- paste0(prefix, res)
    res
}