

#' IterableFragments methods
#'
#' Methods for IterableFragments objects
#'
#' @name IterableFragments-methods
#' @rdname IterableFragments-methods
NULL

setClass("IterableFragments", 
    slots = c(
        cell_names = "character",
        chr_names = "character"
    ),
    prototype = list(
        cell_names = NA_character_,
        chr_names = NA_character_
    )
)

# Get an external pointer to a C++ FragmentsLoader object. This is only called
# when an operation is actively ready to run
setGeneric("iterate_fragments", function(x) standardGeneric("iterate_fragments"))
setMethod("iterate_fragments", "externalptr", function(x) { return (x)} )

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
setGeneric("cellNames", function(x) standardGeneric("cellNames"))
setMethod("cellNames", "IterableFragments", function(x) {
    if(length(x@cell_names) == 1 && all(is.na(x@cell_names))) {
        return(NULL)
    }
    return(x@cell_names)
})

#' Set cell names
#' @param x an IterableFragments object
#' @param value Character vector of new names
#' @details * `cellNames<-` It is only possible to replace names, not add new names.
#' @describeIn IterableFragments-methods Set cell names
setGeneric("cellNames<-", function(x, ..., value) standardGeneric("cellNames<-"))
setMethod("cellNames<-", "IterableFragments", function(x, ..., value) {
    assert_is_character(value)
    if (is.null(cellNames(x))) {
        stop("Assigning new cellNames is not allowed, only renaming")
    }
    assert_len(value, length(cellNames(x)))
    x@cell_names <- value
    return(x)
})

#' Get chromosome names
#' @param x an IterableFragments object
#' @return * `chrNames()`: Character vector of chromosome names, or NULL if none are known
#' @describeIn IterableFragments-methods Set chromosome names
setGeneric("chrNames", function(x) standardGeneric("chrNames"))
setMethod("chrNames", "IterableFragments", function(x) {
    if(length(x@chr_names) == 1 && all(is.na(x@chr_names))) {
        return(NULL)
    }
    return(x@chr_names)
})

#' Set chromosome names
#' @param x an IterableFragments object
#' @param value Character vector of new names
#' @details * `chrNames<-` It is only possible to replace names, not add new names.
#' @describeIn IterableFragments-methods Set chromosome names
setGeneric("chrNames<-", function(x, ..., value) standardGeneric("chrNames<-"))
setMethod("chrNames<-", "IterableFragments", function(x, ..., value) {
    assert_is_character(value)
    if (is.null(chrNames(x))) {
        stop("Assigning new chrNames is not allowed, only renaming")
    }
    assert_len(value, length(chrNames(x)))
    x@chr_names <- value
    return(x)
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

#' Read a 10x fragments file
#' @details Note: no disk operations will take place until the fragments are used in a function
#' @param path Path of 10x fragments file
#' @param comment Skip lines at start of file starting with comment
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
setMethod("iterate_fragments", "FragmentsTsv", function(x) iterate_10x_fragments_cpp(normalizePath(x@path), x@comment))
setMethod("short_description", "FragmentsTsv", function(x) {
    sprintf("Load 10x fragments file from %s", x@path)
})

#' Write to a 10x fragments file
#' @param fragments Input fragments object
#' @param path Output path
#' @param append_5th_column Whether to include 5th column of all 0 for compatibility
#'        with 10x fragment file outputs (defaults to 4 columns chr,start,end,cell)
#' @export
write_fragments_10x <- function(fragments, path, append_5th_column=FALSE) {
    assert_is_file(path, must_exist=FALSE, extension=c(".tsv", ".tsv.gz"))

    write_10x_fragments_cpp(
        normalizePath(path, mustWork=FALSE), 
        iterate_fragments(fragments),
        append_5th_column
    )

    open_fragments_10x(path, comment="")
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


setMethod("iterate_fragments", "PackedMemFragments", function(x) {
    iterate_packed_fragments_cpp(x)
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
    if (compress) {
        res <- write_packed_fragments_cpp(iterate_fragments(fragments))
        do.call(new, c("PackedMemFragments", res))
    } else {
        res <- write_unpacked_fragments_cpp(iterate_fragments(fragments))
        do.call(new, c("UnpackedMemFragments", res))
    }
}

setClass("FragmentsDir",
    contains = "IterableFragments",
    slots = c(
        dir = "character",
        compressed = "logical",
        buffer_size = "integer"
    ),
    prototype = list(
        dir = character(0),
        compressed = TRUE,
        buffer_size = 1024L
    )
)

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
    if (compress)
        write_packed_fragments_file_cpp(iterate_fragments(fragments), dir, buffer_size)
    else
        write_unpacked_fragments_file_cpp(iterate_fragments(fragments), dir, buffer_size)
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

setMethod("iterate_fragments", "FragmentsDir", function(x) {
    if (x@compressed)
        iterate_packed_fragments_file_cpp(x@dir, x@buffer_size)
    else
        iterate_unpacked_fragments_file_cpp(x@dir, x@buffer_size)
})
setMethod("short_description", "FragmentsDir", function(x) {
    sprintf("Read %s fragments from directory %s", 
        if(x@compressed) "compressed" else "uncompressed",
        x@dir
    )
})

setClass("FragmentsHDF5",
    contains = "IterableFragments",
    slots = c(
        path = "character",
        group = "character",
        compressed = "logical",
        buffer_size = "integer"
    ),
    prototype = list(
        path = character(0),
        group = "",
        compressed=TRUE,
        buffer_size = 8192L
    )
)


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
    if (compress)
        write_packed_fragments_hdf5_cpp(iterate_fragments(fragments), path, group, buffer_size, chunk_size)
    else
        write_unpacked_fragments_hdf5_cpp(iterate_fragments(fragments), path, group, buffer_size, chunk_size)
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

setMethod("iterate_fragments", "FragmentsHDF5", function(x) {
    if (x@compressed)
        iterate_packed_fragments_hdf5_cpp(x@path, x@group, x@buffer_size)
    else
        iterate_unpacked_fragments_hdf5_cpp(x@path, x@group, x@buffer_size)
})
setMethod("short_description", "FragmentsHDF5", function(x) {
    sprintf("Read %s fragments from %s, group %s", 
        if(x@compressed) "compressed" else "uncompressed",
        x@path, 
        x@group
    )
})


#' Build a Fragments object from an R data frame or GRanges object
#' @param x An input GRanges, list or data frame. Lists and dataframes 
#'    must have chr, start, end, cell_id. GRanges must have a metadata column
#'    for cell_id
#' @param zero_based_coords Whether to convert the ranges from a 1-based end-inclusive
#'    coordinate system to a 0-based end-exclusive coordinate system. Defaults to true
#'    for GRanges and false for other formats
#'    (see http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/)
#' @return RawFragments object representing the given fragments
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

    chr_ptr <- rep(cumsum(rle(as.integer(x$chr))$lengths), each=2)
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
            seqnames=S4Vectors::Rle(values=chrNames(frags)[order(chr_order)], lengths=chr_counts),
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
    new("ShiftFragments", fragments=fragments, shift_start=as.integer(shift_start), shift_end=as.integer(shift_end),
        cell_names=fragments@cell_names, chr_names=fragments@chr_names
    )
}
setMethod("iterate_fragments", "ShiftFragments", function(x) {
    iterate_shift_cpp(iterate_fragments(x@fragments), x@shift_start, x@shift_end)
})
setMethod("short_description", "ShiftFragments", function(x) {
    c(
        short_description(x@fragments),
        sprintf("Shift start %+dbp, end %+dbp", x@shift_start, x@shift_end)
    )
})

# Select fragments by chromosome
setClass("ChrSelectName",
    contains = "IterableFragments",
    slots = c(
        fragments = "IterableFragments"
    ),
    prototype = list(
        fragments = NULL
    )
)
setMethod("iterate_fragments", "ChrSelectName", function(x) {
    iterate_chr_name_select_cpp(iterate_fragments(x@fragments), x@chr_names)
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
setMethod("iterate_fragments", "ChrSelectIndex", function(x) {
    iterate_chr_index_select_cpp(iterate_fragments(x@fragments), x@chr_index_selection - 1)
})
setMethod("short_description", "ChrSelectIndex", function(x) {
    c(
        short_description(x@fragments),
        sprintf(
            "Select %d chromosomes by index%s", 
            length(x@chr_names),
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
        if (!all(is.na(fragments@chr_names))) {
            assert_true(all(chromosome_selection <= length(fragments@chr_names)))
            new_names <- fragments@chr_names[chromosome_selection]
        }
        return(
            new("ChrSelectIndex", fragments=fragments, chr_index_selection=as.integer(chromosome_selection),
                cell_names=fragments@cell_names, chr_names=new_names)
        )
    } else {
        return(
            new("ChrSelectName", fragments=fragments, chr_names=chromosome_selection,
                cell_names=fragments@cell_names)
        )
    }
}



# Select fragments by cell
setClass("CellSelectName",
    contains = "IterableFragments",
    slots = c(
        fragments = "IterableFragments"
    ),
    prototype = list(
        fragments = NULL
    )
)
setMethod("iterate_fragments", "CellSelectName", function(x) {
    iterate_cell_name_select_cpp(iterate_fragments(x@fragments), x@cell_names)
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
setMethod("iterate_fragments", "CellSelectIndex", function(x) {
    iterate_cell_index_select_cpp(iterate_fragments(x@fragments), x@cell_index_selection - 1)
})
setMethod("short_description", "CellSelectIndex", function(x) {
    c(
        short_description(x@fragments),
        sprintf(
            "Select %d cells by index%s", 
            length(x@chr_names),
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
        if (!all(is.na(fragments@cell_names))) {
            assert_true(all(cell_selection <= length(fragments@cell_names)))
            new_names <- fragments@cell_names[cell_selection]
        }
        return(
            new("CellSelectIndex", fragments=fragments, cell_index_selection=as.integer(cell_selection),
                chr_names=fragments@chr_names, cell_names=new_names)
        )
    } else {
        return(
            new("CellSelectName", fragments=fragments, cell_names=cell_selection,
                chr_names=fragments@chr_names)
        )
    }
}

#' Check if two fragments objects are identical
#' @param fragments1 First IterableFragments to compare
#' @param fragments2 Second IterableFragments to compare
#' @return boolean for whether the fragments objects are identical
fragments_identical <- function(fragments1, fragments2) {
    assert_is(fragments1, "IterableFragments")
    assert_is(fragments2, "IterableFragments")
    fragments_identical_cpp(iterate_fragments(fragments1), iterate_fragments(fragments2))
}

#' Scan through fragments without performing any operations (used for benchmarking)
#' @param fragments Fragments object to scan
#' @return Length 4 vector with fragment count, then sums of chr, starts, and ends
#' @export
scan_fragments <- function(fragments) {
    assert_is(fragments, "IterableFragments")
    scan_fragments_cpp(iterate_fragments(fragments))
}

#' Count fragments by nucleosomal size
#' @param fragments Fragments object
#' @param nucleosome_width Integer cutoff to use as nucleosome width
#' @return List with names subNucleosomal, monoNucleosomal, multiNucleosomal containing the 
#'         count vectors of fragments in each class per cell.
#' @export
nucleosome_counts <- function(fragments, nucleosome_width=147) {
    assert_is(fragments, "IterableFragments")
    assert_wholenumber(nucleosome_width)
    assert_len(nucleosome_width, 1)

    iter <- iterate_fragments(fragments)
    nucleosome_counts_cpp(iter, nucleosome_width)
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