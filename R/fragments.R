

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
#' @return 10x fragments file object
#' @export
open_10x_fragments <- function(path, comment="#") {
    assert_is_file(path, extension=c(".tsv", ".tsv.gz"))
    assert_is_character(comment)
    assert_len(comment, 1)
    
    new("FragmentsTsv", path=path, comment=comment)
}
setMethod("iterate_fragments", "FragmentsTsv", function(x) load_10x_fragments_cpp(normalizePath(x@path), x@comment))
setMethod("short_description", "FragmentsTsv", function(x) {
    sprintf("Load 10x fragments file from %s", x@path)
})

#' Write to a 10x fragments file
#' @param fragments Input fragments object
#' @param path Output path
#' @param append_5th_column Whether to include 5th column of all 0 for compatibility
#'        with 10x fragment file outputs (defaults to 4 columns chr,start,end,cell)
#' @export
write_10x_fragments <- function(fragments, path, append_5th_column=FALSE) {
    assert_is_file(path, must_exist=FALSE, extension=c(".tsv", ".tsv.gz"))

    write_10x_fragments_cpp(
        normalizePath(path, mustWork=FALSE), 
        iterate_fragments(fragments),
        append_5th_column
    )
}

# Read/write in-memory PackedFragments objects
setClass("PackedFragments",
    contains = "IterableFragments",
    slots = c(
        fragments = "list"
    ),
    prototype = list(
        fragments = NULL
    )
)
#' Make a packed_fragments object in memory.
#' @details Memory usage should be about half of the size of a gzip-compressed 10x fragments
#' file holding the same data.
#' @param fragments Input fragments object
#' @return Packed fragments object
#' @export
write_packed_fragments <- function(fragments){
    assert_is(fragments, "IterableFragments")
    packed <- write_packed_fragments_cpp(iterate_fragments(fragments))
    new("PackedFragments", cell_names=packed$cell_names, chr_names=packed$chr_names, fragments=packed$packed_frags)
}
setMethod("iterate_fragments", "PackedFragments", function(x) load_packed_fragments_cpp(x))
setMethod("short_description", "PackedFragments", function(x) {
    "Read packed fragments from memory"
})

# Read/write in-memory RawFragments objects
setClass("RawFragments",
    contains = "IterableFragments",
    slots = c(
        fragments = "list"
    ),
    prototype = list(
        fragments = NULL
    )
)
#' Make a RawFragments object in memory. 
#' @param fragments Input fragments object
#' @return Raw fragments object
#' @export
write_raw_fragments <- function(fragments) {
    assert_is(fragments, "IterableFragments")
    res <- write_unpacked_fragments_cpp(iterate_fragments(fragments))
    
    new("RawFragments", cell_names=res$cell_names, chr_names=res$chr_names, fragments=res$fragments)
}
setMethod("iterate_fragments", "RawFragments", function(x) iterate_unpacked_fragments_cpp(x))
setMethod("short_description", "RawFragments", function(x) {
    "Read raw fragments from memory"
})
#' Build a RawFragments object from an R data frame or GRanges
#' @param x An input GRanges, list or data frame. Lists and dataframes 
#'    must have chr, start, end, cell_id. GRanges must have a metadata column
#'    for cell_id
#' @param convert_to_0_based_coords Whether to convert the ranges from a 1-based
#'    coordinate system to a 0-based coordinate system. 
#'    (see http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/)
#' @return RawFragments object representing the given fragments
RawFragments <- function(x, convert_to_0_based_coords=TRUE) {
    assert_is(x, c("list", "data.frame", "GRanges"))
    assert_is(convert_to_0_based_coords, "logical")
    assert_not_null(x$cell_id)
    if(!is(x, "GRanges")) {
        assert_has_names(x, c("chr", "start", "end", "cell_id"))
        x <- GenomicRanges::GRanges(
            seqnames = as.factor(x$chr),
            ranges = IRanges::IRanges(start = x$start, end=x$end),
            cell_id = as.factor(x$cell_id)
        )
    }
    x <- sort(x)
    x$cell_id <- as.factor(x$cell_id)
    xl <- split(x, GenomicRanges::seqnames(x))
    fragments <- lapply(xl, function(granges) {
        list(
            start = GenomicRanges::start(granges) - 1,
            end = GenomicRanges::end(granges),
            cell_id = as.integer(granges$cell_id) - 1
        )
    })
    names(fragments) <- NULL
    new("RawFragments", cell_names=levels(x$cell_id), chr_names=levels(GenomicRanges::seqnames(x)), fragments=fragments)
}


##### RAW FRAGMENTS V2 TEST
setClass("MemFragments",
    contains = "IterableFragments",
    slots = c(
        fragments = "list",
        compressed = "logical",
        cell_names = "character",
        chr_names = "character"
    ),
    prototype = list(
        fragments = NULL,
        compressed = TRUE
    )
)
#' Make a RawFragments object in memory. 
#' @param fragments Input fragments object
#' @return MemFragments object
#' @export
write_fragments_memory <- function(fragments, compressed=TRUE) {
    assert_is(fragments, "IterableFragments")
    assert_is(compressed, "logical")
    if (compressed)
        res <- write_packed_fragments2_cpp(iterate_fragments(fragments))
    else
        res <- write_unpacked_fragments2_cpp(iterate_fragments(fragments))
    
    new("MemFragments", cell_names=res$cell_names, chr_names=res$chr_names, 
                        fragments=res$fragments, compressed=compressed)
}

setMethod("iterate_fragments", "MemFragments", function(x) {
    if (x@compressed) iterate_packed_fragments2_cpp(x)
    else iterate_unpacked_fragments2_cpp(x)
})
setMethod("short_description", "MemFragments", function(x) {
    sprintf("Read %s fragments from memory", if(x@compressed) "compressed" else "uncompressed")
})
#' Build a RawFragments object from an R data frame or GRanges
#' @param x An input GRanges, list or data frame. Lists and dataframes 
#'    must have chr, start, end, cell_id. GRanges must have a metadata column
#'    for cell_id
#' @param convert_to_0_based_coords Whether to convert the ranges from a 1-based
#'    coordinate system to a 0-based coordinate system. 
#'    (see http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/)
#' @return RawFragments object representing the given fragments
RawFragments2 <- function(x, convert_to_0_based_coords=TRUE) {
    write_raw_fragments2(RawFragments(x, convert_to_0_based_coords))
}

setClass("FragmentsDir",
    contains = "IterableFragments",
    slots = c(
        dir = "character",
        compressed = "logical"
    ),
    prototype = list(
        dir = character(0),
        compressed = TRUE
    )
)

#' @export
write_fragments_dir <- function(fragments, dir, compressed=TRUE) {
    assert_is(fragments, "IterableFragments")
    assert_is(dir, "character")
    assert_is(compressed, "logical")
    if (compressed)
        write_packed_fragments_file_cpp(iterate_fragments(fragments), path.expand(dir))
    else
        write_unpacked_fragments_file_cpp(iterate_fragments(fragments), path.expand(dir))
}

open_fragments_dir <- function(dir, compressed=TRUE) {
    new("FragmentsDir", dir=path.expand(dir), compressed=compressed)
}

setMethod("iterate_fragments", "FragmentsDir", function(x) {
    if (x@compressed)
        iterate_packed_fragments_file_cpp(x@dir)
    else
        iterate_unpacked_fragments_file_cpp(x@dir)
})
setMethod("short_description", "FragmentsDir", function(x) {
    sprintf("Read %s fragments from %s", 
        if(x@compressed) "compressed" else "uncompressed",
        x@dir
    )
})



setClass("FragmentsH5",
    contains = "IterableFragments",
    slots = c(
        path = "character",
        group = "character",
        compressed = "logical"
    ),
    prototype = list(
        path = character(0),
        group = "",
        compressed=TRUE
    )
)

#' @export
write_fragments_h5 <- function(fragments, path, group, chunk_size=250000L, compressed=TRUE) {
    assert_is(fragments, "IterableFragments")
    assert_is(path, "character")
    assert_is(group, "character")
    assert_is(compressed, "logical")
    assert_wholenumber(chunk_size)

    if (compressed)
        write_packed_fragments_hdf5_cpp(iterate_fragments(fragments), path.expand(path), group, as.integer(chunk_size))
    else
        write_unpacked_fragments_hdf5_cpp(iterate_fragments(fragments), path.expand(path), group, as.integer(chunk_size))
}

open_fragments_h5 <- function(path, group="/", compressed=TRUE) {
    if (group == "") group <- "/"
    new("FragmentsH5", path=path, group=group, compressed=compressed)
}

setMethod("iterate_fragments", "FragmentsH5", function(x) {
    if (x@compressed)
        iterate_packed_fragments_hdf5_cpp(path.expand(x@path), x@group)
    else
        iterate_unpacked_fragments_hdf5_cpp(path.expand(x@path), x@group)
})
setMethod("short_description", "FragmentsH5", function(x) {
    sprintf("Read %s fragments from %s, group %s", 
        if(x@compressed) "compressed" else "uncompressed",
        x@path, 
        x@group
    )
})


#### END RAW FRAGMENTS V2 TEST

# Define converters with GRanges if we have the GenomicRanges package available
if (requireNamespace("GenomicRanges", quietly=TRUE)) {
    setAs("IterableFragments", "GRanges", function(from) {
        frags <- write_raw_fragments(from)
        l <- lapply(seq_along(frags@fragments), function(i) {
            if (length(frags@fragments[[i]]$start) == 0) {
                return(GenomicRanges::GRanges())
            }
            GenomicRanges::GRanges(
                seqnames=frags@chr_names[i],
                ranges=IRanges::IRanges(
                    start = frags@fragments[[i]]$start + 1,
                    end = frags@fragments[[i]]$end,
                ),
                cell_id = factor(frags@cell_names[frags@fragments[[i]]$cell_id + 1], frags@cell_names)
            )
        })
        unlist(GenomicRanges::GRangesList(l))
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