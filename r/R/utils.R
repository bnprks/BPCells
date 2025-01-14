# Copyright 2021 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

# Helper-function to get random seed. Adapted from withr:::get_seed
get_seed <- function() {
  if (exists(".Random.seed", globalenv(), mode = "integer", inherits = FALSE)) {
    return(get(".Random.seed", globalenv(), mode = "integer", inherits = FALSE))
  } else {
    return(NULL)
  }
}
# Helper-function to set random seed. Adapted from withr:::restore_seed / withr:::rm_seed
restore_seed <- function(seed) {
  if (is.null(seed)) {
    rm(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    assign(".Random.seed", seed, envir = globalenv(), inherits = FALSE)
  }
}

# Helper-function for pkgdown documentation about genomic ranges inputs
document_granges <- function(
  intro_noun="Genomic regions", 
  position="`chr`, `start`, `end`: genomic position",
  strand=NULL,
  extras=NULL
) {
  if (!is.null(strand) && strand == "default") {
    strand <- "`strand`: +/- or TRUE/FALSE for positive or negative strand"
  }
  if (!is.null(extras)) {
    extras <- sprintf("`%s`: %s", names(extras), as.character(extras))
  }
  bullets <- paste0(sprintf("  - %s", c(position, strand, extras)), collapse="\n")
  sprintf(paste0(
    "%s given as GRanges, data.frame, or list. ",
    "See `help(\"genomic-ranges-like\")` for details on format and coordinate systems. ",
    "Required attributes:\n\n",
    "%s"
  ), intro_noun, bullets)
}

# Function which prints a message using shell echo.
# Useful for printing messages from inside mclapply when running in Rstudio.
log_progress <- function(msg, add_timestamp = TRUE){
  if (add_timestamp) {
    msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S "), msg)
  }
  if (.Platform$GUI == "RStudio") {
    system(sprintf('echo "%s"', paste0(msg, collapse="")))
  } else {
    message(msg)
  }
}

# Helper function to create partial explicit functions
# This builds upon purrr::partial by allowing for nested partial calls, where each partial call
# only does partial application of the arguments that were explicitly provided.
partial_explicit <- function(fn, ...) {
  args <- rlang::enquos(...)
  evaluated_args <- purrr::map(args, rlang::eval_tidy)
  # Fetch the default arguments from the function definition
  default_args <- formals(fn)
  # Keep only explicitly provided arguments that were evaluated
  # where the values are different from the default arguments
  explicitly_passed_args <- evaluated_args[names(evaluated_args) %in% names(default_args) & 
                                            !purrr::map2_lgl(evaluated_args, default_args[names(evaluated_args)], identical)]
  # Return a partially applied version of the function using evaluated arguments
  return(purrr::partial(fn, !!!explicitly_passed_args))
}