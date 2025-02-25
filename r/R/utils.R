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

# Add current timestamp to a character string.
add_timestamp <- function(msg) {
  return(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S "), msg))
}

# Function which prints a message using shell echo.
# Useful for printing messages from inside mclapply when running in Rstudio.
log_progress <- function(msg, add_timestamp = TRUE){
  if (add_timestamp) {
    msg <- add_timestamp(msg)
  }
  if (.Platform$GUI == "RStudio") {
    system(sprintf('echo "%s"', paste0(msg, collapse="")))
  } else {
    message(msg)
  }
}

#' Helper to create partial functions
#'
#' Automatically creates a partial application of the caller
#' function including all non-missing arguments. 
#'
#' @return A `bpcells_partial` object (a function with some extra attributes)
#' @keywords internal
create_partial <- function() {
  env <- rlang::caller_env()
  fn_sym <- rlang::caller_call()[[1]]
  fn <- rlang::caller_fn()

  args <- list()
  for (n in names(formals(fn))) {
    if (rlang::is_missing(env[[n]])) next
    args[[n]] <- env[[n]]
  }

  ret <- do.call(partial_apply, c(fn, args))
  attr(ret, "body")[[1]] <- fn_sym
  return(ret)
}


#' Create partial function calls
#' 
#' Specify some but not all arguments to a function.
#'
#' @param f A function
#' @param ... Named arguments to `f`
#' @param .overwrite (bool) If `f` is already an output from
#'   `partial_apply()`, whether parameter re-definitions should
#'   be ignored or overwrite the existing definitions
#' @param .missing_args_error (bool) If `TRUE`, passing in arguments
#'  that are not in the function's signature will raise an error, otherwise
#'  they will be ignored
#' @return A `bpcells_partial` object (a function with some extra attributes)
#' @keywords internal
partial_apply <- function(f, ..., .overwrite = TRUE, .missing_args_error = TRUE) {
  args <- rlang::list2(...)

  if (is(f, "bpcells_partial")) {
    prev_args <- attr(f, "args")
    for (a in names(prev_args)) {
      if (!(.overwrite && a %in% names(args))) {
        args[[a]] <- prev_args[[a]]
      }
    }
    f <- attr(f, "fn")
    function_name <- attr(f, "body")[[1]]
  } else {
    function_name <- rlang::sym(rlang::caller_arg(f))
  }
  # See which arguments do not exist in f
  missing_args <- which(!names(args) %in% names(formals(f)))
  if (length(missing_args) > 0) {
    if (.missing_args_error) {
      stop(sprintf("Arguments %s are not in the function signature", paste0(names(args)[missing_args], collapse=", ")))
    } else {
      args <- args[-missing_args]}
  }
  partial_fn <- do.call(purrr::partial, c(f, args))
  attr(partial_fn, "body")[[1]] <- function_name
  structure(
    partial_fn,
    class = c("bpcells_partial", "purrr_function_partial", "function"),
    args = args,
    fn = f
  )
}
