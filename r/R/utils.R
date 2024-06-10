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
