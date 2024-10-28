# Copyright 2024 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

# This script is used to automatically update docs-html branch when updating main branch
build_selective_vignette_changes <- function() {
    vignettes_folder <- "r/vignettes/"
    # Get all changes in the last commit
    changed_files <- system("git diff --name-only HEAD~1", intern = TRUE)
    pkgdown::build_articles_index('r')
    # split to only get files in r/vignettes and only get the Rmd files
    changed_vignettes <- changed_files[grepl(vignettes_folder, changed_files)]
    changed_vignettes <- changed_vignettes[grepl("\\.Rmd$", changed_vignettes)]
    for (vgn in changed_vignettes) {
      pkgdown::build_article(gsub("\\.Rmd$", "", gsub("r/vignettes/", "", vgn)), 'r')
    }
}


pkgdown::build_reference('r')
build_selective_vignette_changes()
pkgdown::build_news('r')
pkgdown:::build_sitemap('r')
pkgdown::build_redirects('r')
pkgdown::build_search('r')
pkgdown::build_home('r')
