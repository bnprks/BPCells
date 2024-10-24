# get_vignette_changes <- function() {
#     vignettes_folder <- "r/vignettes/"
#     # Get all changes in the last commit
#     changed_files <- system("git diff --name-only HEAD~1", intern = TRUE)
#     # split to only get files in r/vignettes
#     changed_vignettes <- changed_files[grepl(vignettes_folder, changed_files)]
#     # get all vignettes in the vignettes folder 
#     all_vignettes <- paste0("r/vignettes/", list.files('r/vignettes', recursive = TRUE, pattern = "\\.Rmd$"))
#     # get the ones that are not changed and change the extension so they're not rendered
#     unchanged_vignettes <- setdiff(all_vignettes, changed_vignettes)
#     unchanged_vignettes <- gsub("\\.Rmd$", ".Rmd_", unchanged_vignettes)
#     pkgdown::build_articles('r', lazy = TRUE)
#     # build vignettes
# }
args = commandArgs(TRUE)
# First argument is whether to rebuild the entire site.
rebuild <- as.logical(args[1])
if (rebuild) {
  pkgdown::build_site_github_pages('r')
} else {
  pkgdown::build_reference('r')
  pkgdown::build_articles('r', lazy = TRUE)
  pkgdown::build_news('r')
  pkgdown::build_home('r')
}