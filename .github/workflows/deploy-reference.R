build_selective_vignette_changes <- function() {
    vignettes_folder <- "r/vignettes/"
    # Get all changes in the last commit
    changed_files <- system("git diff --name-only HEAD~1", intern = TRUE)
    # split to only get files in r/vignettes and only get the Rmd files
    changed_vignettes <- changed_files[grepl(vignettes_folder, changed_files)]
    changed_vignettes <- changed_vignettes[grepl("\\.Rmd$", changed_vignettes)]
    for (vgn in changed_vignettes) {
      pkgdown::build_article(gsub("\\.Rmd$", "", gsub("r/vignettes/", "", vgn)), 'r')
    }
}
args = commandArgs(TRUE)
# First argument is whether to rebuild the entire site.
rebuild <- as.logical(args[1])
if (rebuild) {
  pkgdown::build_site_github_pages('r')
} else {
  pkgdown::build_reference('r')
  build_selective_vignette_changes()
  pkgdown::build_news('r')
  pkgdown::build_sitemap('r')
  pkgdown::build_redirects('r')
  pkgdown::build_search('r')
  pkgdown::build_home('r')
}