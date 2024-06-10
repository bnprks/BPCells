#!/usr/bin/env bash
set -eu

Rscript -e "pkgdown::build_site()"