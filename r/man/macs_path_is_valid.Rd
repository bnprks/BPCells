% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/atac_utils.R
\name{macs_path_is_valid}
\alias{macs_path_is_valid}
\title{Test if MACS executable is valid.
If macs_executable is NULL, this function will try to auto-detect MACS from PATH, with preference for MACS3 over MACS2.
If macs_executable is provided, this function will check if MACS can be called.}
\usage{
macs_path_is_valid(macs_executable)
}
\arguments{
\item{macs_executable}{(string) Path to either MACS2/3 executable. Default (\code{NULL}) will autodetect from PATH.}
}
\value{
MACS executable path.
}
\description{
Test if MACS executable is valid.
If macs_executable is NULL, this function will try to auto-detect MACS from PATH, with preference for MACS3 over MACS2.
If macs_executable is provided, this function will check if MACS can be called.
}
\keyword{internal}
