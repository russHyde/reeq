% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{replace_with}
\alias{replace_with}
\title{replace_with}
\usage{
replace_with(x, search_list, return_list, strict = FALSE)
}
\arguments{
\item{x}{A vector.}

\item{search_list}{A vector.}

\item{return_list}{A vector. This should be the same length as
\code{search_list}.}

\item{strict}{Boolean. Should every element of `x` be present
in \code{search_list}?}
}
\value{
A vector of the same length as \code{x}, containing some
  values from \code{x} and some values from \code{return_list} (so the class
  of the vector may be different from that for `x`).
}
\description{
For each element in \code{x}, if it is present in the `search_list`, replace
it with the corresponding value in the `return_list`. Any element that is
not in the `search_list` is returned unaltered.
}
