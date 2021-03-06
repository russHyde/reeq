% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coverage.R
\name{which_k_covered_each_sample}
\alias{which_k_covered_each_sample}
\title{which_k_covered_each_sample}
\usage{
which_k_covered_each_sample(mat, k = 1)
}
\arguments{
\item{mat}{Either a matrix, data.frame, DGEList or
ExpressionSet.}

\item{k}{A single number.}
}
\description{
Returns the row index for those rows in a matrix where all entries in that
row are >= k.
}
