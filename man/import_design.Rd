% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_model.R
\name{import_design}
\alias{import_design}
\title{Obtain a design matrix from a file.}
\usage{
import_design(filepath, dge)
}
\arguments{
\item{filepath}{A tab-separated file that defines the
experimental design for the current experiment. This function fails if the
file does not exist.}

\item{dge}{A DGEList containing the expression data for the
current experiment. This is optional. But providing a DGEList helps sanity
check the design matrix - ensuring the samples in the design are all
present in the DGEList.}
}
\description{
Ensures that the samples mentioned in the rows of the design matrix are all
present in an expression dataset (if provided).
}
\details{
But note, the user should ensure that the DGEList is subset / reordered if
the samples in the DGEList is a superset (or unmatching order) or those in
the returned design matrix.

Here, we assume that the first column in the `.tsv` contains the sample
identifers, as used throughout the _R_ workflow. Each entry in that column
should be a single string that can directly be used as a column name.
All other columns should be numeric, and will be converted to a matrix.
}
