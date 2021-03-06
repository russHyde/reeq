% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logfile_helpers.R
\name{spread_and_rename_numeric_fields}
\alias{spread_and_rename_numeric_fields}
\title{spread_and_rename_numeric_fields}
\usage{
spread_and_rename_numeric_fields(x, fieldnames)
}
\arguments{
\item{x}{A dataframe containing two columns "field" and
"value". The elements in "field" will be renamed according to the mapping
in `fieldnames`}

\item{fieldnames}{A dataframe containing two columns: "expected"
and "output".}
}
\value{
A dataframe. The columns ordering corresponds to the "field"
  column in `x` (one column for each row in `x` and the columns are named
  according to the reformatting defined by `fieldnames`).
}
\description{
spread_and_rename_numeric_fields
}
