% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logfile_helpers.R
\name{parse_colon_separated_lines}
\alias{parse_colon_separated_lines}
\title{parse_colon_separated_lines}
\usage{
parse_colon_separated_lines(x)
}
\arguments{
\item{x}{Vector of colon-separated values.}
}
\value{
`tibble` with two columns: 'field' and 'value'. Any
  leading/trailing whitespace is trimmed off both the field and value.
}
\description{
Split each entry in a vector on the first colon. Strips flanking whitespace
from what remains. Returns the left-hand-side in the column `field` and the
right-hand-side in a column `value`.
}
\details{
Not exported
}
