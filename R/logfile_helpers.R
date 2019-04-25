###############################################################################

# Helper functions for parsing alignment-tool logfiles

###############################################################################

#' parse_numeric_fields
#'
#' @param        x             A single character string. This should be
#'   newline-separated. Key-Value pairs are assumed to be colon-separated
#'   and any line that is, is converted into a key-value pair.
#' @param        fieldnames    A dataframe containing two columns: expected and
#'   output. The \code{expected} gives the fieldnames that are expected to be
#'   present in the text (\code{x}); \code{output} gives the names that these
#'   fields should be converted to in the output dataframe
#'
#' @return       A dataframe. Contains columns field and val (? others) where
#'   the entries of \code{val} are the values found within the text \code{x}
#'   for the fields in \code{field}.
#'
#' @importFrom   tidyr        spread_
#' @importFrom   readr        parse_number
#' @importFrom   stringr      str_subset   str_replace
#' @importFrom   dplyr        mutate_   select_
#'

parse_numeric_fields <- function(x, fieldnames) {
  # `x` should be bare text, not split
  stopifnot(is.character(x) && length(x) == 1)
  stopifnot(all.equal(colnames(fieldnames), c("expected", "output")))

  extract_numeric_lines <- function(x) {
    # split some text into distinct lines, and keep only those lines with a
    # numeric statistic on the RHS of a colon
    x %>%
      strsplit("\n") %>%
      unlist() %>%
      stringr::str_subset(":.*[[:digit:]].*")
  }

  reformat_numeric_lines <- function(x) {
    # drop trailing
    # - parenthesised, percent eg, "some_stat : 123 (0.5%)"
    # - `bp` basepair indicators eg, "Remaining read pairs : 123 bp"
    # - percent signs eg, "percentage_statistic : 98.1%"
    x %>%
      stringr::str_replace("\\([[:graph:]]*%\\)$", "") %>%
      stringr::str_replace("bp[[:blank:]]*$", "") %>%
      stringr::str_replace("%[[:blank:]]*$", "")
  }

  lines_as_df <- x %>%
    extract_numeric_lines() %>%
    reformat_numeric_lines() %>%
    # convert into key-value (string -> string) pairs
    parse_colon_separated_lines() %>%
    # convert into key-value (string -> number) pairs
    dplyr::mutate_(
      value = ~readr::parse_number(value)
    )

  # TODO: return here
  # - logfiles from different tools will require different methods for
  # disambiguating the original fieldnames and for reformatting the
  # disambiguated fieldnames

  lines_as_df
}


#' parse_colon_separated_lines
#'
#' Split each entry in a vector on the first colon. Strips flanking whitespace
#' from what remains. Returns the left-hand-side in the column `field` and the
#' right-hand-side in a column `value`.
#'
#' Not exported
#'
#' @param        x            Vector of colon-separated values.
#'
#' @return       `tibble` with two columns: 'field' and 'value'. Any
#'   leading/trailing whitespace is trimmed off both the field and value.
#'
#' @importFrom   dplyr         mutate_all
#' @importFrom   magrittr      %>%   set_colnames
#' @importFrom   stringr       str_split_fixed
#' @importFrom   tibble        as_tibble
#'

parse_colon_separated_lines <- function(x) {
  if (missing(x)) {
    stop(
      "character vector `x` should be defined in parse_colon_separated_lines"
    )
  }

  # die if any of the input vector lacks a colon
  stopifnot(all(grepl(":", x)))

  # Split on the first colon
  # Join the values into a two-column dataframe
  # Strip all leading or trailing whitespace:
  x %>%
    stringr::str_split_fixed(":", n = 2) %>%
    magrittr::set_colnames(c("field", "value")) %>%
    tibble::as_tibble() %>%
    dplyr::mutate_all(trimws)
}
