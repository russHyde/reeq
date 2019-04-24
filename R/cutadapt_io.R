###############################################################################

# ---- Main `cutadapt` parser

#' Extract the basepair/readpair counts from the 'summary' section in the text
#' from a cutadapt log-file.
#'
#' @param        x             A filename. This should be a cutadapt logfile.
#'   The function dies if this is not an existing file or if more than one
#'   entry is present in \code{x}.
#'
#' @importFrom   readr         read_file
#' @export
#'

import_cutadapt_summary <- function(x) {
  # x is a file path

  # validate
  stopifnot(length(x) == 1)
  stopifnot(file.exists(file.path(x)))

  # read
  log_txt <- readr::read_file(x)

  # parse the cutadapt logfile
  summary_txt <- extract_cutadapt_summary(log_txt)

  # convert the summary-data to a tibble and return
  parse_cutadapt_summary(summary_txt)
}

#' Extract the 'summary' section from the text from a cutadapt log-file.
#'
#' Returns a single string.
#'
#' Not exported.
#'
#' @param        x             The raw text from a cutadapt logfile. This fails
#'   if more than one text entry is present, so the user should not split the
#'   file on newlines or anything like that before calling this function.
#'

extract_cutadapt_summary <- function(x) {
  # `x` is the bare text from a cutadapt logfile

  # return the bare text for the `summary` section of the logfile

  stopifnot(is.character(x) && length(x) == 1)
  # Extract the "Summary" section of the file
  without_header <- gsub("^.*=== Summary ===\n", "", x)
  without_footer <- gsub(
    "\n\n=== First read: Adapter 1 ===.*$", "", without_header
  )

  without_footer
}

#' Takes the 'summary' text from a cutadapt log-file and extracts the number of
#' readpairs / basepairs that were input and which pass/fail various filters.
#'
#' Not exported
#'
#' @param        x             Text from the summary section of a cutadapt log.
#'
#' @importFrom   tibble        tibble
#'

parse_cutadapt_summary <- function(x) {
  # `x` should be newline-stripped / bareline stripped /
  #   have the "=== SUMMARY ===" line removed prior to calling this
  if (missing(x)) {
    stop("`x` should be defined in `parse_cutadapt_summary`")
  }

  fieldnames <- tibble::tibble(
    expected = c(
      "Total read pairs processed", "Read 1 with adapter",
      "Read 2 with adapter",
      "Pairs that were too short", "Pairs written (passing filters)",
      "Total basepairs processed", "Read 1", "Read 2",
      "Total written (filtered)", "Read 1", "Read 2"
    ),
    output = c(
      "rp_input", "r1_with_adapter",
      "r2_with_adapter",
      "rp_too_short", "rp_output",
      "bp_input", "bp_input_r1", "bp_input_r2",
      "bp_output", "bp_output_r1", "bp_output_r2"
    )
  )

  parse_numeric_fields(x, fieldnames)
}

###############################################################################

# ---- Helpers

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
#' @importFrom   dplyr        mutate   mutate_   select_
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
    parse_colon_separated_lines()

  if (!isTRUE(all.equal(lines_as_df$field, fieldnames$expected))) {
    message(setdiff(lines_as_df$field, fieldnames$expected))
    message(setdiff(fieldnames$expected, lines_as_df$field))
    stop(
      "All numeric fields from the text should be in the fieldnames dataframe"
    )
  }

  lines_as_df %>%
    dplyr::mutate(field = fieldnames$output) %>%
    dplyr::mutate_(val = ~readr::parse_number(value)) %>%
    dplyr::select_(.dots = c("field", "val")) %>%
    tidyr::spread_(key_col = "field", value_col = "val")
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
