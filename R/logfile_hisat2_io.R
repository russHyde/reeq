###############################################################################

# Functions for parsing alignment statistics from `hisat2` logfiles

###############################################################################

#' Extract the alignment summary statistics from a hisat2 logfile
#'
#' @param        x             A filepath. This should be a hisat2 logfile. The
#'   function dies if this is not an existing file of if more than one entry is
#'   present in \code{x}.
#'
#' @importFrom   readr         read_file
#'
#' @export
#'

import_hisat2_summary <- function(x) {
  # validate the input is a single filepath
  stopifnot(length(x) == 1)
  stopifnot(file.exists(x))

  log_txt <- readr::read_file(x)

  parse_hisat2_summary(log_txt)
}

#' parse_hisat2_summary
#'
#' @param        x             The text from a hisat2 summary file. It is
#'   assumed that the `--new-summary` flag was used while making this file.
#'   There should be only a single entry in `x`.
#'
#' @importFrom   dplyr         mutate
#' @importFrom   rlang         .data

parse_hisat2_summary <- function(x) {
  # Assumes the input text is from a --new-summary from hisat2
  # - Work out how to distinguish --new-summary from old summary files
  stopifnot(is.character(x))
  stopifnot(length(x) == 1)

  parse_numeric_fields(x) %>%
    spread_and_rename_hisat2_fieldnames() %>%
    dplyr::mutate(align_rate = .data[["align_rate"]] / 100)
}

#' spread_and_rename_hisat2_fieldnames
#'
#' @param        x             A datafame with two columns "field" and "value".
#'   The numeric fields from a hisat2 logfile should have been extracted into
#'   this dataframe.
#'
#' @include      logfile_helpers.R
#'

spread_and_rename_hisat2_fieldnames <- function(x) {
  define_hisat2_summary_renaming <- function() {
    tibble::tribble(
      ~expected, ~output,
      "Total pairs", "rp_input",
      "Aligned concordantly or discordantly 0 time", "rp_zero",
      "Aligned concordantly 1 time", "rp_concordant_once",
      "Aligned concordantly >1 times", "rp_concordant_multi",
      "Aligned discordantly 1 time", "rp_discordant_once",
      "Total unpaired reads", "unpaired",
      "Aligned 0 time", "unpaired_unaligned",
      "Aligned 1 time", "unpaired_once",
      "Aligned >1 times", "unpaired_multi",
      "Overall alignment rate", "align_rate"
    )
  }

  # --

  fieldnames <- define_hisat2_summary_renaming()

  spread_and_rename_numeric_fields(x, fieldnames)
}
