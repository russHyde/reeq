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
#' @include      logfile_helpers.R
#'

parse_cutadapt_summary <- function(x) {
  # `x` should be newline-stripped / bareline stripped /
  #   have the "=== SUMMARY ===" line removed prior to calling this
  if (missing(x)) {
    stop("`x` should be defined in `parse_cutadapt_summary`")
  }

  fieldnames <- define_cutadapt_summary_renaming()

  parse_numeric_fields(x, fieldnames) %>%
    spread_and_rename_cutadapt_fieldnames(fieldnames)
}

###############################################################################

# ---- Helpers

###############################################################################

#' Define the renaming of the fields present in a cutadapt log summary section
#'
#' @importFrom   tibble        tribble
#'

define_cutadapt_summary_renaming <- function() {
  tibble::tribble(
    ~expected, ~output,
    "Total read pairs processed", "rp_input",
    "Read 1 with adapter", "r1_with_adapter",
    "Read 2 with adapter", "r2_with_adapter",
    "Pairs that were too short", "rp_too_short",
    "Pairs that were too long", "rp_too_long",
    "Pairs with too many N", "rp_too_many_n",
    "Pairs written (passing filters)", "rp_output",
    "Total basepairs processed", "bp_input",
    "Read 1", "bp_input_r1",
    "Read 2", "bp_input_r2",
    "Total written (filtered)", "bp_output",
    "Read 1", "bp_output_r1",
    "Read 2", "bp_output_r2"
  )
}

spread_and_rename_cutadapt_fieldnames <- function(x, fieldnames) {
  reformat_field_names <- function(x, fieldnames) {
    # x is a vector of strings
    rows <- match(x, fieldnames$expected)
    newfields <- fieldnames$output[rows]
    # a cutadapt summary contains some duplicated fieldnames
    read_idx <- which(x %in% c("Read 1", "Read 2"))
    read_parent <- read_idx - ifelse(x[read_idx] == "Read 1", 1, 2)
    read_suffix <- c("_r1", "_r2")[ifelse(x[read_idx] == "Read 1", 1, 2)]
    newfields[read_idx] <- paste(
      newfields[read_parent], read_suffix,
      sep = ""
    )
    newfields
  }

  if (!isTRUE(all(x$field %in% fieldnames$expected))) {
    stop(
      "All numeric fields from the text should be in the fieldnames dataframe"
    )
  }

  x %>%
    dplyr::mutate_(
      field = ~reformat_field_names(field, fieldnames)
    ) %>%
    dplyr::select_(.dots = c("field", "value")) %>%
    tidyr::spread_(key_col = "field", value_col = "value")
}
