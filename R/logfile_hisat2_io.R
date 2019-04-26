###############################################################################

# Functions for parsing alignment statistics from `hisat2` logfiles

###############################################################################

#' Extract the alignment summary statistics from a hisat2 logfile
#'
#' @param        x             A filepath. This should be a hisat2 logfile. The
#'   function dies if this is not an existing file of if more than one entry is
#'   present in \code{x}.
#'
#' @export
#'

import_hisat2_summary <- function(x){
  # validate the input is a single filepath
  stopifnot(length(x) == 1)
  stopifnot(file.exists(x))

  log_txt <- readr::read_file(x)

  parse_hisat2_summary(log_txt)
}

parse_hisat2_summary <- function(x){
  # Assumes the input text is from a --new-summary from hisat2
  # - Work out how to distinguish --new-summary from old summary files
  stopifnot(is.character(x))
  stopifnot(length(x) == 1)

  df <- parse_numeric_fields(x) %>%
    dplyr::mutate(field = c(
      "rp_input", "rp_zero", "rp_concordant_once", "rp_concordant_multi",
      "rp_discordant_once", "unpaired", "unpaired_unaligned", "unpaired_once",
      "unpaired_multi", "align_rate"
      ))

  tidyr::spread_(df, key_col = "field", value_col = "value")[df$field] %>%
    dplyr::mutate_(align_rate = ~ align_rate / 100)
}
