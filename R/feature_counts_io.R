###############################################################################

# Functions for importing tables from featureCounts

###############################################################################

#' Import read counts from a set of featureCounts files and combine them into
#' a single `data.frame`
#'
#' The files may be of standard format or of short `featureCounts` format (see
#'   below). Only the `Geneid`, `Length` and counts are kept and the `Geneid`
#'   and `Length` columns are renamed `feature_id` and `length`.
#'
#' Standard format contains: `Geneid`, `Chr`, `Start`, `End`, `Strand`,
#'   `Length`, counts
#'
#' Short format contains: `Geneid`, `Length`, counts
#'
#' @param        files         A vector of file-paths. Each should be an output
#'   from `featureCounts` and contain, at least, the `Geneid`, `Length` and
#'   feature-count columns.
#' @param        sample_ids    A vector of sample IDs corresponding, in order,
#'   to the entries in `files`. If `NULL` ...
#' @param        comment,col_types,progress   Arguments passed through to
#'   `readr::read_tsv` that specify the format of the input files and whether
#'   a progress bar should be shown.
#'
#' @importFrom   magrittr      set_names
#' @importFrom   purrr         map
#'
#' @export

read_feature_counts <- function(
                                files,
                                sample_ids,
                                comment = "#",
                                col_types = "cii",
                                progress = FALSE) {
  files %>%
    purrr::map(
      read_single_feature_counts_file,
      comment = comment, col_types = col_types, progress = progress
    ) %>%
    magrittr::set_names(sample_ids) %>%
    cbind_feature_counts(names_as_colnames = TRUE)
}

###############################################################################

#' Read a single featureCounts file
#'
#' @param        file          A single file, containing the `Geneid`, `Length`
#'   and count columns from a single run of `featureCounts`.
#' @inheritParams   read_feature_counts
#'
#' @importFrom   dplyr         rename
#' @importFrom   magrittr      %>%
#' @importFrom   readr         read_tsv
#' @importFrom   rlang         .data

read_single_feature_counts_file <- function(
                                            file,
                                            comment = "#",
                                            col_types = "cii",
                                            progress = FALSE) {
  # Remove columns that are present in the standard-format featureCounts but
  # which aren't of use in differential expression.
  drop_unrequired_columns <- function(df) {
    df[setdiff(colnames(df), c("Chr", "Start", "End", "Strand"))]
  }

  fcounts <- readr::read_tsv(
    file,
    comment = comment, col_types = col_types, progress = progress
  ) %>%
    drop_unrequired_columns() %>%
    dplyr::rename(
      feature_id = .data[["Geneid"]],
      length = .data[["Length"]]
    )

  stopifnot(is_valid_fcounts_df(fcounts))

  fcounts
}

###############################################################################
