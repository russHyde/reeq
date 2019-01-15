###############################################################################

#' is_single_string
#'
#' @noRd

is_single_string <- function(.x) {
  is.character(.x) && length(.x) == 1
}

###############################################################################

#' is_nonempty_df - checks if an R object is a data-frame and has at least one
#' column and at least one row
#'
#' @param        .x            Any R object
#' @return       Boolean: is `.x` a non-empty data-frame
#'
#' @noRd

is_nonempty_df <- function(.x) {
  if (missing(.x)) {
    stop(".x should be defined in is_nonempty_df")
  }
  is.data.frame(.x) &&
    nrow(.x) > 0 &&
    ncol(.x) > 0
}

#' is_nonempty_list - checks if an R object is a list and is non-empty
#'
#' @param        .x            Any R object
#' @return       Boolean: is .x a non-empty list?
#'
#' @noRd

is_nonempty_list <- function(.x) {
  if (missing(.x)) {
    stop(".x should be defined in is_nonempty_list")
  }
  is.list(.x) &&
    length(.x) > 0
}

###############################################################################

#' - `fcounts_df` should have first two cols "feature_id" and "length".
#' - Both length and all cols other than "feature_id" (of which there should be
#' >= 1) should be numeric.
#' - This only tests the structure of the data.frame; each count column
#' should also correspond to a sample id in sample_df, this is tested later.
#'
#' @param        df            Some data-frame
#'
#' @noRd

is_valid_fcounts_df <- function(df) {
  is_nonempty_df(df) &&
    ncol(df) >= 3 &&
    all(
      colnames(df)[1:2] == c("feature_id", "length")
    ) &&
    is.numeric(df[["length"]]) &&
    all(
      vapply(df[, -(1:2)], is.numeric, logical(1))
    )
}
