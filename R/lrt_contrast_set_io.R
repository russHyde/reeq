#' parse_lrt_table
#'
#' Each likelihood-ratio test is performed over a set of contrasts. This
#' function extracts the names of the contrasts over which each LRT should be
#' performed from a data.frame. The collection of contrasts should be defined
#' as a ';'-separated string in the input data.frame.
#'
#' A tibble where the contrast sets are present in a list-column is returned.
#'
#' The number of siblings of each LRT-node is also determined.
#'
#' The `level` should be 1 for the LRT-nodes comprised of the largest
#' collection of contrasts and then increments for each proper subset of these
#' nodes. level-1 nodes have no parent, so set parent = NA for these rows.
#'
#' The tibble that is returned is ordered by increasing 'level', to try to
#' ensure that parent nodes occur before children.
#'
#' We allow a node to have at most one parent. For use in nested-LRT analysis,
#' any node that has a missing (NA) parent must be at level-1.
#'
#' @param   x   A data-frame containing columns "lrt_name", "level", "parent",
#'   "contrast_set". The `contrast_set` column must be a vector of
#'   ';'-separated characters
#'
#' @importFrom   dplyr      add_count   rename
#' @importFrom   magrittr   %>%
#' @importFrom   rlang      .data
#'
#' @export
#'
parse_lrt_table <- function(x) {
  # TODO: for each level-k (k>=2) node, ensure it's parent is present at level
  # k-1.
  #

  reqd_columns <- c("lrt_name", "level", "parent", "contrast_set")
  stopifnot(all(reqd_columns %in% colnames(x)))

  # Convert the ';'-separated character column "contrast_set" into a list of
  # character vectors
  x$contrast_set <- strsplit(x$contrast_set, split = ";")

  # Then order the rows so that 'level' is non-decreasing (we want parents to
  # occur before children) and so the number of nodes that share the same
  # parent node is indicated (n_siblings).
  x %>%
    dplyr::arrange(.data[["level"]]) %>%
    dplyr::add_count(.data[["parent"]]) %>%
    dplyr::rename(n_siblings = .data[["n"]])
}
