###############################################################################

# utility functions

###############################################################################

#' replace_with
#'
#' For each element in \code{x}, if it is present in the `search_list`, replace
#' it with the corresponding value in the `return_list`. Any element that is
#' not in the `search_list` is returned unaltered.
#'
#' @param        x             A vector.
#' @param        search_list   A vector.
#' @param        return_list   A vector. This should be the same length as
#'   \code{search_list}.
#' @param        strict        Boolean. Should every element of `x` be present
#'   in \code{search_list}?
#'
#' @return       A vector of the same length as \code{x}, containing some
#'   values from \code{x} and some values from \code{return_list} (so the class
#'   of the vector may be different from that for `x`).
#'
#' @export

replace_with <- function(x, search_list, return_list, strict = FALSE) {
  stopifnot(
    length(search_list) == length(return_list)
  )
  if (strict && !all(x %in% search_list)) {
    stop("When `strict`==FALSE, all members of `x` should be in `search_list`")
  }
  if (length(x) == 0) {
    return(x)
  }
  ifelse(
    x %in% search_list,
    return_list[match(x, search_list)],
    x
  )
}
