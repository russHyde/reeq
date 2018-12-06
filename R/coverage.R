###############################################################################
#' Checks the validity of input to k-coverage functions
#'
#' Helper function - it just checks that a numeric matrix or dataframe is
#'   passed in and that a number is passed in; ensures validity of arguments to
#'   the k-coverage functions
#'
#' @param        mat           Some numeric matrix or all-numeric data.frame
#' @param        k             A single number
#'
#' @rdname       helper_check_validity_k_covered

check_validity_k_covered <- function(
                                     mat,
                                     k) {
  # `mat` should be a matrix or a data.frame (containing all numeric columns)
  stopifnot(
    is.matrix(mat) ||
      (is.data.frame(mat) && all(vapply(mat, is.numeric, logical(1))))
  )
  # `k` should be numeric
  stopifnot(is.numeric(k) && length(k) == 1)
  return(TRUE)
}

#' runner_k_covered
#'
#' Runs a function over a matrix/data.frame/DGEList/ExpressionSet. Returns a
#' subset of the row indices for the matrix-like datastructure.
#'
#' @param        mat           Either a matrix, data.frame, DGEList or
#'   ExpressionSet.
#' @param        k             A single number.
#' @param        lambda        Some function lambda(mat, k).
#'
#' @importFrom   Biobase       ExpressionSet   exprs
#' @importFrom   edgeR         DGEList   getCounts
#' @importFrom   magrittr      %>%   set_names
#' @importFrom   methods       is
#'
#' @rdname       helper_k_covered_docs

runner_k_covered <- function(
                             mat,
                             k,
                             lambda) {
  m <- if (methods::is(mat, "DGEList")) {
    edgeR::getCounts(mat)
  } else if (methods::is(mat, "ExpressionSet")) {
    Biobase::exprs(mat)
  } else {
    mat
  }

  check_validity_k_covered(m, k)

  lambda(m, k) %>%
    magrittr::set_names(NULL)
}

###############################################################################

#' which_k_covered_each_sample
#'
#' Returns the row index for those rows in a matrix where the row-sum is >= k.
#'
#' @inheritParams   runner_k_covered
#'
#' @export
#'
which_k_covered_each_sample <- function(
                                        mat,
                                        k = 1) {
  lambda <- function(m, k) {
    which(rowSums(m >= k) == ncol(m))
  }
  runner_k_covered(mat, k, lambda)
}

###############################################################################

#' which_k_covered_across_samples
#'
#' Returns the row indices for those rows in a matrix where every element in
#' that row is >= k
#'
#' @inheritParams   runner_k_covered
#'
#' @export
#'
which_k_covered_across_samples <- function(
                                           mat,
                                           k = 2 * ncol(mat)) {
  lambda <- function(m, k) {
    which(rowSums(m) >= k)
  }
  runner_k_covered(mat, k, lambda)
}

###############################################################################
