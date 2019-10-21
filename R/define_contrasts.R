#' define a contrasts matrix based on a vector of string definitions
#'
#' There is a slight problem with limma::makeContrasts.
#' `makeContrasts(contrasts = some_named_vector, levels = design)` does not set
#' the names of the resulting contrasts to match the names in
#' `some_named_vector`.
#'
#' Here we convert a vector of string-based contrast definitions and the names
#' for those contrasts into a contrasts matrix.
#'
#' @param        .contrasts    A character vector. All words within the entries
#'   should be columns in the design matrix. Any numbers will be treated as
#'   numbers. Examples, c("0", "1 + coef1", "coef2 - 2*coef1").
#' @param        .names        A character vector of the same length as
#'   `.contrasts`. These names will be the column names of the returned matrix.
#' @param        levels        A design matrix. The column names of this matrix
#'   will become the row names of the returned matrix
#'
#' @export
#'

define_contrasts <- function(.contrasts, .names, levels) {
  stopifnot(is.character(.contrasts))
  stopifnot(is.character(.names))
  # TODO: allow .design as a (solely numeric) data.frame
  stopifnot(is.matrix(levels))

  m <- limma::makeContrasts(
    contrasts = .contrasts,
    levels = levels
  )

  dimnames(m) <- list(
    colnames(levels),
    .names
  )

  m
}
