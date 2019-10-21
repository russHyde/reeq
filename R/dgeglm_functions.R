###############################################################################

#' Obtain fitted values for a given combination of design coefficients
#'
#' @param        x             A fit from running an edgeR pipeline, eg, a
#'   DGEGLM object.
#' @param        contrasts     A matrix defining linear combination of design
#'   coefficients. eg, this might define the combination of design coefs that
#'   gives the fitted value for a particular (set of) treatment group. This
#'   should be specified in the format of a contrasts matrix (a row for each
#'   design coefficient, a column for each computed value). If rownames are
#'   specified, they should be identical to those used in the corresponding
#'   design matrix.
#' @param        center       Boolean: should the mean value be subtracted from
#'   each row before returning the results matrix?
#'
#' @export
compute_fitted_values <- function(x, contrasts, center = FALSE) {
  stopifnot(methods::is(x, "DGEGLM"))

  if (!is.null(rownames(contrasts))) {
    stopifnot(
      all(rownames(contrasts) == colnames(x[["coefficients"]]))
    )
  }

  fitted_values <- log2(exp(1)) * (x[["coefficients"]] %*% contrasts)

  if (center) {
    fitted_values - rowMeans(fitted_values)
  } else {
    fitted_values
  }
}

###############################################################################
