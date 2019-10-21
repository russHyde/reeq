#' get_lrt
#'
#' @importClassesFrom   edgeR   DGEGLM
#' @importFrom   edgeR         glmLRT
#' @noRd
#'
get_lrt <- function(fit, testable_features, contrast) {
  stopifnot(is(fit, "DGEGLM"))
  stopifnot(all(testable_features %in% rownames(fit)))
  sub_fit <- fit[testable_features, ]
  edgeR::glmLRT(sub_fit, contrast = contrast)
}

#' get_sig_features
#'
#' @importFrom   limma         decideTests
#' @noRd
#'
get_sig_features <- function(lrt, p_threshold) {

  # TODO: expand to use limma::decideTests.MArrayLM
  stopifnot(is(lrt, "DGEExact") || is(lrt, "DGELRT"))
  rows <- which(limma::decideTests(lrt, p.value = p_threshold) != 0)
  rownames(lrt)[rows]
}

#' get_sig_test_summary
#'
#' @param        fit           A DGEGLM from a edgeR model fit.
#' @param        testable_features   Which rows of the DGEGLM should be
#'   considered?
#' @param        contrast      A matrix or vector defining the contrasts over
#'   which the LRT test should summarise.
#' @param        p_threshold   A single float given the FDR threshold under
#'   which a given feature is considered 'significant'.
#'
#' @importFrom   edgeR         topTags
#'
#' @export

get_sig_test_summary <- function(fit,
                                 testable_features = row.names(fit),
                                 contrast,
                                 p_threshold = 0.05) {
  lrt <- get_lrt(fit, testable_features, contrast)
  sig_features <- get_sig_features(lrt, p_threshold)
  top_table <- edgeR::topTags(lrt, n = Inf)

  structure(
    list(
      features = testable_features,
      sig_features = sig_features,
      lrt = lrt,
      top_table = top_table,
      num_features = length(testable_features),
      num_sig_features = length(sig_features),
      p_threshold = p_threshold
    ),
    class = "sig_test_summary"
  )
}
