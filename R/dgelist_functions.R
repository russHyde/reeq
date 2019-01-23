###############################################################################

#' filter_by_read_count
#'
#' @param        dge           An `edgeR::DGEList`.
#' @param        threshold     A single numeric value. This value decides which
#'   rows are kept in the `DGEList` that is returned. The value of the other
#'   arguments dictates how this threshold is applied: if `use_sample_average`
#'   is TRUE, rows that have an average value >= this threshold are kept;
#'   otherwise, rows that have a value >= this threshold in at least a fraction
#'   `fraction_of_samples` of the columns are kept. If the remaining args are
#'   left at default values, only rows where every column has a value >= this
#'   threshold are kept. Default: 0.
#' @param        count_type    A choice of `counts` (the default) or `cpm`.
#'   Should the current filter / threshold be applied to the raw counts or to
#'   the counts-per-million?
#' @param        fraction_of_samples   What fraction of the input samples (ie,
#'   columns) should meet the `threshold` for the filter to have passed for a
#'   given feature (ie, row)? Default: 1.
#' @param        use_sample_average    Should rows be filtered based on whether
#'   the average of the value across a row is >= the `threshold`? If this is
#'   `TRUE` then the `fraction_of_samples` value is ignored. Default: `FALSE`.
#' @param        keep_lib_sizes   Boolean. Indicates whether library-sizes
#'   should be recomputed for each sample after any features have been filtered
#'   out. Passed through to `edgeR` subsetting function.
#'
#' @importFrom   methods       is
#' @importFrom   edgeR         cpm
#' @importFrom   edgeR         [.DGEList
#'
#' @importClassesFrom   edgeR   DGEList
#'
#' @export
#'

filter_by_read_count <- function(dge,
                                 threshold = 0,
                                 count_type = c("counts", "cpm"),
                                 fraction_of_samples = 1,
                                 use_sample_average = FALSE,
                                 keep_lib_sizes = FALSE) {
  stopifnot(methods::is(dge, "DGEList"))
  stopifnot(is.numeric(threshold))
  stopifnot(is.logical(use_sample_average))

  count_type <- match.arg(count_type)

  count_matrix <- if (count_type == "cpm") {
    edgeR::cpm(dge)
  } else {
    dge$counts
  }

  keep_rows <- if (use_sample_average) {
    which_k_covered_across_samples(
      mat = count_matrix,
      k = threshold * ncol(dge)
    )
  } else {
    which_k_covered(
      mat = count_matrix,
      k = threshold,
      fraction_of_samples = fraction_of_samples
    )
  }

  dge[keep_rows, keep.lib.sizes = keep_lib_sizes]
}

###############################################################################