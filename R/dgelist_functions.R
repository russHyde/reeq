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

#' Add new columns from a table of feature-data into the feature-data
#' (`$genes`) of a `DGEList`
#'
#' If the `DGEList$genes` entry is undefined, the feature data is added as the
#' genes entry; otherwise the new columns are appended to the RHS of the
#' existing genes entry.
#'
#' Ensures that the feature-id / gene-id matches up between the provided
#' data-frame and the existing `DGEList`
#'
#' The column named by `feature_id` must be present in the `annotations` data
#' frame. If that column is present in the $genes entry of `dge`, matching will
#' be done on that column, otherwise, matching will be done on the rownames of
#' the DGEList.
#'
#' @param        dge           A DGEList
#' @param        annotations   A data-frame to either use as the `genes` entry
#'   or to append to the `genes` entry of the DGEList. Will match its
#'   `feature_id` column against the feature ordering in the DGEList.
#' @param        feature_id    Which column of `annotations` contains the
#'   IDs for the features. The contents of this column will be matched against
#'   a column of `dge$genes` - this is preferentially a column of the same name
#'   but, if that colname is absent from `dge$genes`, then the rownames of the
#'   DGEList will be used.
#'
#' @return       a DGEList
#'
#' @export

append_feature_annotations <- function(dge,
                                       annotations,
                                       feature_id = "feature_id") {
  stopifnot(is(dge, "DGEList"))
  stopifnot(
    is.data.frame(annotations) && feature_id %in% colnames(annotations)
  )

  # We preferentially match on the feature_id of dge$genes if it exists;
  #   otherwise we match the feature_id column of annotations with the rownames
  #   of dge.
  feature_target <- if (
    !is.null(dge$genes) && feature_id %in% colnames(dge$genes)
  ) {
    dge[["genes"]][[feature_id]]
  } else {
    rownames(dge)
  }

  # Get a reordering of `annotations` to match the row-order of `dge`.
  m <- match(feature_target, annotations[[feature_id]])

  # Extract all novel columns from the annotations dataset, reordering the rows
  #   to match the dge dataset. Append these novel columns to the dge dataset's
  #   `genes` dataframe (initialising if necessary).
  genes <- if (is.null(dge$genes)) {
    annotations[m, , drop = FALSE]
  } else {
    cols_to_add <- setdiff(
      colnames(annotations),
      colnames(dge$genes)
    )
    cbind(
      dge$genes,
      annotations[m, cols_to_add, drop = FALSE]
    )
  }
  rownames(genes) <- rownames(dge)

  dge$genes <- genes
  dge
}

###############################################################################

#' Add (or replace) an offset matrix in DGEList
#'
#' Trivial implementation, but useful in pipelines.
#'
#' @param        dge           A DGEList.
#' @param        offset        A matrix of the same dimensions as the counts in
#'   `dge`.
#'
#' @return       A DGEList.
#' @export

set_offset <- function(dge, offset) {
  stopifnot(methods::is(dge, "DGEList"))
  stopifnot(is.numeric(offset) && is.matrix(offset))
  stopifnot(all(dim(offset) == dim(dge)))

  dge$offset <- offset
  dge
}

###############################################################################

#' cqn_dgelist
#'
#' @param        x             A DGEList object
#' @param        length_column   Which column of x$genes contains the
#'   gene-lengths.
#' @param        gc_column     Which column of x$genes contains the GC content?
#' @param        lib_size_column   Which column of x$samples contains the
#'   library sizes?
#' @param        ...           Further arguments to pass along to `cqn::cqn`
#'
#' @importFrom   cqn           cqn
#' @importFrom   edgeR         getCounts
#'
#' @export

cqn_dgelist <- function(x,
                        length_column = "length",
                        gc_column = "gc_percent",
                        lib_size_column = "lib.size",
                        ...) {
  stopifnot(is(x, "DGEList"))
  stopifnot(
    length_column %in% colnames(x[["genes"]]) &&
      gc_column %in% colnames(x[["genes"]]) &&
      lib_size_column %in% colnames(x[["samples"]])
  )

  cqn::cqn(
    counts = edgeR::getCounts(x),
    x = x[["genes"]][[gc_column]],
    lengths = x[["genes"]][[length_column]],
    sizeFactors = x[["samples"]][[lib_size_column]]
  )
}

###############################################################################
