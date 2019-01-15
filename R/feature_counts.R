#' Combine a list of feature-counts data-frames into an `edgeR` `DGEList`
#'
#' @param        fcounts_df    A `data.frame` containing the feature-counts
#'   (ie, read number) for each feature (gene; rows) in each sample (columns
#'   3:END). The first two columns should be `feature_id` and `length`.
#'   `length` should be a numeric column. The column names are taken to be the
#'   working IDs for the samples (except `feature_id` and `length`) and each
#'   column name must be present once in the `id_column` (may be the
#'   `rownames`) of `sample_df`.
#'
#' @param        sample_df     A `data.frame` containing the sample-level
#'   annotation data, eg, the phenotypic info, the file locations, for each
#'   sample. This must either contain `id_column` as a column name or have
#'   defined `rownames` if `id_column` is `NULL`. The samples defined in the
#'   `id_column` of `sample_df` must be a superset of the samples present in
#'   the colnames of `fcounts_df`. The rows of `sample_df` that correspond to
#'   the samples in `fcounts_df` will be passed into the `DGEList::samples`
#'   slot, possibly being reordered to match the ordering of the samples in
#'   `fcounts_df`.
#'
#' @param        id_column     A single string giving the column of `sample_df`
#'   that gives the working IDs for the samples in the experiment. If `NULL`,
#'   the IDs are taken from the `rownames` of `sample_df` (which must be
#'   defined in that case).
#'
#' @param        ...           Additional parameters for passing to `DGEList`.
#'   Note that these should be ordered according to the sample-ordering in
#'   `fcounts_df` (if this differs from that in `sample_df`).
#'
#' @importFrom   edgeR         DGEList
#' @importFrom   magrittr      %>%   set_rownames   set_colnames
#'
#' @export

feature_counts_to_dgelist <- function(
                                      fcounts_df,
                                      sample_df,
                                      id_column = NULL,
                                      ...) {
  stopifnot(is_valid_fcounts_df(fcounts_df))

  stopifnot(is_nonempty_df(sample_df))

  sample_df_ids <- if (is.null(id_column)) {
    rownames(sample_df)
  } else {
    stopifnot(
      is_single_string(id_column) &&
        id_column %in% colnames(sample_df)
    )
    sample_df[[id_column]]
  }
  if (any(duplicated(sample_df_ids))) {
    stop("Can't uniquely map sample IDs between `fcounts_df` and `sample_df`")
  }

  fcounts_df_ids <- colnames(fcounts_df)[-(1:2)]
  stopifnot(all(fcounts_df_ids %in% sample_df_ids))

  sample_df_reorder <- match(
    fcounts_df_ids,
    sample_df_ids
  )

  counts <- as.matrix(
    fcounts_df[, -(1:2)]
  ) %>%
    magrittr::set_rownames(fcounts_df$feature_id) %>%
    magrittr::set_colnames(fcounts_df_ids)

  edgeR::DGEList(
    counts = counts,
    samples = data.frame(
      sample_df[
        sample_df_reorder,
        seq_along(sample_df),
        drop = FALSE
      ],
      row.names = sample_df_ids[sample_df_reorder],
      stringsAsFactors = FALSE
    ),
    genes = data.frame(
      fcounts_df[, 1:2],
      row.names = fcounts_df$feature_id,
      stringsAsFactors = FALSE
    ),
    ...
  )
}

###############################################################################

#' Combine a list of feature-count `data.frame`s into a single feature-count
#' `data.frame`
#'
#' The input `data.frame`s should have columns including `feature_id`, `length`
#' and a set of file paths. The file path column contains the read count for a
#' given feature in that sample. Returns a `data.frame` with the count column
#' for each sample juxtaposed and with `feature_id` and `length` columns.
#'
#' @param        count_list    A `list` of `data.frame`s. Each `data.frame`
#'   should have columns `feature_id`, `length` and a final column that
#'   contains the read counts for some sample.
#'
#' @param        names_as_colnames   Boolean. If `TRUE`, the column IDs in the
#' output `data.frame` are taken from the `list` names in `count_list`.
#'
#' @importFrom   magrittr      %>%
#' @importFrom   stats         setNames
#'
#' @export

cbind_feature_counts <- function(
                                 count_list,
                                 names_as_colnames = FALSE) {
  # - count_list is a list of data.frames
  # - Each data.frame should have a feature_id and a length column

  # - The output should be a data.frame with feature_id, length columns and then
  # the counts from each data.frame in count_list in turn

  # Input should be a list
  if (missing(count_list) || !is_nonempty_list(count_list)) {
    stop("Input to merge_feature_counts should be a list of data.frames")
  }

  is_df_valid <- function(df) {
    is_valid_fcounts_df(df) && ncol(df) == 3
  }

  # Each entry should be a data.frame
  if (!all(vapply(count_list, is_df_valid, logical(1)))) {
    stop(
      paste(
        "Each entry in count_list should be a data.frame with three cols:",
        "feature_id, length and a count column (the latter can have any name",
        "but must be a numeric column)"
      )
    )
  }

  # If names_as_colnames = TRUE, there should be names in the count_list
  if (names_as_colnames && is.null(names(count_list))) {
    stop("names(count_list) should be defined if names_as_colnames = TRUE")
  }

  # Reset the count column to be the same as the name of the corresponding
  # entry in the count_list
  if (names_as_colnames) {
    count_list <- Map(
      function(df, cname) {
        stats::setNames(df, c("feature_id", "length", cname))
      },
      count_list,
      names(count_list)
    )
  }

  # Merge the data.frames together
  Reduce(
    function(x, y) {
      merge(x, y, by = c("feature_id", "length"))
    },
    count_list
  ) %>%
    as.data.frame(stringsAsFactors = FALSE)
}

###############################################################################
