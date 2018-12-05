#' feature_counts_to_dgelist
#'
#' @param        fcounts_df    A data.frame containing the featureCounts (ie,
#'   read number) for each gene (rows) in each sample (columns 3:END). The
#'   first two columns should be "Geneid" and "Length". "Length" should be a
#'   numeric column. The column names are
#'   taken to be the working ids for the samples (except Geneid and Length) and
#'   each column name must be present once in the `id_column` (may be the
#'   rownames) of `sample_df`.
#'
#' @param        sample_df     A data.frame containing the sample-level
#'   annotation data, eg, the phenotypic info, the file locations, for each
#'   sample. This must either contain `id_column` as a colname or have defined
#'   rownames if `id_column` is NULL. The samples defined in the `id_column` of
#'   `sample_df` must be a superset of the samples present in the colnames of
#'   `counts.df`. The rows of sample_df that correspond to the samples in
#'   `fcounts_df` will be passed into the DGEList::samples slot, possibly being
#'   reordered to match the ordering of the samples in `fcounts_df`.
#'
#' @param        id_column     A single string giving the column of sample_df
#'   that gives the working ids for the samples in the experiment. If NULL, the
#'   ids are taken from the rownames of `sample_df` (which must be defined in
#'   that case).
#'
#' @param        ...           Additional parameters for passing to DGEList.
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
    stop("Can't uniquely map sample.ids between fcounts_df and sample_df")
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
    magrittr::set_rownames(fcounts_df$Geneid) %>%
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
      row.names = fcounts_df$Geneid,
      stringsAsFactors = FALSE
    ),
    ...
  )
}
