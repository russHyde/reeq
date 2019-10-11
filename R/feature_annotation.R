###############################################################################

#' Obtain the GC percentage for a set of genes from a biomaRt database
#'
#' @param        feature_ids   A set of gene IDs.
#' @param        mart          A biomaRt `Mart` Dataset object.
#' @param        feature_column   The column of the `mart` dataset that
#'   contains IDs of the same type as the `feature_ids`.
#' @param        gc_column     The column of the `mart` dataset that contains
#'   the GC percentage values for the genes / features.
#'
#' @importFrom   biomaRt       getBM
#' @importFrom   dplyr         mutate
#' @importFrom   rlang         .data
#'
#' @export

get_gc_percent <- function(
                           feature_ids,
                           mart,
                           feature_column = "ensembl_gene_id",
                           gc_column = "percentage_gc_content") {
  if (missing(feature_ids) || missing(mart)) {
    stop("`feature_ids` and `mart` should be defined in `get_gc_percent`")
  }

  bm <- biomaRt::getBM(
    attributes = c(feature_column, gc_column),
    filters = feature_column,
    values = feature_ids,
    mart = mart
  ) %>%
    magrittr::set_colnames(c("feature_id", "gc_percent"))

  data.frame(
    feature_id = feature_ids,
    stringsAsFactors = FALSE
  ) %>%
    merge(bm, by = "feature_id", all.x = TRUE) %>%
    dplyr::mutate(gc_percent = as.numeric(.data[["gc_percent"]]))
}

###############################################################################
