#' Obtain a design matrix from a file.
#'
#' Ensures that the samples mentioned in the rows of the design matrix are all
#' present in an expression dataset (if provided) and reorders the rows of the
#' design to make this happen.
#'
#' Here, we assume that the first column in the `.tsv` contains the sample
#' identifers, as used throughout the _R_ workflow. Each entry in that column
#' should be a single string that can directly be used as a column name.
#' All other columns should be numeric, and will be converted to a matrix.
#'
#' @param        filepath      A tab-separated file that defines the
#'   experimental design for the current experiment. This function fails if the
#'   file does not exist.
#'
#' @param        dge           A DGEList containing the expression data for the
#'   current experiment. This is optional. But providing a DGEList helps sanity
#'   check the design matrix - ensuring the samples in the design are all
#'   present in the correct order etc.
#'
#' @importFrom   readr         read_tsv   cols   stop_for_problems
#' @importFrom   magrittr      set_rownames
#' @importFrom   methods       is
#'
#' @export

import_design <- function(filepath, dge) {
  # assumes the sample identifier is the first column

  stopifnot(file.exists(filepath))

  df <- readr::read_tsv(
    filepath,
    comment = "#", col_types = readr::cols()
  )
  readr::stop_for_problems(df)

  design <- df[, -1] %>%
    as.matrix() %>%
    magrittr::set_rownames(df[[1]])

  # check that each sample mentioned in the design is present in the expression
  # dataset
  if (!missing(dge)) {
    stopifnot(methods::is(dge, "DGEList"))
    stopifnot(all(rownames(design) %in% colnames(dge)))
    design <- design[colnames(dge), ]
  }

  design
}