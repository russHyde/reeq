#' Export a data structure to a (set of) file(s)
#'
#' @param        x             Some datastructure. Typically a
#'   `sig_test_summary` as created by `reeq::get_sig_test_summary`.
#' @param        ...           Further arguments. Currently unused.
#'
#' @export

export <- function(x, ...) {
  UseMethod("export")
}

#' Default function for exporting a data-structure.
#'
export.default <- function(x, ...) {
  stop("Not implemented")
}

#' Export a `sig_test_summary` object to files
#'
#' Write the (unfiltered) top-table and significant gene lists from a
#' `sig_test_summary` object to some files. The file format is
#' `<output_dir>/<test_name>.top_tags.tsv` or
#' `<output_dir>/<test_name>.sig_features.<p_threshold>.tsv`, where the p-value
#' threshold is obtained from inside the `sig_test_summary` object and is
#' formatted to 4 decimal places.
#'
#' @inheritParams   export
#'
#' @param        output_dir    Which directory should the results be output to
#' @param        test_name     A name for the hypothesis test that is
#'   summarised in the dataset `x`.
#'
#' @importFrom   readr         write_tsv
#'

export.sig_test_summary <- function(x,
                                    output_dir,
                                    test_name,
                                    ...) {
  p_threshold <- sprintf("p%.4f", x$p_threshold)

  # obtain the top-table for all features
  tt_file <- paste(
    test_name, "top_tags", "tsv",
    sep = "."
  )
  tt_path <- file.path(output_dir, tt_file)
  tt_data <- as.data.frame(x$top_table)

  # obtain the list of all significant features
  gene_file <- paste(
    test_name, "sig_features", p_threshold, "tsv",
    sep = "."
  )
  gene_path <- file.path(output_dir, gene_file)
  gene_data <- x$sig_features

  # write out files containing the results
  readr::write_tsv(tt_data, tt_path)
  write(gene_data, gene_path, ncolumns = 1)
}
