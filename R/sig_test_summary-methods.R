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

export.default <- function(x, ...) {
  stop("Not implemented")
}

#'

export.sig_test_summary <- function(x,
                                    output_dir,
                                    test_name,
                                    test_type,
                                    ...) {
  p_threshold <- sprintf("p%.4f", x$p_threshold)

  # obtain the top-table for all features
  tt_file <- paste(
    test_name, test_type, "top_tags", "tsv",
    sep = "."
  )
  tt_path <- file.path(output_dir, tt_file)
  tt_data <- as.data.frame(x$top_table)

  # obtain the list of all significant features
  gene_file <- paste(
    test_name, test_type, "sig_features", p_threshold, "tsv",
    sep = "."
  )
  gene_path <- file.path(output_dir, gene_file)
  gene_data <- x$sig_features

  # write out files containing the results
  readr::write_tsv(tt_data, tt_path)
  write(gene_data, gene_path, ncolumns = 1)
}
