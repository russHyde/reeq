###############################################################################

# Helper functions and utils for use in `testthat` tests

###############################################################################

.df <- function(...) data.frame(..., stringsAsFactors = FALSE)

###############################################################################

dge <- function(m, add_genes = FALSE) {
  if (is.null(dimnames(m))) {
    features <- paste0("g", seq(nrow(m)))
    samples <- paste0("s", seq(ncol(m)))
    dimnames(m) <- list(features, samples)
  }
  gene_df <- if (add_genes) {
    .df(
      feature_id = rownames(m), length = 10, row.names = rownames(m)
    )
  } else {
    NULL
  }
  edgeR::DGEList(counts = m, genes = gene_df)
}

###############################################################################
