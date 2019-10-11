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

get_dge1 <- function(add_genes = FALSE) {
  # nolint start
  dge(matrix(1:10, nrow = 5), add_genes = add_genes)
  # nolint end
}

get_dge2 <- function() {
  my_dge <- get_dge1(TRUE)
  my_dge$genes$gc_percent <- seq(0.1, 0.9, length.out = nrow(my_dge))
  my_dge
}

###############################################################################
