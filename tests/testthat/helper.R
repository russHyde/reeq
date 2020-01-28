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

expect_equal_dgelrt <- function(object,
                                expected,
                                tolerance = 1e-8,
                                info = NULL) {
  # capture object and it's label
  act <- testthat::quasi_label(rlang::enquo(object), arg = "object")

  # call expect
  ns <- union(names(object), names(expected))
  for (n in ns) {
    testthat::expect(
      n %in% names(object) &&
        all.equal(object[[n]], expected[[n]], tolerance = tolerance),
      failure_message = paste0(
        "Failed on ", n
      ),
      info = info
    )
  }

  invisible(act$val)
}
###############################################################################
