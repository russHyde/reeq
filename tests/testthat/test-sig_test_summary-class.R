context("Tests for hypothesis-test functions")

test_that("significance summary for an edgeR dataset", {
  n_genes <- 100
  n_samples <- 6

  counts <- matrix(
    1000 * seq(n_genes * n_samples),
    nrow = n_genes,
    ncol = n_samples,
    dimnames = list(
      paste0("gene_", seq(n_genes)),
      paste0("patient_", seq(n_samples))
    )
  )
  groups <- rep(c("A", "B"), each = n_samples / 2)
  design <- model.matrix(~ -1 + groups)
  contrasts <- limma::makeContrasts(
    b_vs_a = "groupsB - groupsA",
    levels = design
  )
  dge <- edgeR::DGEList(counts)
  fit <- edgeR::glmFit(edgeR::estimateGLMCommonDisp(dge), design = design)

  expect_equal_dgelrt(
    get_sig_test_summary(fit, contrast = contrasts)$lrt,
    edgeR::glmLRT(fit, contrast = contrasts),
    info = "glmLRT and get_sig_test_summary should give the same LRT results",
    tolerance = 1e-8
  )
})
