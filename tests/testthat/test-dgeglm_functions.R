###############################################################################

context("Tests for `DGEGLM` manipulation functions")

###############################################################################

test_that("compute_fitted_values on a DGEGLM", {
  # User provides:
  # - a DGEGLM
  # - a contrasts matrix / df
  # - center: should the rows in the returned values be mean-subtracted

  # Note that edgeR stores the coefficient matrix in base-e, but users want
  # fitted values with-respect-to base-2

  # Expectations:
  # - Dimnames of output: rownames from rownames of DGEList, colnames from
  # colnames of contrasts matrix
  # - Mean, Sum and Select for the coefficients can be defined as a contrast
  # and return the expected values
  # - If center = TRUE, then rowSums ~ 0
  # - Rownames of contrasts matrix should be identical to colnames of the
  # `coefficients` entry of the DGEList

  n_groups <- 2
  n_samples_per_group <- 3
  n_samples <- n_groups * n_samples_per_group
  n_genes <- 100
  samples <- paste0("s", seq(n_samples))
  genes <- paste0("g", seq(n_genes))
  group <- factor(rep(letters[seq(n_groups)], each = n_samples_per_group))
  design <- model.matrix(~ -1 + group)
  counts <- matrix(
    rpois(n_genes * n_samples, lambda = 1000),
    nrow = n_genes,
    dimnames = list(genes, samples)
  )
  dge <- edgeR::DGEList(counts)
  fit <- edgeR::glmFit(edgeR::estimateGLMCommonDisp(dge), design = design)

  design_coefs_base2 <- log2(exp(fit[["coefficients"]]))

  contrasts <- as.matrix(data.frame(
    "sum_of_groups" = c(1, 1),
    "group_a" = c(1, 0),
    "group_b" = c(0, 1),
    "difference_of_groups" = c(-1, 1)
  ))

  expected <- as.matrix(data.frame(
    "sum_of_groups" = rowSums(design_coefs_base2),
    "group_a" = design_coefs_base2[, 1],
    "group_b" = design_coefs_base2[, 2],
    "difference_of_groups" = design_coefs_base2[, 2] - design_coefs_base2[, 1],
    row.names = row.names(fit)
  ))

  expect_equal(
    object = compute_fitted_values(fit, contrasts, center = FALSE),
    expected = expected,
    info = paste(
      "get fitted values for a set of contrasts, without centering the result"
    )
  )

  expect_equal(
    object = compute_fitted_values(fit, contrasts, center = TRUE),
    expected = expected - rowMeans(expected),
    info = "get fitted values for a set of contrasts, with row centering"
  )

  expect_error(
    object = compute_fitted_values(
      fit,
      contrasts = matrix(
        1,
        nrow = 1 + nrow(contrasts), ncol = 1,
        dimnames = list(c(rownames(contrasts), "zyx"), "a")
      )
    ),
    info = paste(
      "number of rows in contrasts should match number of coefficients in fit"
    )
  )

  expect_error(
    object = compute_fitted_values(
      fit,
      contrasts = magrittr::set_rownames(
        contrasts, rev(letters[seq(nrow(contrasts))])
      )
    ),
    info = paste(
      "If defined, the rownames of contrasts should match colnames of",
      "coefficients in fit"
    )
  )

  expect_error(
    object = compute_fitted_values(
      # mimics the relevant bits of the structure of a DGE fit, but is not a
      # DGE fit object:
      list(
        coefficients = matrix(
          1,
          nrow = n_genes, ncol = ncol(design),
          dimnames = list(genes, colnames(design))
        )
      ),
      contrasts = contrasts
    ),
    info = "x should be a DGE fit object (eg, DGEGLM) in compute_fitted_values"
  )
})
