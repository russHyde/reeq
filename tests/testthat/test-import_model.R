###############################################################################

context(
  paste(
    "model design, contrasts, contrast-collection, lrt-tree and",
    "plot-comparisons can be imported"
  )
)

###############################################################################

test_that("a design matrix can be imported from .tsv", {
  # .tsv contains:
  # - any number >= 2 of columns
  # - a sample ID column with values matching the sample IDs in the DGEList

  # Constraints
  # - .tsv file must exist
  # - DGEList may be provided, and it is:
  #   - it should have sample IDs
  #   - rownames of the design should be a subset of the samples in the DGEList

  # happy path:
  # - sample IDs ("s1", "s2") are present in the DGEList
  # - sample IDs occur in the same order in 1st column of `design.tsv`

  # TODO:
  # - non-numeric entries in the design file -> error
  # - a design file that causes readr::read_tsv to fail

  design_file <- file.path("model_tables", "design.tsv")

  my_dge <- get_dge1(add_genes = TRUE)
  expected_design <- matrix(
    c(1, 1, 0, 1),
    nrow = 2,
    dimnames = list(c("s1", "s2"), c("intercept", "f2"))
  )

  expect_equivalent(
    import_design(design_file, dge = my_dge),
    expected_design,
    info = "import a design from a file & validate against a DGEList"
  )

  expect_equivalent(
    import_design(design_file),
    expected_design,
    info = "when no DGEList is provided, just import the design matrix"
  )

  expect_error(
    import_design(
      design_file, dge = magrittr::set_colnames(my_dge, c("X", "Y"))
    ),
    info = "sample IDs of the design should overlap with those in DGEList"
  )

  expect_equal(
    import_design(design_file, my_dge),
    import_design(design_file, my_dge[, rev(colnames(my_dge))]),
    info = paste(
      "the sample-order in the design matrix is determined by the file, not",
      "the DGEList"
    )
  )

  expect_error(
    import_design(file.path("NOT-A-FILE"), dge = my_dge),
    info = "design .tsv should exist"
  )

  # sample names are correct, but class should be DGEList
  df <- data.frame(s1 = 1, s2 = 2)
  expect_error(
    import_design(design_file, dge = df),
    info = "if `dge` is provided, it should be a DGEList"
  )
})

###############################################################################

test_that("a contrast matrix can be imported from a .tsv", {

  contrasts_file <- file.path("model_tables", "contrasts.tsv")

  treatments <- rep(letters[1:3], each = 4)
  design <- model.matrix(~ treatments) %>%
    magrittr::set_colnames(c("intercept", "b", "c"))

  expect_equal(
    import_contrasts(contrasts_file, design),
    matrix(
      c(0, 1, 0, 0, -1, 1),
      nrow = 3,
      dimnames = list(c("intercept", "b", "c"), c("b_vs_a", "c_vs_b"))
    ),
    info = "contrasts-matrix can be read from a .tsv"
  )

  expect_error(
    import_contrasts("NOT-A-FILE", design),
    info = "contrasts file should exist"
  )

  expect_error(
    import_contrasts(contrasts_file),
    info = "design must be defined in import_contrasts()"
  )

  expect_error(
    import_contrasts(
      contrasts_file, magrittr::set_colnames(design, letters[26:24])
    ),
    info = "model coefficient names should match between contrasts and design"
  )
})
