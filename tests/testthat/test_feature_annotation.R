###############################################################################

context("Tests for feature-annotation functions")

###############################################################################

test_that("get_gc_percent: invalid input", {
  expect_error(
    get_gc_percent(),
    info = "no input to `get_gc_percent`"
  )
})

test_that("get_gc_percent: string input", {
  # Genes that are missing from the database do not return a row in a getBM
  # query
  genes <- c("not_a_gene", "valid_gene")
  gc_percent_string <- c(NA, "38.5")
  gc_percent_real <- c(NA, 38.5)

  mock_get_bm <- mockery::mock(
    .df(
      ensembl_gene_id = "valid_gene", "percentage_gc_content" = "38.50"
    )
  )
  mock_mart <- structure(list(), class = "Mart")

  testthat::with_mock(
    getBM = mock_get_bm, {
      expect_equal(
        get_gc_percent(
          feature_ids = genes, mart = mock_mart,
          feature_column = "ensembl_gene_id"
        ),
        .df(
          feature_id = genes,
          gc_percent = gc_percent_real
        ),
        info = "parse gc-percentage from a character-column"
      )
    },
    .env = "biomaRt"
  )

  mockery::expect_called(mock_get_bm, 1)
  mockery::expect_args(
    mock_get_bm,
    1,
    attributes = c("ensembl_gene_id", "percentage_gc_content"),
    filters = "ensembl_gene_id",
    values = genes,
    mart = mock_mart
  )
})

test_that("get_gc_percent: numeric input", {
  # Genes that are missing from the database do not return a row in a getBM
  # query
  genes <- c("not_a_gene", "valid_gene")
  gc_percent_real <- c(NA, 38.5)

  mock_get_bm <- mockery::mock(
    .df(
      ensembl_gene_id = "valid_gene", "percentage_gene_gc_content" = 38.50
    )
  )
  mock_mart <- structure(list(), class = "Mart")

  testthat::with_mock(
    getBM = mock_get_bm, {
      expect_equal(
        get_gc_percent(
          feature_ids = genes, mart = mock_mart,
          feature_column = "ensembl_gene_id",
          gc_column = "percentage_gene_gc_content"
        ),
        .df(
          feature_id = genes,
          gc_percent = gc_percent_real
        ),
        info = "parse gc-percentage from a numeric-column"
      )
    },
    .env = "biomaRt"
  )

  mockery::expect_called(mock_get_bm, 1)
  mockery::expect_args(
    mock_get_bm,
    1,
    attributes = c("ensembl_gene_id", "percentage_gene_gc_content"),
    filters = "ensembl_gene_id",
    values = genes,
    mart = mock_mart
  )
})
