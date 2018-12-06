###############################################################################

context("Tests for coverage functions in `reeq`")

###############################################################################

.DF <- function(...) data.frame(..., stringsAsFactors = FALSE)

### ======================================================================= ###
#   Functions for QC-filtering on coverage
### ======================================================================= ###

test_that("check_validity_k_covered", {
  # 2 x 2 identity matrix
  I2 <- diag(rep(1, 2))
  # 2 x 2 all-ones upper-triangular
  T2 <- matrix(c(1, 0, 1, 1), nrow = 2)

  # `mat` should be a matrix or a data.frame
  expect_error(
    object = check_validity_k_covered(mat = "Not a matrix"),
    info = "mat should be a matrix in which_k_covered_*"
  )

  # If `mat` is a data.frame, it should have all-numeric columns
  expect_error(
    object = check_validity_k_covered(
      mat = data.frame(a = 1:2, b = letters[1:2]),
      k = 1
    ),
    info = "If mat is a data.frame, then it should be all-numeric"
  )
  # `k` should be a number
  expect_error(
    object = check_validity_k_covered(
      mat = I2,
      k = "Not a number"
    ),
    info = "k should be a number in which_k_covered_*"
  )
  expect_error(
    object = check_validity_k_covered(
      mat = I2,
      k = 1:2
    ),
    info = "k should be a single number in which_k_covered_*"
  )
})

###############################################################################

test_that("which_k_covered_each_sample", {
  # Assumes: genes are indexed on rows, samples on columns

  # 2 x 2 identity matrix
  I2 <- diag(rep(1, 2))
  # 2 x 2 all-ones upper-triangular
  T2 <- matrix(c(1, 0, 1, 1), nrow = 2)

  # Validity checks are done in .check_validity_k_covered

  # Return all rows where coverage >= k in every column

  # I2 should return empty sequence for k > 0
  expect_equal(
    object = which_k_covered_each_sample(mat = I2, k = 1),
    expected = integer(0),
    info = paste(
      "2x2 identity: all rows contain some zeros => no rows pass when k > 0"
    )
  )

  # I2 should return empty sequence for k > 0, even when passed as a data.frame
  expect_equal(
    object = which_k_covered_each_sample(mat = data.frame(I2), k = 1),
    expected = integer(0),
    info = "2x2 identity as data.frame"
  )

  # I2 should return c(1, 2) for k <= 0
  expect_equal(
    object = which_k_covered_each_sample(mat = I2, k = 0),
    expected = c(1, 2),
    info = paste(
      "2x2 identity: all rows contain some zeros => all rows pass when k <= 0"
    )
  )

  # T2 should return empty sequence for k > 1
  expect_equal(
    object = which_k_covered_each_sample(mat = T2, k = 2),
    expected = integer(0),
    info = "[[1 1], [0 1]]: all rows should fail when k > 1"
  )

  # T2 should return c(1) for k == 1
  expect_equal(
    object = which_k_covered_each_sample(mat = T2, k = 1),
    expected = 1,
    info = "[[1 1], [0 1]]: first row should pass when k == 1"
  )

  # T2 should return c(1, 2) for k <= 0
  expect_equal(
    object = which_k_covered_each_sample(mat = I2, k = 0),
    expected = c(1, 2),
    info = "[[1 1], [0 1]]: all rows should pass when k <= 0"
  )

  #
  dge <- edgeR::DGEList(counts = I2)
  expect_equal(
    object = which_k_covered_each_sample(dge, k = 1),
    expected = which_k_covered_each_sample(I2, k = 1),
    info = paste(
      "`which_k_covered...` should give the same results on DGEList as on",
      "the matrix of counts within it."
    )
  )

  #
  eset <- Biobase::ExpressionSet(assayData = I2)
  expect_equal(
    object = which_k_covered_each_sample(eset, k = 1),
    expected = which_k_covered_each_sample(I2, k = 1),
    info = paste(
      "`which_k_covered...` should give the same results on ExpressionSet as",
      "on the matrix of intensities within it."
    )
  )
})

###############################################################################

test_that("which_k_covered_across_samples", {
  # Assumes: genes are indexed on rows, samples on columns

  # 2 x 2 identity matrix
  I2 <- diag(rep(1, 2))
  # 2 x 2 all-ones upper-triangular
  T2 <- matrix(c(1, 0, 1, 1), nrow = 2)

  # Validity checks are done in check-validity.which_kCovered

  # Return all rows where coverage >= k across the summed columns

  # I2 should return empty sequence for k > 1
  expect_equal(
    object = which_k_covered_across_samples(mat = I2, k = 1.1),
    expected = integer(0),
    info = "2x2 identity: all rows sum to 1 => no rows pass when k > 1"
  )

  # I2 should return empty sequence for k > 1, even when passed as a data.frame
  expect_equal(
    object = which_k_covered_across_samples(mat = data.frame(I2), k = 1.1),
    expected = integer(0),
    info = "2x2 identity as a data.frame"
  )

  # I2 should return c(1, 2) for k <= 1
  expect_equal(
    object = which_k_covered_across_samples(mat = I2, k = 1),
    expected = c(1, 2),
    info = "2x2 identity: all rows sum to 1 => all rows pass when k <= 1"
  )

  # T2 should return empty sequence for k > 2
  expect_equal(
    object = which_k_covered_across_samples(mat = T2, k = 3),
    expected = integer(0),
    info = "[[1 1], [0 1]]: all rows should fail when k > 2"
  )

  # T2 should return c(1) for k == 2
  expect_equal(
    object = which_k_covered_across_samples(mat = T2, k = 2),
    expected = 1,
    info = "[[1 1], [0 1]]: first row should pass when k == 2"
  )

  # T2 should return c(1, 2) for k <= 1
  expect_equal(
    object = which_k_covered_across_samples(mat = T2, k = 1),
    expected = 1:2,
    info = "[[1 1], [0 1]]: all rows should pass when k <= 1"
  )

  dge <- edgeR::DGEList(counts = I2)
  expect_equal(
    object = which_k_covered_across_samples(dge, k = 1),
    expected = which_k_covered_across_samples(I2, k = 1),
    info = paste(
      "which_k_covered... should give the same results on DGEList as on the",
      "matrix of counts within it."
    )
  )
})

###############################################################################
