###############################################################################

context("Tests for coverage functions in `reeq`")

###############################################################################

get_i2 <- function() {
  # 2 x 2 identity matrix
  diag(rep(1, 2))
}

get_t2 <- function() {
  # 2 x 2 all-ones upper-triangular
  matrix(c(1, 0, 1, 1), nrow = 2)
}

### ======================================================================= ###
#   Functions for QC-filtering on coverage
### ======================================================================= ###

test_that("check_validity_k_covered", {
  # `mat` should be a matrix or a data.frame
  expect_error(
    object = check_validity_k_covered(mat = "Not a matrix"),
    info = "`mat` should be a `matrix` in `which_k_covered_*`"
  )

  # If `mat` is a data.frame, it should have all-numeric columns
  expect_error(
    object = check_validity_k_covered(
      mat = data.frame(a = 1:2, b = letters[1:2]),
      k = 1
    ),
    info = "If `mat` is a `data.frame`, then it should be all-`numeric`"
  )

  # `k` should be a number
  expect_error(
    object = check_validity_k_covered(
      mat = get_i2(),
      k = "Not a number"
    ),
    info = "k should be a number in which_k_covered_*"
  )

  expect_error(
    object = check_validity_k_covered(
      mat = get_i2(),
      k = 1:2
    ),
    info = "k should be a single number in which_k_covered_*"
  )
})

###############################################################################

test_that("which_k_covered_each_sample", {
  # Assumes: genes are indexed on rows, samples on columns

  # get_i2 returns the 2 x 2 identity matrix

  # get_t2() returns the 2 x 2 all-ones upper-triangular matrix

  # Validity checks are done in .check_validity_k_covered

  # Return all rows where coverage >= k in every column

  # i2 should return empty sequence for k > 0
  expect_equal(
    object = which_k_covered_each_sample(mat = get_i2(), k = 1),
    expected = integer(0),
    info = paste(
      "2x2 identity: all rows contain some zeros => no rows pass when k > 0"
    )
  )

  # i2 should return empty sequence for k > 0, even when passed as a data.frame
  expect_equal(
    object = which_k_covered_each_sample(mat = data.frame(get_i2()), k = 1),
    expected = integer(0),
    info = "2x2 identity as data.frame"
  )

  # i2 should return c(1, 2) for k <= 0
  expect_equal(
    object = which_k_covered_each_sample(mat = get_i2(), k = 0),
    expected = c(1, 2),
    info = paste(
      "2x2 identity: all rows contain some zeros => all rows pass when k <= 0"
    )
  )

  expect_equal(
    object = which_k_covered_each_sample(mat = get_t2(), k = 2),
    expected = integer(0),
    info = paste(
      "[[1 1], [0 1]]: no rows have coverage >= k in all columns when k > 1"
    )
  )

  expect_equal(
    object = which_k_covered_each_sample(mat = get_t2(), k = 1),
    expected = 1,
    info = paste(
      "[[1 1], [0 1]]: only the first row has coverage >= k in all columns",
      "when k = 1"
    )
  )

  expect_equal(
    object = which_k_covered_each_sample(mat = get_t2(), k = 0),
    expected = c(1, 2),
    info = "[[1 1], [0 1]]: all rows have coverage >= k when k <= 0"
  )

  #
  dge <- edgeR::DGEList(counts = get_i2())
  expect_equal(
    object = which_k_covered_each_sample(dge, k = 1),
    expected = which_k_covered_each_sample(get_i2(), k = 1),
    info = paste(
      "`which_k_covered...` should give the same results on DGEList as on",
      "the matrix of counts within it."
    )
  )

  #
  eset <- Biobase::ExpressionSet(assayData = get_i2())
  expect_equal(
    object = which_k_covered_each_sample(eset, k = 1),
    expected = which_k_covered_each_sample(get_i2(), k = 1),
    info = paste(
      "`which_k_covered...` should give the same results on ExpressionSet as",
      "on the matrix of intensities within it."
    )
  )
})

###############################################################################

test_that("which_k_covered_across_samples", {
  # Assumes: genes are indexed on rows, samples on columns

  # get_i2() returns the 2 x 2 identity matrix

  # get_t2() returns the 2 x 2 all-ones upper-triangular

  # Validity checks are done in check-validity.which_kCovered

  # Return all rows where coverage >= k across the summed columns

  expect_equal(
    object = which_k_covered_across_samples(mat = get_i2(), k = 1.1),
    expected = integer(0),
    info = "2x2 identity: no rows have coverage-sum >= k when k > 1"
  )

  expect_equal(
    object = which_k_covered_across_samples(
      mat = data.frame(get_i2()), k = 1.1
    ),
    expected = integer(0),
    info = "2x2 identity as a data.frame"
  )

  expect_equal(
    object = which_k_covered_across_samples(mat = get_i2(), k = 1),
    expected = c(1, 2),
    info = "2x2 identity: all rows have coverage-sum >= k when k <= 1"
  )

  expect_equal(
    object = which_k_covered_across_samples(mat = get_t2(), k = 3),
    expected = integer(0),
    info = "[[1 1], [0 1]]: no rows have coverage-sum >= k when k > 2"
  )

  expect_equal(
    object = which_k_covered_across_samples(mat = get_t2(), k = 2),
    expected = 1,
    info = paste(
      "[[1 1], [0 1]]: only the first row has coverage-sum >= k when k == 2"
    )
  )

  expect_equal(
    object = which_k_covered_across_samples(mat = get_t2(), k = 1),
    expected = 1:2,
    info = "[[1 1], [0 1]]: all rows have coverage-sum >= k when k <= 1"
  )

  dge <- edgeR::DGEList(counts = get_i2())
  expect_equal(
    object = which_k_covered_across_samples(dge, k = 1),
    expected = which_k_covered_across_samples(get_i2(), k = 1),
    info = paste(
      "which_k_covered... should give the same results on DGEList as on the",
      "matrix of counts within it."
    )
  )
})

###############################################################################
