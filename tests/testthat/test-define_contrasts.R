###############################################################################

context("tests for define_contrasts()")

###############################################################################

get_design1 <- function() {
  n_samples <- 10
  n_design_coefs <- 5

  samples <- paste0("s", seq(n_samples))
  design_coefs <- paste0("d", seq(n_design_coefs))

  matrix(
    rbinom(n = n_samples * n_design_coefs, size = 1, prob = 0.25),
    nrow = n_samples,
    dimnames = list(samples, design_coefs)
  )
}

test_that("define_contrasts: valid input", {
  # Invariants for define_contrasts:
  # - colnames of output should match .names
  # - rownames of output should match colnames of .design
  # - a numeric value in a .contrasts string should give a constant column
  # - output should be a matrix

  # User specifies:
  # - contrasts as a vector of strings
  # - names of those contrasts as a vector of strings

  # TODO:
  # - user can specify a design as a dataframe
  # - user can specify contrasts as a named vector (leaving .names NULL)

  design <- get_design1()
  design_coefs <- colnames(design)

  contrast_df <- tibble::tribble(
    ~name, ~definition,
    "zero", "0",
    "one", "1",
    "average", "(d1 + d2 + d3 + d4 + d5) / 5",
    "delta_2_1", "d2 - d1",
    "same_as_zero", "d1 - d1"
  )
  contrast_matrix <- as.matrix(data.frame(
    "zero" = 0,
    "one" = 1,
    "average" = 1 / 5,
    delta_2_1 = c(-1, 1, 0, 0, 0),
    "same_as_zero" = 0,
    row.names = design_coefs
  ))

  expect_equal(
    object = define_contrasts(contrast_df$definition, contrast_df$name, design),
    expected = contrast_matrix,
    info = paste(
      "define_contrasts() results should: be a matrix with dimnames",
      "(design-coefs, contrast-coefs)"
    )
  )
})

test_that("define_contrasts: invalid input", {
  # - .contrasts should be a vector of strings
  # - .names should be a vector of strings
  # - .names and .contrasts should be the same length
  # - all words in a .contrasts string should be a design column
  # - .design should be a matrix

  # TODO:
  # - .design could be a numeric-only data.frame

  design <- get_design1()

  expect_error(
    define_contrasts(1, "contrast1", design),
    info = ".contrasts should be a vector of strings"
  )

  expect_error(
    define_contrasts("d2 - d1", 3, design),
    info = ".names should be a vector of strings"
  )

  expect_error(
    define_contrasts("d2 - d1", c("a", "b"), design),
    info = ".names and .contrasts should be the same length"
  )

  expect_error(
    define_contrasts(c("0", "d2 - d1"), "a", design),
    info = ".names and .contrasts should be the same length"
  )

  expect_error(
    define_contrasts(c("not - a + valid - 2*contrast"), "a", design),
    info = "all words in the contrast string should be in the design colnames"
  )

  expect_error(
    define_contrasts("d2 - d1", "a", as.data.frame(design)),
    info = ".design should be a matrix"
  )
})
