###############################################################################

context("Tests for utility functions")

###############################################################################

test_that("replace with", {
  expect_equal(
    replace_with(c(), letters, LETTERS),
    c(),
    info = "empty vector has nothing to replace"
  )

  expect_equal(
    replace_with("absent_from_letters", letters, LETTERS),
    "absent_from_letters",
    info = "don't replace a value if it is absent from search-list"
  )

  expect_equal(
    replace_with(c("a", "b", "a"), letters, LETTERS),
    c("A", "B", "A"),
    info = "replace_with where all of `x` are in the search_list"
  )

  expect_equal(
    replace_with(c("a", "A", NA, "d", "jkl"), letters, LETTERS),
    c("A", "A", NA, "D", "jkl"),
    info = "replace_with where some values are absent from the search_list"
  )

  expect_error(
    replace_with("some_value", letters[1:3], LETTERS[1:4]),
    info = "search_list and return_list should be the same length"
  )

  expect_error(
    replace_with(c("a", "b", "c", "D"), letters, LETTERS, strict = TRUE),
    info = paste(
      "when strict=TRUE, every element in `x` should be in the `search_list`"
    )
  )
})
