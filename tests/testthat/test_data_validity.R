###############################################################################

context("Tests for data-validity functions in `reeq`")

###############################################################################

test_that("is_single_string - helper function", {
  expect_true(
    object = is_single_string("abc"),
    info = "a single string"
  )
  expect_false(
    object = is_single_string(1),
    info = "a number, not a string"
  )
  expect_false(
    object = is_single_string(c("abc", "def")),
    info = "more than one string"
  )
})

###############################################################################

test_that("is_nonempty_df", {
  expect_error(
    is_nonempty_df(),
    info = "no input to is_nonempty_df"
  )

  expect_equal(
    object = is_nonempty_df("Not a data.frame"),
    expected = FALSE,
    info = "String input to is_nonempty_df"
  )

  expect_equal(
    object = is_nonempty_df(data.frame()),
    expected = FALSE,
    info = "data.frrame with no cols and no rows is empty"
  )

  expect_equal(
    object = is_nonempty_df(data.frame(a = character(0))),
    expected = FALSE,
    info = "data.frame with no rows is empty"
  )

  expect_equal(
    object = is_nonempty_df(data.frame(row.names = letters[1:2])),
    expected = FALSE,
    info = "data.frame with no cols is empty"
  )

  expect_equal(
    object = is_nonempty_df(data.frame(a = 1:3)),
    expected = TRUE,
    info = "data.frame with some entries is nonempty"
  )
})

###############################################################################

test_that("is_nonempty_list", {
  expect_error(
    object = is_nonempty_list(),
    info = "no input to is_nonempty_list"
  )

  expect_equal(
    object = is_nonempty_list("Not a list"),
    expected = FALSE,
    info = "Non-list input to is_nonempty_list"
  )

  expect_equal(
    object = is_nonempty_list(list()),
    expected = FALSE,
    info = "Empty list input to is_nonempty_list"
  )

  expect_equal(
    object = is_nonempty_list(list(1, 2, 3)),
    expected = TRUE,
    info = "Unnamed list with 3 entries"
  )

  expect_true(
    object = is_nonempty_list(list(NULL)),
    info = ".x is a list of NULL in input to is_nonempty_list"
  )

  expect_equal(
    object = is_nonempty_list(data.frame(a = 1:3)),
    expected = TRUE,
    info = paste(
      "`data.frame`s are `list`s, so should be valid input to",
      "`is_nonempty_list`"
    )
  )
})

###############################################################################

test_that("is_valid_fcounts_df", {
  expect_true(
    object = is_valid_fcounts_df(
      data.frame(
        feature_id = c("ENSG01234567890"),
        Length = 1234,
        some_sample = 9876,
        stringsAsFactors = FALSE
      )
    ),
    info = paste(
      "`df` with `feature_id`, `Length` and at-least-one sample is a valid",
      "feature-counts `df`"
    )
  )

  expect_false(
    object = is_valid_fcounts_df(
      data.frame(
        Geneid = "abc",
        Length = 1234,
        some_sample = 9876,
        stringsAsFactors = FALSE
      )
    ),
    info = paste(
      "feature-counts `data.frame`s should use `feature_id`, not `Geneid`",
      "to indicate feature"
    )
  )
})
