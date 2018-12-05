###############################################################################

context("Tests for data-validity functions in `reeq`")

###############################################################################

test_that(".is_single_string - helper function", {
  expect_true(
    object = .is_single_string("abc"),
    info = "a single string"
  )
  expect_false(
    object = .is_single_string(1),
    info = "a number, not a string"
  )
  expect_false(
    object = .is_single_string(c("abc", "def")),
    info = "more than one string"
  )
})

###############################################################################
