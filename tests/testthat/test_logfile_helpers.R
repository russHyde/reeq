###############################################################################

context("Tests for logfile-helper functions")

###############################################################################

test_that("spread_and_rename_numeric_fields", {

  expect_error(
    spread_and_rename_numeric_fields(
      x = tibble::tibble(field = letters[1:3], value = 1:3),
      fieldnames = tibble::tibble(
        wrong = 1:3, columns = 4:6
      )
    ),
    info = "fieldnames should have columns (expected,output)"
  )

  expect_error(
    spread_and_rename_numeric_fields(
      x = tibble::tibble(wrong = 1:3, columns = 4:6),
      fieldnames = tibble::tibble(
        expected = "Some Field", output = "some_field"
      )
    ),
    info = "x should have columns (field,value)"
  )

  expect_error(
    spread_and_rename_numeric_fields(
      x = tibble::tibble(
        field = rep("some_field", 2), value = 1
      ),
      fieldnames = tibble::tibble(
        expected = "some_field", output = "some_field"
      )
    ),
    info = "There should be no duplicates in x$field"
  )
})

###############################################################################

test_that("parse_colon_separated_lines", {
  expect_error(
    object = parse_colon_separated_lines(),
    info = "`x` should be defined in parse_colon_separated_lines"
  )

  expect_error(
    object = parse_colon_separated_lines("abc 123"),
    info = "each entry of `x` should contain a colon"
  )

  expect_equal(
    object = parse_colon_separated_lines(c()),
    expected = tibble::tibble(
      field = character(0),
      value = character(0)
    ),
    info = paste(
      "parse_colon_separated_lines on an empty vector returns an empty,",
      "2-column tibble"
    )
  )

  expect_equal(
    object = parse_colon_separated_lines(
      c("abc:123", "def:456", "ghi:some_value")
    ),
    expected = tibble::tibble(
      field = c("abc", "def", "ghi"),
      value = c("123", "456", "some_value")
    ),
    info = "parse_colon_separated_lines without whitespace"
  )

  expect_equal(
    object = parse_colon_separated_lines(
      c(
        "  abc:   123", "\t\tdef:456", "ghi:   some_value\t\t"
      )
    ),
    expected = tibble::tibble(
      field = c("abc", "def", "ghi"),
      value = c("123", "456", "some_value")
    ),
    info = paste(
      "parse_colon_separated_lines strips trailing/leading whites"
    )
  )

  expect_equal(
    object = parse_colon_separated_lines(
      c("spacey fieldname: 123", "field: spacey value")
    ),
    expected = tibble::tibble(
      field = c("spacey fieldname", "field"),
      value = c("123", "spacey value")
    ),
    info = paste(
      "parse_colon_separated_lines does not strip field/value-embedded",
      "whitespace"
    )
  )
})

###############################################################################

test_that("parse_numeric_fields from a colon-separated logfile text", {
  expect_error(
    object = parse_numeric_fields(),
    info = paste(
      "`parse_numeric_fields` requires a character string (x), and a",
      "data-frame (fieldnames)"
    )
  )

  expect_equal(
    object = parse_numeric_fields(
      x = paste(
        "   my_field_name   : 123,456.45",
        "drop_the_trailing_percent : 12345 (0.01%)",
        "omit_any_nonnumeric_lines : not_a_number",
        "drop_the_trailing_bp : 987 bp   ",
        "drop_the_percent_sign    : 10%",
        "\tOverall alignment rate: 98.51%",
        sep = "\n"
      )
    ),
    expected = tibble::tribble(
      ~field, ~value,
      "my_field_name", 123456.45,
      "drop_the_trailing_percent", 12345,
      "drop_the_trailing_bp", 987,
      "drop_the_percent_sign", 10,
      "Overall alignment rate", 98.51
    ),
    info = paste(
      "extract numeric fields from the summary section of a cutadapt log"
    )
  )

  # expect_error(
  #  object = parse_numeric_fields(
  #    x = "observed_field1 : 123\nobserved_field2 : 456",
  #    fieldnames = tibble::tibble(
  #      expected = c("observed_field1"),
  #      output = c("observed_field1")
  #    )
  #  ),
  #  info = paste(
  #    "all numeric fields in the cutadapt log-summary should have a",
  #    "corresponding entry in fieldnames (so they can be reformatted)"
  #  )
  # )
})
##############################################################################
