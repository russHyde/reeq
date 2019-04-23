###############################################################################

context("Tests for cutadapt logfile importing in `reeq`")

###############################################################################

# cutadapt parser tests

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
    expected = tibble::data_frame(
      field = character(0),
      value = character(0)
    ),
    info = paste(
      "parse_colon_separated_lines on an empty vector returns an empty,",
      "2-column data_frame"
    )
  )

  expect_equal(
    object = parse_colon_separated_lines(
      c("abc:123", "def:456", "ghi:some_value")
    ),
    expected = tibble::data_frame(
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
    expected = tibble::data_frame(
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
    expected = tibble::data_frame(
      field = c("spacey fieldname", "field"),
      value = c("123", "spacey value")
    ),
    info = paste(
      "parse_colon_separated_lines does not strip field/value-embedded",
      "whitespace"
    )
  )
})

test_that("cutadapt summary-text parser", {
  expect_error(
    object = parse_cutadapt_summary(),
    info = "parse_cutadapt_summary should be provided some text as argument"
  )

  example <-
    "Total read pairs processed:            114,728
  Read 1 with adapter:                   2,923 (2.5%)
  Read 2 with adapter:                   3,833 (3.3%)
Pairs that were too short:                  57 (0.0%)

Pairs written (passing filters):       114,671 (100.0%)
Total basepairs processed:    17,254,556 bp
  Read 1:     8,628,890 bp
  Read 2:     8,625,666 bp
Total written (filtered):     17,217,852 bp (99.8%)
  Read 1:     8,612,728 bp
  Read 2:     8,605,124 bp"

  expect_equal(
    object = parse_cutadapt_summary(example),
    expected = tibble::data_frame(
      rp_input = 114728, r1_with_adapter = 2923, r2_with_adapter = 3833,
      rp_too_short = 57, rp_output = 114671,
      bp_input = 17254556, bp_input_r1 = 8628890, bp_input_r2 = 8625666,
      bp_output = 17217852, bp_output_r1 = 8612728, bp_output_r2 = 8605124
    ),
    info = "Single cutadapt-summary parsing example"
  )
})
