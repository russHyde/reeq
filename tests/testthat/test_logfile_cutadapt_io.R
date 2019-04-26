###############################################################################

context("Tests for cutadapt logfile importing in `reeq`")

###############################################################################

# cutadapt parser tests

###############################################################################

test_that("cutadapt summary extractor", {
  prefix <- "This is cutadapt 1.13 with Python 3.6.7
Blah blah blah

=== Summary ===
"

  summary_section <- "Total read pairs processed:         14,605,217
  Read 1 with adapter:                 482,313 (3.3%)
  Read 2 with adapter:               1,012,028 (6.9%)
Pairs that were too short:               8,529 (0.1%)
Pairs with too many N:                   3,241 (0.0%)
Pairs written (passing filters):    14,593,447 (99.9%)

Total basepairs processed: 2,205,292,673 bp
  Read 1: 1,103,006,677 bp
  Read 2: 1,102,285,996 bp
Total written (filtered):  2,191,574,660 bp (99.4%)
  Read 1: 1,093,901,386 bp
  Read 2: 1,097,673,274 bp"

  suffix <- "

=== First read: Adapter 1 ===

Sequence: AAAAAAAAAAAAA..; Type: regular 3'; Length: 76; Trimmed: 173335 times.

Blah blah blah
"

  text <- paste0(prefix, summary_section, suffix, sep = "\n")

  expect_equal(
    extract_cutadapt_summary(text),
    summary_section,
    info = "Extract the summary section from the text in a cutadapt log"
  )
})


###############################################################################

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
    expected = tibble::tibble(
      rp_input = 114728, r1_with_adapter = 2923, r2_with_adapter = 3833,
      rp_too_short = 57, rp_output = 114671,
      bp_input = 17254556, bp_input_r1 = 8628890, bp_input_r2 = 8625666,
      bp_output = 17217852, bp_output_r1 = 8612728, bp_output_r2 = 8605124
    ),
    info = "Single cutadapt-summary parsing example"
  )
})
