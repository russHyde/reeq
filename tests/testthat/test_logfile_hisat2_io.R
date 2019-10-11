###############################################################################

context("Tests for the hisat2-logfile parsing functions")

###############################################################################

test_that("import data from a hisat2 logfile", {
  # hard-coded test file
  file <- "logfile_hisat2/paired_end_alignment.bam"

  expect_error(
    object = import_hisat2_summary(c(file, file)),
    info = "input to import_hisat2_summary should be a single string"
  )

  expect_error(
    object = import_hisat2_summary("Not a file"),
    info = "input to import_hisat2_summary should be a filepath"
  )

  expected <- tibble::tibble(
    rp_input = 14593447,
    rp_zero = 418319,
    rp_concordant_once = 13006296,
    rp_concordant_multi = 1145116,
    rp_discordant_once = 23716,
    unpaired = 836638,
    unpaired_unaligned = 434111,
    unpaired_once = 344675,
    unpaired_multi = 57852,
    align_rate = 0.9851
  )

  expect_equivalent(
    object = as.data.frame(import_hisat2_summary(file)),
    expected = as.data.frame(expected),
    info = "import hisat2 summary statistics from a file"
  )
})

###############################################################################

test_that("parse numeric fields from a hisat2 logfile's text", {
  expect_error(
    object = parse_hisat2_summary(),
    info = "Input to parse_hisat2_summary should be a string"
  )

  expect_error(
    object = parse_hisat2_summary(c("Two", "Strings")),
    info = "Input to parse_hisat2_summary should be a single string"
  )

  # nolint start

  example <-
    "HISAT2 summary stats:
	Total pairs: 114671
		Aligned concordantly or discordantly 0 time: 8256 (7.20%)
		Aligned concordantly 1 time: 90414 (78.85%)
		Aligned concordantly >1 times: 15517 (13.53%)
		Aligned discordantly 1 time: 484 (0.42%)
	Total unpaired reads: 16512
		Aligned 0 time: 12573 (76.14%)
		Aligned 1 time: 3078 (18.64%)
		Aligned >1 times: 861 (5.21%)
	Overall alignment rate: 94.52%"

  # nolint end

  object <- parse_hisat2_summary(example)

  expected <- tibble::tibble(
    rp_input = 114671,
    rp_zero = 8256,
    rp_concordant_once = 90414,
    rp_concordant_multi = 15517,
    rp_discordant_once = 484,
    unpaired = 16512,
    unpaired_unaligned = 12573,
    unpaired_once = 3078,
    unpaired_multi = 861,
    align_rate = 0.9452
  )

  expect_equal(
    object = as.data.frame(object),
    expected = as.data.frame(expected),
    info = "Single hisat2 text-parsing example (hacked into data.frame to get
    the test to work)"
  )
})
