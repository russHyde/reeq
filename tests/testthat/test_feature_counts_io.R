###############################################################################

context("Tests for importing `featureCounts` datasets")

###############################################################################

test_that("read_single_feature_counts_file - from file", {
  expect_equal(
    read_single_feature_counts_file("feature_counts/short_1.tsv"),
    tibble::tibble(
      feature_id = c("ENSG00000223972", "ENSG00000227232"),
      length = c(1735L, 1351L),
      "some/file/name.bam" = c(0L, 2L)
    ),
    info = "read a single valid short-format featureCounts file"
  )
})

test_that("read_single_feature_counts_file - mocking read_tsv", {
  m <- mockery::mock(
    tibble::tibble(
      Geneid = c("abc", "def", "ghi"),
      Length = c(123L, 234L, 345L),
      counts = c(1L, 0L, 10L),
      non_numeric = c("shouldn't", "be", "possible")
    )
  )
  with_mock(
    read_tsv = m, {
      expect_error(
        read_single_feature_counts_file("malformed_contents_file"),
        info = "All non-ID columns should be numeric"
      )
    },
    .env = "readr"
  )
})

###############################################################################

test_that("read_feature_counts - matching features", {
  file_contents <- tibble::tibble(
    Geneid = c("abc", "def"),
    Length = c(123L, 234L),
    some_file_name = c(1L, 10L)
  )

  m <- mockery::mock(file_contents, cycle = TRUE)

  with_mock(
    read_tsv = m, {
      expect_equivalent(
        object = read_feature_counts(
          files = c("some_file", "some_other_file"),
          sample_ids = c("sample1", "sample2")
        ),
        expected = file_contents %>%
          dplyr::transmute(
            feature_id = Geneid, length = Length, sample1 = some_file_name,
            sample2 = some_file_name
          ),
        info = paste(
          "import and combine several short-format featureCounts files that",
          "have the same features"
        )
      )
    },
    .env = "readr"
  )
})

###############################################################################
