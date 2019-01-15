###############################################################################

context("Tests for importing `featureCounts` datasets")

###############################################################################

test_that("read_single_feature_counts_file", {
  expect_equal(
    read_single_feature_counts_file("feature_counts/short_1.tsv"),
    tibble::tibble(
      feature_id = c("ENSG00000223972", "ENSG00000227232"),
      length = c(1735L, 1351L),
      "some/file/name.bam" = c(0L, 2L)
    ),
    info = "read a single valid short-format featureCounts file"
  )

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

test_that("read_feature_counts", {

})

###############################################################################
