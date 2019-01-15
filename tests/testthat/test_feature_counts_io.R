###############################################################################

context("Tests for importing `featureCounts` datasets")

###############################################################################

setup_fcounts_contents <- function() {
  tibble::tibble(
    Geneid = c("abc", "def"),
    Length = c(123L, 234L),
    some_file_name = c(1L, 10L)
  )
}

###############################################################################

test_that("read_single_feature_counts_file - from file", {
  expect_equal(
    read_single_feature_counts_file("feature_counts/short_1.tsv"),
    tibble::tibble(
      feature_id = c("ENSG00000223972", "ENSG00000227232"),
      length = c(1735L, 1351L),
      "some/file/name.bam" = c(0L, 2L)
    ),
    info = "read a single valid short-format `featureCounts` file"
  )

  expect_equal(
    read_single_feature_counts_file(
      "feature_counts/standard_1.tsv",
      col_types = "cccccii"
    ),
    tibble::tibble(
      feature_id = c("ENSG00000223972", "ENSG00000227232"),
      length = c(1735L, 1351L),
      "some/file/name.bam" = c(0L, 2L)
    ),
    info = "read a single valid standard-format `featureCounts` file"
  )
})

###############################################################################

test_that("read_single_feature_counts_file - mocking read_tsv", {
  m <- mockery::mock(
    dplyr::mutate(
      setup_fcounts_contents(),
      non_numeric = c("shouldn't be", "possible")
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
  file_contents <- setup_fcounts_contents()

  m <- mockery::mock(file_contents, cycle = TRUE)

  with_mock(
    read_tsv = m, {
      expect_equivalent(
        object = read_feature_counts(
          files = c("file1", "file2"), sample_ids = c("sample1", "sample2")
        ),
        expected = dplyr::transmute(
          file_contents,
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

test_that("read_feature_counts - non-matching features", {
  file1 <- setup_fcounts_contents()
  file2 <- dplyr::mutate(file1, Geneid = c("different", "IDs"))
  file3 <- dplyr::mutate(file1, Length = c(10000L, 100000L))

  m1 <- mockery::mock(file1, file2)
  m2 <- mockery::mock(file1, file3)

  with_mock(
    read_tsv = m1, {
      expect_error(
        read_feature_counts(c("file1", "file2"), c("sample1", "sample2")),
        info = paste(
          "Gene IDs should match when collapsing several feature-counts files"
        )
      )
    },
    .env = "readr"
  )

  with_mock(
    read_tsv = m2, {
      expect_error(
        read_feature_counts(c("file1", "file3"), c("sample1", "sample3")),
        info = paste(
          "Gene Lengths should match when collapsing several feature-counts",
          "files"
        )
      )
    },
    .env = "readr"
  )
})
