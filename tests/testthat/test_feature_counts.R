###############################################################################

context("Tests for feature-count manipulating functions")

###############################################################################
# helpers:

.df <- function(...) data.frame(..., stringsAsFactors = FALSE)

###############################################################################

test_that("feature_counts_to_dgelist: invalid input", {
  fcounts_df1 <- .df(Geneid = letters[1:3], Length = 1:3, id1 = 11:13)
  samp_df1 <- .df(id = paste0("id", 1:3))

  # fcounts_df should be:
  # - a non-empty data.frame
  # - with Geneid ...
  # - and Length as the first two column names,
  # - and numeric entries elsewhere
  expect_error(
    object = feature_counts_to_dgelist(
      fcounts_df = .df(),
      sample_df = samp_df1,
      id_column = "id"
    ),
    info = "fcounts_df should be a non-empty data.frame"
  )
  expect_error(
    object = feature_counts_to_dgelist(
      fcounts_df = .df(Length = 1:10, id1 = 11:20),
      sample_df = samp_df1,
      id_column = "id"
    ),
    info = "fcounts_df should have Geneid as first column"
  )
  expect_error(
    object = feature_counts_to_dgelist(
      fcounts_df = .df(Geneid = letters[1:10], id1 = 11:20),
      sample_df = samp_df1,
      id_column = "id"
    ),
    info = "fcounts_df should have Length as second column"
  )
  expect_error(
    object = feature_counts_to_dgelist(
      fcounts_df = .df(
        Geneid = letters[1:10],
        Length = letters[1:10],
        id1 = 11:20
      ),
      sample_df = samp_df1,
      id_column = "id"
    ),
    info = "fcounts_df::Length should be numerical"
  )
  expect_error(
    object = feature_counts_to_dgelist(
      fcounts_df = .df(Geneid = letters[1:10], Length = 1:10),
      sample_df = samp_df1,
      id_column = "id"
    ),
    info = "Should be at least one column of counts in fcounts_df"
  )

  expect_error(
    object = feature_counts_to_dgelist(
      fcounts_df = .df(
        Geneid = letters[1:10],
        Length = 1:10,
        id1 = c("NOT", "A", "COUNT")
      ),
      sample_df = samp_df1,
      id_column = "id"
    ),
    info = "All columns other than Geneid and Length should be counts"
  )
  # `sample_df` should be
  # - a `data.frame`
  # - with either sample.ids defined in the rownames ...
  # - ... or in a column that is indexed by the id_column
  # - and every one of the ids in the columns of fcounts_df should be present
  # as a sample.id in sample_df
  expect_error(
    object = feature_counts_to_dgelist(
      fcounts_df = fcounts_df1,
      sample_df = .df(),
      id_column = "id"
    ),
    info = "Sample.df should be a non-empty data.frame"
  )

  expect_error(
    object = feature_counts_to_dgelist(
      fcounts_df = fcounts_df1,
      sample_df = .df(id = "NOT_id1"),
      id_column = "id"
    ),
    info = "Each sample.id used in fcounts_df should be present in the
  id_column of sample_df"
  )

  expect_error(
    object = feature_counts_to_dgelist(
      fcounts_df = fcounts_df1,
      sample_df = .df(row.names = "NOT_id1", a = 1, b = "some_file.txt"),
      id_column = NULL
    ),
    info = "Each sample.id used in fcounts_df should be present in the
  row.names of sample_df, if id_column = NULL"
  )

  expect_error(
    object = feature_counts_to_dgelist(
      fcounts_df = fcounts_df1,
      sample_df = .df(id = "id1", a = 1, b = "some_file.txt"),
      id_column = "not-a-column"
    ),
    info = "If specified, id_column should be a column of sample_df"
  )

  expect_error(
    object = feature_counts_to_dgelist(
      fcounts_df = fcounts_df1,
      sample_df = .df(id = "id1", a = 1, b = "some_file.txt"),
      id_column = NULL
    ),
    info = paste(
      "If id_column is NULL, rownames of sample_df are used as sample-IDs;",
      "so the rownames can't be NULL"
    )
  )

  expect_error(
    object = feature_counts_to_dgelist(
      fcounts_df = fcounts_df1,
      sample_df = .df(id = rep("id1", 2), a = 1, b = "some_file.txt"),
      id_column = "id"
    ),
    info = paste(
      "Sample IDs shouldn't be repeated in feature_counts_to_dgelist"
    )
  )
})

###############################################################################

#' @importFrom   edgeR         DGEList
#' @importFrom   magrittr      %>%   set_rownames
#'
test_that("feature_counts_to_dgelist: valid input", {
  fcounts_df1 <- .df(Geneid = letters[1:3], Length = 1:3, id1 = 11:13)
  fcounts_df2 <- .df(
    Geneid = letters[1:3],
    Length = 1:3,
    id1 = 11:13,
    id2 = 20:22,
    id3 = 31:33
  )
  samp_df1 <- .df(id = paste0("id", 1:3))

  object <- feature_counts_to_dgelist(
    fcounts_df = fcounts_df1,
    sample_df = .df(id = "id1"),
    id_column = "id"
  )
  expect <- edgeR::DGEList(
    counts = matrix(
      11:13,
      nrow = 3,
      dimnames = list(letters[1:3], "id1")
    ),
    samples = .df(id = "id1"),
    genes = .df(
      row.names = letters[1:3],
      Geneid = letters[1:3],
      Length = 1:3
    )
  )
  expect_equal(
    object = object,
    expected = expect,
    info = "single entry"
  )

  expect_equal(
    object = feature_counts_to_dgelist(
      fcounts_df = fcounts_df2,
      sample_df = .df(id = paste0("id", 1:3)),
      id_column = "id"
    ),
    expected = edgeR::DGEList(
      counts = matrix(
        c(11:13, 20:22, 31:33),
        nrow = 3,
        dimnames = list(letters[1:3], paste0("id", 1:3))
      ),
      samples = magrittr::set_rownames(samp_df1, samp_df1$id),
      genes = .df(
        row.names = letters[1:3],
        Geneid = letters[1:3],
        Length = 1:3
      )
    )
  )
})

###############################################################################
