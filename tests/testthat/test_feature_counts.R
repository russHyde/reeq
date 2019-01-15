###############################################################################

context("Tests for feature-count manipulating functions")

###############################################################################

test_that("feature_counts_to_dgelist: invalid input", {
  fcounts_df1 <- .df(feature_id = letters[1:3], length = 1:3, id1 = 11:13)
  samp_df1 <- .df(id = paste0("id", 1:3))

  # fcounts_df should be:
  # - a non-empty data.frame
  # - with feature_id ...
  # - and length as the first two column names,
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
      fcounts_df = .df(length = 1:10, id1 = 11:20),
      sample_df = samp_df1,
      id_column = "id"
    ),
    info = "fcounts_df should have feature_id as first column"
  )
  expect_error(
    object = feature_counts_to_dgelist(
      fcounts_df = .df(feature_id = letters[1:10], id1 = 11:20),
      sample_df = samp_df1,
      id_column = "id"
    ),
    info = "fcounts_df should have length as second column"
  )
  expect_error(
    object = feature_counts_to_dgelist(
      fcounts_df = .df(
        feature_id = letters[1:10],
        length = letters[1:10],
        id1 = 11:20
      ),
      sample_df = samp_df1,
      id_column = "id"
    ),
    info = "fcounts_df::length should be numerical"
  )
  expect_error(
    object = feature_counts_to_dgelist(
      fcounts_df = .df(feature_id = letters[1:10], length = 1:10),
      sample_df = samp_df1,
      id_column = "id"
    ),
    info = "Should be at least one column of counts in fcounts_df"
  )

  expect_error(
    object = feature_counts_to_dgelist(
      fcounts_df = .df(
        feature_id = letters[1:10],
        length = 1:10,
        id1 = c("NOT", "A", "COUNT")
      ),
      sample_df = samp_df1,
      id_column = "id"
    ),
    info = "All columns other than feature_id and length should be counts"
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
  fcounts_df1 <- .df(feature_id = letters[1:3], length = 1:3, id1 = 11:13)
  fcounts_df2 <- .df(
    feature_id = letters[1:3],
    length = 1:3,
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
      feature_id = letters[1:3],
      length = 1:3
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
        feature_id = letters[1:3],
        length = 1:3
      )
    )
  )
})

###############################################################################

test_that("cbind_feature_counts", {

  # - Input should be a list of data.frames, each containing three columns
  # (feature_id, length, <count.column>) in that order
  # - count.column may or may not have a sample-specific id as header
  # - and the user may want to rename the cols anyway (eg, if they are
  # long-winded filenames)

  # Invalid input
  # - Missing input
  expect_error(
    object = cbind_feature_counts(),
    info = "Missing input to cbind_feature_counts"
  )
  # - Not a list
  expect_error(
    object = cbind_feature_counts("Not a list"),
    info = "Input to cbind_feature_counts should be a list"
  )
  # - An empty list
  expect_error(
    object = cbind_feature_counts(list()),
    info = "Empty list input to cbind_feature_counts"
  )
  # - A list with non-data.frame contents
  expect_error(
    object = cbind_feature_counts(list("Not a data.frame")),
    info = "Input should be a list of data.frames"
  )
  # - A list with data.frames that don't have feature_id column
  expect_error(
    object = cbind_feature_counts(list(.df(
      A = 1,
      length = 2,
      Count = 3
    ))),
    info = "1st column should be 'feature_id' in each data.frame"
  )
  # - A list with data.frames that don't have length column
  expect_error(
    object = cbind_feature_counts(list(.df(
      feature_id = 1,
      B = 2,
      Count = 3
    ))),
    info = "2nd column should be 'length' in each data.frame"
  )
  # - A list with data.frames that don't have counts column (too few columns)
  expect_error(
    object = cbind_feature_counts(list(.df(
      feature_id = 1,
      length = 2
    ))),
    info = "Missing count column - too few columns"
  )
  # - A list with data.frames that don't have counts column (non-numeric col)
  expect_error(
    object = cbind_feature_counts(list(.df(
      feature_id = 1,
      length = 2,
      not.county = "A"
    ))),
    info = "Missing count column - non-numeric column"
  )
  # - A list with data.frames that have more than three cols
  expect_error(
    object = cbind_feature_counts(list(.df(
      feature_id = 1,
      length = 2,
      county = 123,
      extra = "ABC"
    ))),
    info = "Too many columns in input"
  )
  # - names_as_colnames = TRUE, but no names in the input list
  expect_error(
    object = cbind_feature_counts(
      list(.df(
        feature_id = 1,
        length = 2,
        count.col = 345
      )),
      names_as_colnames = TRUE
    ),
    info = paste(
      "If names_as_colnames = TRUE, there should be names in the count_list"
    )
  )

  # Valid input
  # - Function uses `merge` semantics over the feature_id and length
  # so there doesn't need to be exactly the same features in each input
  # data.frame
  # - names_as_colnames determines what the column names should be in the
  # output

  # Single data.frame input
  x <- .df(
    feature_id = c(1, 2, 3),
    length = c(10, 30, 30),
    count.col = c(123, 234, 345)
  )
  expect_equal(
    object = cbind_feature_counts(
      list(x),
      names_as_colnames = FALSE
    ),
    expected = x,
    info = "Single data.frame input, names_as_colnames = FALSE"
  )
  expect_equal(
    object = cbind_feature_counts(
      list(a = x),
      names_as_colnames = TRUE
    ),
    expected = setNames(x, c("feature_id", "length", "a")),
    info = "Single data.frame input, names_as_colnames = TRUE"
  )

  # Two data.frames, same features
  x <- .df(
    feature_id = c(1, 2),
    length = c(10, 30),
    count.col1 = c(123, 234)
  )
  y <- .df(
    feature_id = c(1, 2),
    length = c(10, 30),
    count.col2 = c(111, 222)
  )
  expect_equal(
    object = cbind_feature_counts(
      list(x, y),
      names_as_colnames = FALSE
    ),
    expected = .df(
      feature_id = c(1, 2),
      length = c(10, 30),
      count.col1 = c(123, 234),
      count.col2 = c(111, 222)
    ),
    info = "Two data.frame input with names.as.conames = FALSE"
  )

  # Two data.frames, names_as_colnames, count.colname is the same in input
  expect_equal(
    object = cbind_feature_counts(
      Map(
        function(df) {
          setNames(df, c("feature_id", "length", "count"))
        },
        list(a = x, b = y)
      ),
      names_as_colnames = TRUE
    ),
    expected = .df(
      feature_id = c(1, 2),
      length = c(10, 30),
      a = c(123, 234),
      b = c(111, 222)
    ),
    info = "Two data.frames, names_as_colnames = TRUE, same count colname"
  )

  # read_tsv returns a tbl_df
  # But, if T is a tbl_df, and the third column stores a numeric entry,
  #   then is.numeric(T[, 3]) is FALSE
  x <- dplyr::tbl_df(
    .df(feature_id = c(1, 2), length = c(2, 3), count.col = c(10, 20))
  )
  expect_equal(
    object = cbind_feature_counts(list(x = x), names_as_colnames = TRUE),
    expected = data.frame(
      feature_id = c(1, 2),
      length = c(2, 3),
      x = c(10, 20)
    ),
    info = "tbl_df as input - this should fail if is.numeric(x[, 3]) is tested"
  )
})

###############################################################################
