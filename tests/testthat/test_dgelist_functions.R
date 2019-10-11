###############################################################################

context("Tests for `DGEList` manipulation functions")

###############################################################################

test_that("filter_by_read_count: Valid input", {
  expect_equal(
    object = filter_by_read_count(get_dge1()),
    expected = get_dge1(),
    info = paste(
      "with no other args set, `filter_by_read_count(dge)` should not modify",
      "`dge`"
    )
  )

  # Filter by counts: all samples should have at least `k` counts for each
  # passing feature
  expect_equal(
    object = filter_by_read_count(get_dge1(), threshold = 3),
    expected = dge(get_dge1()$counts[3:5, ]),
    info = paste(
      "[[1:5], [6:10]] : Keep only rows where every entry is >= 3 (ie, 3:5)"
    )
  )

  # Filter by `cpm`: all samples should have a `cpm` of at least `k` for each
  # passing feature (in the input dataset)
  # - first col of dge1 sums to 15, second to 40
  # - let's keep features if all samples have cpm > (1.5 / 15) * 10^6 = 10^5
  # - the only entry in dge1 with cpm < 10^5 is [1, 1]
  expect_equal(
    object = filter_by_read_count(
      get_dge1(),
      threshold = 1.0e5, count_type = "cpm"
    ),
    expected = dge(get_dge1()$counts[2:5, ]),
    info = paste(
      "[[1:5], [6:10]] : Keep only rows where cpm >= 10^5 (ie, 2:5)"
    )
  )

  # Filter by counts in a fraction of samples: for any passing feature, at
  # least X samples must have coverage at level `k`
  expect_equal(
    object = filter_by_read_count(
      get_dge1(),
      threshold = 9, fraction_of_samples = 0.4
    ),
    expected = dge(get_dge1()$counts[4:5, ]),
    info = paste(
      "[[1:5], [6:10]] : Keep only rows where >= 40% samples have >= 9",
      "coverage (ie, 4:5)"
    )
  )

  # Filter by average count: any passing feature should have an average count
  # of at least `k`
  expect_equal(
    object = filter_by_read_count(
      get_dge1(),
      threshold = 4.5, use_sample_average = TRUE
    ),
    expected = dge(get_dge1()$counts[2:5, ]),
    info = paste(
      "[[1:5], [6:10]] : Keep only rows where average count >= 4.5 (ie, 2:5)"
    )
  )

  # Ensure as `complete` a DGEList is returned as is passed in:
  # - if a `genes` entry is present in the input, then it should be in the
  # output..
  expect_equal(
    object = filter_by_read_count(
      get_dge1(add_genes = TRUE),
      threshold = 4
    ),
    expected = dge(
      get_dge1()$counts[4:5, ],
      add_genes = TRUE
    ),
    info = paste(
      "Output DGEList should contain all entries that input DGEList contains",
      "(ie, counts, genes, samples)"
    )
  )
})

test_that("filter_by_read_count: Invalid input", {
  expect_error(
    filter_by_read_count("NOT A DGEList"),
    info = "`filter_by_read_count` expects a DGEList as input"
  )

  expect_error(
    filter_by_read_count(get_dge1(), threshold = "Not a number"),
    info = "`threshold` argument should be a single real number"
  )

  expect_error(
    filter_by_read_count(get_dge1(), count_type = "invalid_type"),
    info = "`count_type` argument should be one of `counts` or `cpm`"
  )

  expect_error(
    filter_by_read_count(get_dge1(), fraction_of_samples = c(0.5, 0.6)),
    info = "`fraction_of_samples` should be a single value in [0, 1]"
  )

  expect_error(
    filter_by_read_count(
      get_dge1(),
      threshold = 1, use_sample_average = 1
    ),
    info = "`use_sample_average` should be Boolean"
  )
})


###############################################################################

# Checklist:
# - [+] no feature-data in DGEList
# - [+] feature_id column present in both DGEList and feature_id
# - [+] the DGEList does not have a feature_id column (use rownames of the
#   former)
# - [] feature-data of DGEList contains more features than the new feature-data
# - [] a tibble can be passed as the new feature-data
# -

test_that("A genes data-frame can be added, and appended-to", {
  # Features are "g1" to "g5"; columns are "feature_id" and "length"
  dge_without_genes <- get_dge1(add_genes = FALSE)
  dge_with_genes <- get_dge1(add_genes = TRUE)

  # Feature data can be added to a DGEList that has no `genes` dataframe

  features <- .df(
    feature_id = paste0("g", 5:1),
    some_annotation = 200
  )

  expect_is(
    append_feature_annotations(dge_without_genes, features, "feature_id"),
    "DGEList"
  )

  expect_equal(
    object = append_feature_annotations(
      dge_without_genes, features, "feature_id"
    )[["genes"]],
    expected = .df(
      feature_id = paste0("g", 1:5),
      some_annotation = 200,
      row.names = paste0("g", 1:5)
    ),
    info = "no feature-data in DGEList"
  )

  # Feature data can be appended to an existing `genes` entry:
  features1 <- .df(
    my_id = paste0("g", c(3, 4, 2, 5, 1)),
    my_annotation = 200
  )
  expect_equal(
    object = append_feature_annotations(
      dge_with_genes, features1, "my_id"
    )[["genes"]],
    expected = .df(
      feature_id = paste0("g", 1:5),
      length = 10,
      my_id = paste0("g", 1:5),
      my_annotation = 200,
      row.names = paste0("g", 1:5)
    ),
    info = paste(
      "feature IDs are matched on DGE rownames if feature-ID column is",
      "absent from DGE$genes"
    )
  )

  # nolint start
  dge_with_mismatching_rownames_and_feature_ids <- get_dge1(add_genes = TRUE)
  dge_with_mismatching_rownames_and_feature_ids$genes <- dplyr::mutate(
    dge_with_mismatching_rownames_and_feature_ids$genes,
    feature_id = letters[1:5]
  )
  # nolint end

  features2 <- .df(
    feature_id = letters[5:1],
    annotation = 1234
  )

  expect_equal(
    object = append_feature_annotations(
      dge_with_mismatching_rownames_and_feature_ids, features2, "feature_id"
    )[["genes"]],
    expected = .df(
      feature_id = letters[1:5],
      length = 10,
      annotation = 1234,
      row.names = paste0("g", 1:5)
    ),
    info = paste(
      "feature-IDs are matched by common feature_id column if it is present",
      "in both DGEList$genes and feature data (preferentialy over rownames of",
      "$genes)"
    )
  )
})

###############################################################################

# Checklist:
# - [+] returns a DGEList
# - [+] error if the dimensions of the offset don't match the counts matrix
# - [+] can add a dim-matched offset to a DGEList
# - [+] error if not a DGEList
# - [+] error if not a numeric matrix
test_that("An offset matrix can be added to a DGEList", {
  dge <- get_dge1()
  new_offset <- edgeR::getCounts(dge) * -1

  expect_is(
    set_offset(dge, new_offset),
    "DGEList"
  )

  expect_equal(
    edgeR::getCounts(set_offset(dge, new_offset)),
    edgeR::getCounts(dge),
    info = "set_offset should not affect the counts matrix"
  )

  expect_error(
    set_offset(dge, cbind(new_offset, new_offset)),
    info = "offset matrix should have same number of cols as dge"
  )

  expect_error(
    set_offset(dge, rbind(new_offset, new_offset)),
    info = "offset matrix should have same number of rows as dge"
  )

  expect_equal(
    edgeR::getOffset(set_offset(dge, new_offset)),
    new_offset,
    info = "can add a dimension-matched offset matrix to a DGEList"
  )

  expect_error(
    set_offset(list(), new_offset),
    info = "input to set_offset should be a DGEList"
  )

  expect_error(
    set_offset(dge, "not a numeric matrix"),
    info = "`offset` should be a numeric matrix"
  )
})

###############################################################################

test_that("user can run cqn on a dgelist", {
  dge <- get_dge2()

  expect_error(
    cqn_dgelist(x = "Not a DGEList"),
    info = "Input to `cqn_dgelist` should be a DGEList"
  )

  expect_error(
    cqn_dgelist(x = dge, length_column = "not a column"),
    info = "length_column should be defined in x$genes"
  )

  expect_error(
    cqn_dgelist(x = dge, gc_column = "not a column"),
    info = "gc_column should be defined in x$genes"
  )

  expect_error(
    cqn_dgelist(x = dge, lib_size_column = "not a column"),
    info = "lib_size_column should be defined in x$genes"
  )

  mock_cqn <- mockery::mock(NULL)
  mockery::stub(cqn_dgelist, "cqn::cqn", mock_cqn)

  result <- cqn_dgelist(dge)
  expect_null(result)
  mockery::expect_called(mock_cqn, 1)
})
