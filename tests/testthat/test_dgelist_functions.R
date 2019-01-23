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
