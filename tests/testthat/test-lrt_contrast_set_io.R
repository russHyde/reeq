context(
  paste(
    "Tests for importing a hierarchical tree of contrast sets for use in",
    "likelihood ratio tests"
  )
)

# ? should parse_lrt_table behave differently when the LRT-sets are going to
# be used in nested vs unnested analysis of contrast-sets

test_that(
  "user can parse an LRT-set containing a single contrast-collection",
  {

    # ... one contrast collection with one contrast within it
    table <- tibble::tibble(
      lrt_name = "all", parent = NA, level = 1,
      contrast_set = "treatment_vs_control"
    )
    expected <- tibble::tibble(
      lrt_name = "all", parent = NA, level = 1,
      contrast_set = list("treatment_vs_control"), n_siblings = 1
    )
    expect_equal_tbl(
      parse_lrt_table(table),
      expected,
      info = "a single contrast-collection containing a single contrast"
    )

    # ... one contrast collection with multiple contrasts within it
    table <- tibble::tibble(
      lrt_name = "all", parent = NA, level = 1,
      contrast_set = "treatment1;treatment2"
    )
    expected <- tibble::tibble(
      lrt_name = "all", parent = NA, level = 1,
      contrast_set = list(c("treatment1", "treatment2")), n_siblings = 1
    )
    expect_equal_tbl(
      parse_lrt_table(table),
      expected,
      info = paste(
        "a single contrast-collection containing multiple ';'-separated",
        "contrast names"
      )
    )
  }
)

test_that(
  "user can parse an LRT-set containing multiple contrast-collections",
  {
    table <- tibble::tribble(
      ~lrt_name, ~parent, ~level, ~contrast_set,
      "drug1", NA, 1, "treatment1",
      "drug2", NA, 1, "treatment2"
    )
    expected <- tibble::tibble(
      lrt_name = c("drug1", "drug2"),
      parent = c(NA, NA),
      level = c(1, 1),
      contrast_set = list("treatment1", "treatment2"),
      n_siblings = c(2, 2)
    )
    expect_equal_tbl(
      parse_lrt_table(table),
      expected,
      info = paste(
        "two non-overlapping contrast collections, each containing a single",
        "contrast"
      )
    )

    # Multiple contrast collections, level-2 contrast sets are nested inside the
    # level-1 contrast set
    table <- tibble::tribble(
      ~lrt_name, ~parent, ~level, ~contrast_set,
      "all", NA, 1, "treatment1;treatment2",
      "drug1", "all", 2, "treatment1",
      "drug2", "all", 2, "treatment2"
    )

    expected <- tibble::tibble(
      lrt_name = c("all", "drug1", "drug2"),
      parent = c(NA, "all", "all"),
      level = c(1, 2, 2),
      contrast_set = list(
        c("treatment1", "treatment2"), "treatment1", "treatment2"
      ),
      n_siblings = c(1, 2, 2)
    )
    expect_equal_tbl(
      parse_lrt_table(table),
      expected,
      info = paste(
        "nested contrast collections; level-2 nodes have a common level-1",
        "parent"
      )
    )

    # Multiple contrast collections, level-2 contrast sets have distinct
    # parent nodes in level-1
    table <- tibble::tribble(
      ~lrt_name, ~parent, ~level, ~contrast_set,
      "p1", NA, 1, "cA;cB",
      "p2", NA, 1, "cC;cD",
      "c1", "p1", 2, "cA",
      "c2", "p2", 2, "cC"
    )
    expected <- tibble::tibble(
      lrt_name = c("p1", "p2", "c1", "c2"),
      parent = c(NA, NA, "p1", "p2"),
      level = c(1, 1, 2, 2),
      contrast_set = list(c("cA", "cB"), c("cC", "cD"), "cA", "cC"),
      n_siblings = c(2, 2, 1, 1)
    )
    expect_equal_tbl(
      parse_lrt_table(table),
      expected,
      info = paste(
        "nodes in lower levels that do not share parents are not siblings"
      )
    )
  }
)

test_that("User can parse an LRT table with contrast-collections as a list", {
  table <- tibble::tibble(
    lrt_name = c("drug1", "drug2"),
    parent = NA,
    level = 1,
    contrast_set = list("treatment1", "treatment2")
  )
  expected <- dplyr::mutate(table, n_siblings = 2)
  expect_equal_tbl(
    parse_lrt_table(table),
    expected,
    info = paste(
      "two contrast collections, passed as a list (with single entries)"
    )
  )

  table <- tibble::tibble(
    lrt_name = c("drug1", "drug2"),
    parent = NA,
    level = 1,
    contrast_set = list(
      c("treatment1", "treatment3"), c("treatment2", "treatment4")
    )
  )
  expected <- dplyr::mutate(table, n_siblings = 2)
  expect_equal_tbl(
    parse_lrt_table(table),
    expected,
    info = paste(
      "two contrast collections, passed as a list (with multiple entries)"
    )
  )
})

test_that("parse_lrt_table should be idempotent", {
  # That is, running parse_lrt_table on a table once should give the same
  # result as if you ran it multiple times

  table <- tibble::tibble(
    lrt_name = c("drug1", "drug2"),
    parent = NA,
    level = 1,
    contrast_set = list("treatment1", "treatment2")
  )
  expect_equal_tbl(
    parse_lrt_table(parse_lrt_table(table)),
    parse_lrt_table(table),
    info = "parse_lrt_table should be idempotent"
  )
})

test_that("parents come before children in the output dataframe", {
  table <- tibble::tribble(
    ~lrt_name, ~parent, ~level, ~contrast_set,
    "node2", "all", 2, "contrast2",
    "all", NA, 1, "contrast2;contrast1"
  )
  expected <- tibble::tibble(
    lrt_name = c("all", "node2"),
    parent = c(NA, "all"),
    level = c(1, 2),
    contrast_set = list(c("contrast2", "contrast1"), c("contrast2")),
    n_siblings = c(1, 1)
  )
  expect_equal_tbl(
    parse_lrt_table(table),
    expected,
    info = paste(
      "in the data-frame returned by parse_lrt_table, the rows are ordered by",
      "increasing 'level' so that parent-nodes occur before their children"
    )
  )
})

test_that("inappropriate args for parse_lrt_table", {
  table <- tibble::tibble(
    lrt_name = "all", parent = NA, level = 1,
    contrast_set = "treatment_vs_control"
  )
  for (i in seq_len(ncol(table))) {
    column_name <- colnames(table)[i]
    expect_error(
      parse_lrt_table(table[-(i)]),
      info = sprintf(
        "`x` should have column '%s' in parse_lrt_table", column_name
      )
    )
  }

  # contrast set should be a vector of semi-colon separated strings or a list
  # of character-vectors
  table <- tibble::tibble(
    lrt_name = "all", parent = NA, level = 1,
    contrast_set = 123
  )
  expect_error(
    parse_lrt_table(table),
    info = "contrast_set should be a character vector (or list thereof)"
  )

  table <- tibble::tibble(
    lrt_name = letters[1:3], parent = c(NA, "a", "a"), level = c(1, 2, 2),
    contrast_set = list(123, "x", factor("a"))
  )
  expect_error(
    parse_lrt_table(table),
    info = "if contrast_set is a list, each entry should be a character-vector"
  )
})
