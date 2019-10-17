###############################################################################

context("Tests for methods on `sig_test_summary` objects")

###############################################################################

test_that("export - invalid input", {
  expect_error(
    object = export("NOT A sig_test_summary"),
    info = "`reeq::export` is only defined for `sig_test_summary` objects"
  )
})

test_that("export - creates the expected files", {
  # TODO: explicit constructor for `sig_test_summary` objects

  # `export.sig_test_summary` requires that `p_threshold`, `sig_features` and
  #   `top_table` are present in the `sig_test_summary` that is passed in.

  input <- structure(
    list(
      # 4 significant digits of the p-threshold are present in the output
      #   "significant features" file
      p_threshold = 0.001111,
      top_table = data.frame(x = 1:3, y = 4:6, row.names = LETTERS[1:3]),
      sig_features = "A"
    ),
    class = "sig_test_summary"
  )

  out_dir <- file.path(tempdir(), "export1")
  dir.create(out_dir)
  top_table_file <- file.path(
    out_dir, "some_contrast.top_tags.tsv"
  )
  sig_features_file <- file.path(
    out_dir, "some_contrast.sig_features.p0.0011.tsv"
  )

  expect_silent(
    export(
      input,
      output_dir = out_dir,
      test_name = "some_contrast"
    )
  )

  expect_true(
    file.exists(top_table_file),
    info = "`export` should make a top-tags file"
  )
  expect_true(
    file.exists(sig_features_file),
    info = "`export` should make a significant-features file"
  )
  expect_equal(
    object = scan(sig_features_file, what = character()),
    expected = input$sig_features,
    info = "contents of significant-features file should match the input"
  )
})
