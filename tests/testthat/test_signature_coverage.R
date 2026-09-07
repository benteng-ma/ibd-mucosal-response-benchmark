test_that("coverage audit is bounded and separates uncomputable mappings", {
  x <- readr::read_csv(root_file("metadata", "signature_coverage_matrix.csv"), show_col_types = FALSE)
  expect_equal(nrow(x), 11 * 22)
  expect_true(all(is.na(x$coverage_fraction) | (x$coverage_fraction >= 0 & x$coverage_fraction <= 1)))
  expect_true(any(x$mapping_status == "NOT_COMPUTABLE_NO_DIRECT_PLATFORM_CATALOG"))
  expect_true(any(x$primary_coverage_pass %in% TRUE))
})
