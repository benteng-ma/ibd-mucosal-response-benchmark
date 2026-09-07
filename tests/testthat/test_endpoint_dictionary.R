test_that("endpoint dictionary remains explicit and unique", {
  x <- readr::read_csv(root_file("metadata", "endpoint_dictionary.csv"), show_col_types = FALSE)
  expect_equal(nrow(x), 15)
  expect_equal(anyDuplicated(x$endpoint_family), 0)
  expect_true(all(c("RELAPSE_AFTER_WITHDRAWAL", "FLARE_ON_CONTINUED_THERAPY",
                    "MUCOSAL_HEALING", "INDUCTION_CLINICAL_RESPONSE") %in% x$endpoint_family))
})
