test_that("dataset manifest has stable unique keys and required provenance", {
  x <- readr::read_csv(root_file("metadata", "dataset_manifest.csv"), show_col_types = FALSE)
  required <- c("dataset_id", "repository", "accession", "trial_name", "therapy_class",
                "sample_compartment", "assay", "platform", "limitations")
  expect_true(all(required %in% names(x)))
  expect_equal(anyDuplicated(x$dataset_id), 0)
  expect_gte(nrow(x), 20)
  expect_false(any(is.na(x$dataset_id) | x$dataset_id == ""))
})
