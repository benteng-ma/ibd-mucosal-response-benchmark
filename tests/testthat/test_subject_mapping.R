test_that("subject and sample manifests preserve explicit identifiers", {
  skip_if_not(file.exists(root_file("metadata", "subject_manifest.csv")), "Individual-level Phase 0 manifests are intentionally excluded from the public payload")
  subject <- readr::read_csv(root_file("metadata", "subject_manifest.csv"), show_col_types = FALSE)
  sample <- readr::read_csv(root_file("metadata", "sample_manifest.csv"), show_col_types = FALSE)
  expect_equal(anyDuplicated(paste(subject$dataset_id, subject$subject_id)), 0)
  expect_false(any(is.na(subject$subject_id) | subject$subject_id == ""))
  expect_false(any(is.na(sample$sample_id) | sample$sample_id == ""))
  expect_true(any(grepl("UNLINKED::", subject$subject_id, fixed = TRUE)))
})
