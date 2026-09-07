test_that("known reuploads, subsets and same-trial tissues are encoded", {
  x <- readr::read_csv(root_file("metadata", "data_lineage_matrix.csv"), show_col_types = FALSE)
  expect_true(any(x$relationship_class == "SAME_BIOLOGICAL_SAMPLE_REUPLOAD"))
  expect_true(any(x$relationship_class == "SAME_ACCESSION_SUBSET"))
  expect_true(any(x$dataset_a == "GSE112366_UNITI2_ILEUM" & x$dataset_b == "GSE207022_UNITI2_RECTUM"))
  comparison <- jsonlite::read_json(root_file("results", "lineage", "act1_overlap_aggregate.json"))
  expect_equal(comparison$n_exact_uncompressed_cel_hash_matches, 23)
  expect_true(comparison$all_files_found)
})
