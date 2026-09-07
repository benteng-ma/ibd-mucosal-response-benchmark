test_that("Public repository is portable and separated from submission material", {
  expect_true(file.exists(root_file("README.md")))
  expect_true(file.exists(root_file("LICENSE")))
  expect_false(dir.exists(root_file("submission")))
  expect_false(file.exists(root_file("metadata", "manuscript_identity_and_disclosures.json")))
})
