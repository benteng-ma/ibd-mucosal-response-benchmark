test_that("machine-readable decision passes schema and Phase 0 ceiling", {
  x <- jsonlite::read_json(root_file("reports", "phase0_decision.json"), simplifyVector = TRUE)
  expect_identical(x$decision, "GO_CROSS_THERAPY_ENDPOINT_BENCHMARK")
  expect_identical(x$audit_cutoff, "2026-08-31")
  expect_gte(x$evidence_summary$reconstructable_main_grid_signatures, 8)
  expect_gte(x$evidence_summary$independent_trials_in_transportability_grid, 4)
  expect_true(grepl("no predictive performance", x$highest_justified_evidence, ignore.case = TRUE))
  expect_length(x$unmet_criteria, 0)
})
