test_that("evaluation grid is frozen, lineage-aware and performance-free", {
  x <- readr::read_csv(root_file("metadata", "evaluation_grid.csv"), show_col_types = FALSE)
  expect_equal(nrow(x), 11 * 25)
  forbidden <- c("auroc", "auc", "auprc", "p_value", "performance", "prediction", "score")
  expect_false(any(tolower(names(x)) %in% forbidden))
  expect_gte(length(unique(x$trial_id[x$transportability_eligible %in% TRUE])), 4)
  expect_true(any(x$evaluation_role == "SAME_THERAPY_SAME_ENDPOINT_EXTERNAL"))
  expect_true(any(x$evaluation_role == "CROSS_THERAPY_TRANSPORT_STRESS_TEST"))
  leaked <- x$primary_validation_eligible %in% TRUE & x$lineage_class != "INDEPENDENT_TRIAL_OR_COHORT"
  expect_false(any(leaked))
})
