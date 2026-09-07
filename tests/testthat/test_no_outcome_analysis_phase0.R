test_that("Phase 0 source tree contains no outcome-analysis implementation", {
  files <- c(list.files(root_file("R"), pattern = "[.]R$", full.names = TRUE),
             list.files(root_file("scripts"), pattern = "[.](R|py)$", full.names = TRUE))
  files <- files[!grepl("phase1", basename(files), ignore.case = TRUE)]
  text <- paste(vapply(files, function(path) paste(readLines(path, warn = FALSE), collapse = "\n"), character(1)), collapse = "\n")
  forbidden_calls <- c("DESeq2::", "limma::", "edgeR::", "ComBat(", "combat_seq(",
                       "randomForest(", "xgboost(", "glmnet(", "roc(", "t.test(", "wilcox.test(")
  expect_false(any(vapply(forbidden_calls, grepl, logical(1), x = text, fixed = TRUE)))
  phase0 <- yaml::read_yaml(root_file("config", "phase0.yaml"))
  expect_false(phase0$outcome_analysis_allowed)
  expect_false(phase0$performance_analysis_allowed)
  expect_false(dir.exists(root_file("results", "differential_expression")))
})
