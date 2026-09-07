assert_phase0_boundary <- function(config) {
  forbidden <- c(
    "outcome_analysis_allowed", "differential_expression_allowed",
    "signature_scoring_allowed", "performance_analysis_allowed",
    "model_training_allowed", "cross_cohort_combat_allowed"
  )
  values <- unlist(config[forbidden], use.names = TRUE)
  if (any(values %in% TRUE)) stop("Phase 0 boundary violation")
  invisible(TRUE)
}

