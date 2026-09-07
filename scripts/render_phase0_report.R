required <- c(
  "reports/novelty_audit.md", "reports/dataset_audit.md",
  "reports/data_lineage_audit.md", "reports/endpoint_harmonization_audit.md",
  "reports/signature_reconstructability_audit.md", "reports/evidence_ceiling.md",
  "reports/risk_register.md", "reports/phase0_decision.md",
  "reports/phase0_decision.json"
)
stopifnot(all(file.exists(required)))
message("All Phase 0 report sources are present; no Phase 1 analysis executed.")
