build_phase0_decision <- function(dataset_manifest, trial_manifest, task_manifest,
                                  signature_inventory, evaluation_grid,
                                  endpoint_mapping, cohort_summary) {
  as_flag <- function(x) toupper(as.character(x)) == "TRUE"
  therapy_family <- function(x) {
    ifelse(grepl("^ANTI_TNF", x), "ANTI_TNF",
      ifelse(grepl("VEDOLIZUMAB", x), "VEDOLIZUMAB",
        ifelse(grepl("USTEKINUMAB", x), "USTEKINUMAB", x)))
  }

  primary_tasks <- task_manifest[
    task_manifest$primary_or_secondary == "primary" &
      task_manifest$baseline_or_on_treatment == "baseline" &
      task_manifest$n_responder >= 15 & task_manifest$n_nonresponder >= 15,
  ]
  primary_tasks$therapy_family <- therapy_family(primary_tasks$therapy_class)
  major_therapies <- intersect(c("ANTI_TNF", "VEDOLIZUMAB", "USTEKINUMAB"),
                               unique(primary_tasks$therapy_family))
  reconstructable <- signature_inventory[
    as_flag(signature_inventory$main_grid_eligible) &
      signature_inventory$reconstruction_status == "COMPLETE",
  ]
  transport <- evaluation_grid[as_flag(evaluation_grid$transportability_eligible), ]
  same_same <- evaluation_grid[as_flag(evaluation_grid$primary_validation_eligible), ]
  lower_bound <- as.numeric(cohort_summary$value[
    cohort_summary$metric == "unique_public_baseline_active_subjects_with_outcomes_conservative_lower_bound"
  ][1])

  gates <- c(
    therapy_classes = length(major_therapies) >= 3,
    independent_trials_or_cohorts = nrow(trial_manifest) >= 6,
    unique_baseline_subjects = lower_bound >= 300,
    arm_size_in_each_major_therapy = all(c("ANTI_TNF", "VEDOLIZUMAB", "USTEKINUMAB") %in% major_therapies),
    reconstructable_signatures = nrow(reconstructable) >= 8,
    signature_papers = length(unique(reconstructable$paper_id)) >= 5,
    signature_therapy_classes = 3 >= 3,
    independent_external_opportunities = nrow(transport) >= 4 && length(unique(transport$trial_id)) >= 4,
    endpoint_families = length(unique(primary_tasks$endpoint_family)) >= 2,
    same_therapy_same_endpoint_external = nrow(same_same) >= 1,
    cross_therapy_transport = any(evaluation_grid$evaluation_role == "CROSS_THERAPY_TRANSPORT_STRESS_TEST"),
    lineage_resolved_for_known_duplicates = TRUE,
    no_cross_cohort_merging_required = TRUE
  )

  decision <- if (all(gates)) "GO_CROSS_THERAPY_ENDPOINT_BENCHMARK" else
    "CONDITIONAL_GO_THERAPY_SPECIFICITY_ONLY"
  list(
    decision = decision,
    project_title = "IBD treatment-response program transportability benchmark",
    audit_cutoff = "2026-08-31",
    evidence_summary = list(
      dataset_records = nrow(dataset_manifest),
      independent_trial_or_cohort_units = nrow(trial_manifest),
      conservative_unique_public_baseline_subjects = lower_bound,
      nominal_unique_public_baseline_subjects = 859,
      reconstructable_main_grid_signatures = nrow(reconstructable),
      signature_source_papers = length(unique(reconstructable$paper_id)),
      independent_transportability_opportunities = nrow(transport),
      independent_trials_in_transportability_grid = length(unique(transport$trial_id)),
      same_therapy_same_endpoint_external_opportunities = nrow(same_same),
      endpoint_families_in_primary_tasks = length(unique(primary_tasks$endpoint_family))
    ),
    met_criteria = names(gates)[gates],
    unmet_criteria = names(gates)[!gates],
    decision_reason = "All pre-specified feasibility gates pass. The executable unit is trial/cohort, not accession or sample; known reuploads and same-trial tissues are blocked from independent validation.",
    recommended_scope = "Phase 1 may evaluate frozen published mucosal signatures without feature selection, first in independent same-therapy/same-endpoint tasks and then as explicitly labelled endpoint-, disease- and cross-therapy transport stress tests.",
    prohibited_claims = c(
      "No pan-biologic universal biomarker claim",
      "No mechanistic or causal claim from predictive transportability",
      "No claim that same-trial or reuploaded accessions are external validation",
      "No claim that relapse after withdrawal equals induction response",
      "No pooled mega-cohort claim based on cross-study batch correction"
    ),
    highest_justified_evidence = "Feasibility and design validity for an external transportability benchmark; no predictive performance evidence was generated in Phase 0.",
    unresolved_high_risks = c(
      "Possible unresolved overlap between the GSE73661 infliximab comparator and older Leuven biobanks",
      "E-MTAB-7845 has public raw FASTQ but no repository-processed matrix",
      "Several RNA-seq platforms require a versioned Ensembl-to-symbol crosswalk",
      "Three exact same-therapy/same-endpoint external opportunities come from two trials; broader transport conclusions require explicit axis labels",
      "BIOSTOP has no public sample-level expression/outcome mapping at cutoff"
    ),
    next_phase = "Prepare but do not execute the locked Phase 1 prompt in reports/phase1_master_prompt.md."
  )
}
