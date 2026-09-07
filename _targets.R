library(targets)

tar_option_set(packages = c("readr", "dplyr", "yaml", "jsonlite"))

source("R/utils.R")
source("R/build_phase0_decision.R")

list(
  tar_target(phase0_config, yaml::read_yaml("config/phase0.yaml")),
  tar_target(dataset_manifest, readr::read_csv("metadata/dataset_manifest.csv", show_col_types = FALSE)),
  tar_target(trial_manifest, readr::read_csv("metadata/trial_manifest.csv", show_col_types = FALSE)),
  tar_target(task_manifest, readr::read_csv("metadata/task_manifest.csv", show_col_types = FALSE)),
  tar_target(signature_inventory, readr::read_csv("literature/signature_inventory.csv", show_col_types = FALSE)),
  tar_target(evaluation_grid, readr::read_csv("metadata/evaluation_grid.csv", show_col_types = FALSE)),
  tar_target(endpoint_mapping, readr::read_csv("metadata/endpoint_mapping.csv", show_col_types = FALSE)),
  tar_target(cohort_summary, readr::read_csv("metadata/cohort_summary.csv", show_col_types = FALSE)),
  tar_target(
    phase0_decision,
    build_phase0_decision(dataset_manifest, trial_manifest, task_manifest,
                          signature_inventory, evaluation_grid, endpoint_mapping,
                          cohort_summary)
  ),
  tar_target(
    phase0_decision_json,
    {
      jsonlite::write_json(phase0_decision, "reports/phase0_decision.json",
                           auto_unbox = TRUE, pretty = TRUE)
      "reports/phase0_decision.json"
    },
    format = "file"
  )
)
