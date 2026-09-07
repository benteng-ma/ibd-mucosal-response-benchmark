parse_trial_metadata <- function(path = "metadata/trial_manifest.csv") {
  x <- readr::read_csv(path, show_col_types = FALSE)
  stopifnot(!anyDuplicated(x$trial_id))
  x
}
