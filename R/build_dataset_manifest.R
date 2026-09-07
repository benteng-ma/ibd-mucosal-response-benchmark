build_dataset_manifest <- function(path = "metadata/dataset_manifest.csv") {
  x <- readr::read_csv(path, show_col_types = FALSE)
  stopifnot(!anyDuplicated(x$dataset_id), all(c("accession", "trial_name", "platform") %in% names(x)))
  x
}
