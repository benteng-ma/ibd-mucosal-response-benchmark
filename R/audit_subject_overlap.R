audit_subject_overlap <- function(path = "metadata/data_lineage_matrix.csv") {
  x <- readr::read_csv(path, show_col_types = FALSE)
  stopifnot(all(c("relationship_class", "phase0_action") %in% names(x)))
  x
}
