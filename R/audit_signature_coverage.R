audit_signature_coverage <- function(path = "metadata/signature_coverage_matrix.csv") {
  x <- readr::read_csv(path, show_col_types = FALSE)
  stopifnot(all(is.na(x$coverage_fraction) | (x$coverage_fraction >= 0 & x$coverage_fraction <= 1)))
  x
}
