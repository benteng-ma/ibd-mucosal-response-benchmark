query_geo <- function(accession = NULL, path = "metadata/dataset_manifest.csv") {
  x <- readr::read_csv(path, show_col_types = FALSE)
  if (is.null(accession)) x[x$repository == "GEO", ] else x[x$accession %in% accession, ]
}
