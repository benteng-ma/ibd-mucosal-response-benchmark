query_arrayexpress <- function(accession = NULL, path = "metadata/dataset_manifest.csv") {
  x <- readr::read_csv(path, show_col_types = FALSE)
  x <- x[x$repository == "BioStudies/ArrayExpress", ]
  if (is.null(accession)) x else x[x$accession %in% accession, ]
}
