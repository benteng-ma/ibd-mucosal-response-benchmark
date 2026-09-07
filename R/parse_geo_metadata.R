parse_geo_metadata <- function(accession, directory = "data/raw/geo_metadata") {
  path <- file.path(directory, paste0(accession, "_samples.csv"))
  if (!file.exists(path)) stop("Missing parsed GEO metadata: ", path)
  readr::read_csv(path, show_col_types = FALSE)
}
