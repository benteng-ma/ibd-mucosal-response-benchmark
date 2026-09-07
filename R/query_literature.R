query_literature <- function(path = "literature/search_log.csv") {
  readr::read_csv(path, show_col_types = FALSE)
}
