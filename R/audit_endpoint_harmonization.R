audit_endpoint_harmonization <- function(dictionary = "metadata/endpoint_dictionary.csv",
                                         mapping = "metadata/endpoint_mapping.csv") {
  list(
    dictionary = readr::read_csv(dictionary, show_col_types = FALSE),
    mapping = readr::read_csv(mapping, show_col_types = FALSE)
  )
}
