source("R/query_literature.R")
log <- query_literature()
novelty <- readr::read_csv("literature/novelty_matrix.csv", show_col_types = FALSE)
stopifnot(nrow(log) == 57, !any(toupper(novelty$exact_duplicate) == "TRUE"))
