PHASE0_ROOT <- normalizePath(testthat::test_path("..", ".."), winslash = "/", mustWork = TRUE)
root_file <- function(...) file.path(PHASE0_ROOT, ...)
