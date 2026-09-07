dir.create("results/logs", recursive = TRUE, showWarnings = FALSE)
required <- c("targets", "testthat", "yaml", "jsonlite", "readr", "dplyr")
status <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
writeLines(c(capture.output(sessionInfo()), "", paste(names(status), status, sep = "=")),
           "results/logs/session_info.txt")
if (!all(status)) stop("Missing packages: ", paste(names(status)[!status], collapse = ", "))
