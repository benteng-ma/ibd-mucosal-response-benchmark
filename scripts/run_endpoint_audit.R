source("R/audit_endpoint_harmonization.R")
x <- audit_endpoint_harmonization()
stopifnot(nrow(x$dictionary) == 15, nrow(x$mapping) >= 20)
