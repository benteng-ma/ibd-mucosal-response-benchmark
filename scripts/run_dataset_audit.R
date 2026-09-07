source("R/build_dataset_manifest.R")
x <- build_dataset_manifest()
stopifnot(nrow(x) >= 20)
