if (!requireNamespace("renv", quietly = TRUE)) stop("renv is required")
renv::restore(prompt = FALSE)
message("Phase 0 private renv library is ready.")
