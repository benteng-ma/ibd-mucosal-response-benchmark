source("R/audit_subject_overlap.R")
x <- audit_subject_overlap()
stopifnot(any(grepl("REUPLOAD", x$relationship_class)), any(grepl("SUBSET", x$relationship_class)))
