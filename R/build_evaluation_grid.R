build_evaluation_grid <- function(path = "metadata/evaluation_grid.csv") {
  x <- readr::read_csv(path, show_col_types = FALSE)
  forbidden <- c("auroc", "auc", "auprc", "p_value", "performance", "prediction")
  stopifnot(!any(tolower(names(x)) %in% forbidden))
  x
}
