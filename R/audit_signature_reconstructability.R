audit_signature_reconstructability <- function(inventory = "literature/signature_inventory.csv",
                                                genes = "literature/signature_gene_table.csv") {
  inv <- readr::read_csv(inventory, show_col_types = FALSE)
  tbl <- readr::read_csv(genes, show_col_types = FALSE)
  list(inventory = inv, genes = tbl)
}
