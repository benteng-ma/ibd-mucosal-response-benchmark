test_that("main-grid signatures are finite, sourced and directional", {
  inv <- readr::read_csv(root_file("literature", "signature_inventory.csv"), show_col_types = FALSE)
  genes <- readr::read_csv(root_file("literature", "signature_gene_table.csv"), show_col_types = FALSE)
  main <- inv[toupper(inv$main_grid_eligible) == "TRUE", ]
  expect_gte(nrow(main), 8)
  expect_gte(length(unique(main$paper_id)), 5)
  expect_true(all(main$reconstruction_status == "COMPLETE"))
  expect_true(all(main$direction_available %in% TRUE))
  expect_true(all(main$signature_id %in% genes$signature_id))
  expect_false(any(is.na(genes$gene_symbol) | genes$gene_symbol == ""))
})
