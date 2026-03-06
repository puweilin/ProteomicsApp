# =============================================================================
# Tests: ProteomicsGSVA
# =============================================================================

test_that("ProteomicsGSVA initialization validates input type", {
  expect_error(
    ProteomicsGSVA$new(data_manager = "not_a_manager"),
    "must be a ProteomicsDataManager"
  )
})

test_that("ProteomicsGSVA initialization sets empty fields", {
  mgr <- create_test_data_manager()
  gsva <- ProteomicsGSVA$new(mgr)

  expect_type(gsva$gsva_matrices, "list")
  expect_length(gsva$gsva_matrices, 0)
  expect_type(gsva$gsva_results, "list")
  expect_length(gsva$gsva_results, 0)
  expect_type(gsva$diff_pathways, "list")
  expect_length(gsva$diff_pathways, 0)
})

test_that("clean_pathway_names() works correctly", {
  mgr <- create_test_data_manager()
  gsva <- ProteomicsGSVA$new(mgr)

  clean_fn <- gsva$.__enclos_env__$private$clean_pathway_names
  result <- clean_fn("GOBP_CELL_CYCLE")
  expect_equal(result, "Cell Cycle")

  result2 <- clean_fn("KEGG_APOPTOSIS")
  expect_equal(result2, "Apoptosis")
})

test_that("clean_pathway_df() adds Pathway_ID column", {
  mgr <- create_test_data_manager()
  gsva <- ProteomicsGSVA$new(mgr)

  test_df <- data.frame(Pathway = "GOBP_CELL_CYCLE", logFC = 1.5,
                        stringsAsFactors = FALSE)

  clean_fn <- gsva$.__enclos_env__$private$clean_pathway_df
  result <- clean_fn(test_df)

  expect_true("Pathway_ID" %in% colnames(result))
  expect_equal(result$Pathway_ID, "GOBP_CELL_CYCLE")
  expect_equal(result$Pathway, "Cell Cycle")
})

test_that("clear_gsva_cache() removes cached files", {
  mgr <- create_test_data_manager()
  gsva <- ProteomicsGSVA$new(mgr)

  # Write a dummy cache file
  dummy_file <- file.path(gsva$gsva_cache_dir, "test_cache.rds")
  saveRDS("dummy", dummy_file)
  expect_true(file.exists(dummy_file))

  gsva$clear_gsva_cache()

  remaining <- list.files(gsva$gsva_cache_dir)
  expect_length(remaining, 0)
})

test_that("GSVA cache dir is created on init", {
  mgr <- create_test_data_manager()
  gsva <- ProteomicsGSVA$new(mgr)

  expect_true(dir.exists(gsva$gsva_cache_dir))
})

test_that("plot_pathway_volcano() errors without data", {
  mgr <- create_test_data_manager()
  gsva <- ProteomicsGSVA$new(mgr)

  expect_error(gsva$plot_pathway_volcano(db = "GOBP"),
               "No differential results")
})

test_that("plot_pathway_heatmap() errors without data", {
  mgr <- create_test_data_manager()
  gsva <- ProteomicsGSVA$new(mgr)

  expect_error(gsva$plot_pathway_heatmap(db = "GOBP"),
               "No GSVA matrix")
})

test_that("get_gsva_matrix() errors without data", {
  mgr <- create_test_data_manager()
  gsva <- ProteomicsGSVA$new(mgr)

  expect_error(gsva$get_gsva_matrix("GOBP"),
               "No GSVA matrix")
})
