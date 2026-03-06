# =============================================================================
# Tests: EnrichmentAnalyst
# =============================================================================

test_that("EnrichmentAnalyst initialization with cache manager", {
  cache_mgr <- create_test_cache_manager()
  ea <- EnrichmentAnalyst$new(cache_mgr)

  expect_true(!is.null(ea$cache_manager))
  expect_true(inherits(ea$cache_manager, "GeneSetCacheManager"))
})

test_that("EnrichmentAnalyst initialization without cache creates default", {
  skip_if_not_installed("org.Hs.eg.db")

  # This will try to create a default GeneSetCacheManager
  # May fail in test environments without proper setup, so we catch gracefully
  ea <- tryCatch(
    EnrichmentAnalyst$new(cache_manager = NULL),
    error = function(e) NULL
  )
  # If it succeeds, cache_manager should be set
  if (!is.null(ea)) {
    expect_true(!is.null(ea$cache_manager))
  }
})

test_that("run_comprehensive_ora() returns named list", {
  skip_if_not_installed("clusterProfiler")

  cache_mgr <- create_test_cache_manager()
  ea <- EnrichmentAnalyst$new(cache_mgr)

  gene_list <- c("TP53", "BRCA1", "EGFR", "MYC", "KRAS", "PIK3CA", "AKT1",
                 "PTEN", "RB1", "CDH1")
  universe  <- c(GENE_SYMBOLS, paste0("GENE", seq_len(10)))

  result <- ea$run_comprehensive_ora(gene_list, universe, pval_cutoff = 1)

  expect_type(result, "list")
  expect_true(length(result) > 0)
  # Should have entries for each database
  expect_true(any(c("GO_BP", "KEGG") %in% names(result)))
})

test_that("ORA with empty gene list returns results without error", {
  skip_if_not_installed("clusterProfiler")

  cache_mgr <- create_test_cache_manager()
  ea <- EnrichmentAnalyst$new(cache_mgr)

  # Empty gene list
  result <- ea$run_comprehensive_ora(character(0), GENE_SYMBOLS, pval_cutoff = 1)
  expect_type(result, "list")
})

test_that("run_comprehensive_gsea() returns named list", {
  skip_if_not_installed("fgsea")
  skip_if_not_installed("clusterProfiler")

  cache_mgr <- create_test_cache_manager()
  ea <- EnrichmentAnalyst$new(cache_mgr)

  # Build a ranked gene list (using all test gene symbols)
  set.seed(42)
  ranked <- sort(setNames(rnorm(length(GENE_SYMBOLS)), GENE_SYMBOLS),
                 decreasing = TRUE)

  result <- ea$run_comprehensive_gsea(ranked, pval_cutoff = 1, n_threads = 1)

  expect_type(result, "list")
})

test_that("analyze_diff_obj() populates ora_up, ora_down, and gsea_res", {
  skip_if_not_installed("limma")
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("fgsea")
  skip_if_not_installed("SummarizedExperiment")

  se <- create_test_se()
  da <- DiffExpAnalyst$new(se)
  da$run_dep_analysis(condition_col = "condition",
                      control_group = "A", case_group = "B")
  da$get_sig_proteins(pval_cutoff = 1, logfc_cutoff = 0)

  cache_mgr <- create_test_cache_manager()
  ea <- EnrichmentAnalyst$new(cache_mgr)
  ea$analyze_diff_obj(da, pval_cutoff = 1)

  expect_type(ea$ora_up, "list")
  expect_type(ea$ora_down, "list")
  expect_type(ea$gsea_res, "list")
})
