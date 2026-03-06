# =============================================================================
# Tests: fgsea_to_gseaResult conversion utility
# =============================================================================

test_that("fgsea_to_gseaResult returns gseaResult S4 object", {
  skip_if_not_installed("fgsea")
  skip_if_not_installed("clusterProfiler")

  # Build mock fgsea result
  fgsea_res <- data.frame(
    pathway      = c("GO:0006915", "GO:0007049"),
    pval         = c(0.001, 0.05),
    padj         = c(0.01, 0.1),
    log2err      = c(0.5, 0.3),
    ES           = c(0.6, -0.4),
    NES          = c(1.8, -1.2),
    size         = c(50, 30),
    stringsAsFactors = FALSE
  )
  fgsea_res$leadingEdge <- list(c("TP53", "BCL2", "BAX"), c("MYC", "RB1"))

  gene_list <- sort(setNames(rnorm(100), paste0("GENE", 1:100)), decreasing = TRUE)
  gene_sets <- list(
    "GO:0006915" = c("TP53", "BCL2", "BAX", "CASP3"),
    "GO:0007049" = c("MYC", "RB1", "CDH1")
  )

  result <- fgsea_to_gseaResult(fgsea_res, gene_list, gene_sets)

  expect_true(inherits(result, "gseaResult"))
})

test_that("fgsea_to_gseaResult result slot has required columns", {
  skip_if_not_installed("fgsea")
  skip_if_not_installed("clusterProfiler")

  fgsea_res <- data.frame(
    pathway = c("PathA"),
    pval    = c(0.01),
    padj    = c(0.05),
    log2err = c(0.5),
    ES      = c(0.5),
    NES     = c(1.5),
    size    = c(20),
    stringsAsFactors = FALSE
  )
  fgsea_res$leadingEdge <- list(c("TP53", "EGFR"))

  gene_list <- sort(setNames(rnorm(50), paste0("G", 1:50)), decreasing = TRUE)
  gene_sets <- list(PathA = c("TP53", "EGFR", "MYC"))

  result <- fgsea_to_gseaResult(fgsea_res, gene_list, gene_sets)
  res_df <- result@result

  expect_true(all(c("ID", "Description", "NES", "pvalue", "p.adjust",
                     "core_enrichment") %in% colnames(res_df)))
})

test_that("fgsea_to_gseaResult .sign column reflects NES direction", {
  skip_if_not_installed("fgsea")
  skip_if_not_installed("clusterProfiler")

  fgsea_res <- data.frame(
    pathway = c("Up_path", "Down_path"),
    pval    = c(0.01, 0.02),
    padj    = c(0.05, 0.1),
    log2err = c(0.5, 0.3),
    ES      = c(0.6, -0.5),
    NES     = c(1.8, -1.5),
    size    = c(20, 15),
    stringsAsFactors = FALSE
  )
  fgsea_res$leadingEdge <- list(c("A", "B"), c("C", "D"))

  gene_list <- sort(setNames(rnorm(50), paste0("G", 1:50)), decreasing = TRUE)
  gene_sets <- list(Up_path = c("A", "B"), Down_path = c("C", "D"))

  result <- fgsea_to_gseaResult(fgsea_res, gene_list, gene_sets)
  res_df <- result@result

  expect_true(".sign" %in% colnames(res_df))
  expect_equal(res_df[res_df$ID == "Up_path", ".sign"], "activated")
  expect_equal(res_df[res_df$ID == "Down_path", ".sign"], "suppressed")
})

test_that("fgsea_to_gseaResult returns NULL for empty input", {
  skip_if_not_installed("fgsea")

  fgsea_res <- data.frame(
    pathway = character(0), pval = numeric(0), padj = numeric(0),
    log2err = numeric(0), ES = numeric(0), NES = numeric(0), size = integer(0)
  )
  fgsea_res$leadingEdge <- list()

  gene_list <- sort(setNames(rnorm(10), paste0("G", 1:10)), decreasing = TRUE)
  gene_sets <- list(PathA = c("G1", "G2"))

  result <- fgsea_to_gseaResult(fgsea_res, gene_list, gene_sets)
  expect_null(result)
})

test_that("fgsea_to_gseaResult term2name mapping works", {
  skip_if_not_installed("fgsea")
  skip_if_not_installed("clusterProfiler")

  fgsea_res <- data.frame(
    pathway = c("GO:0006915"),
    pval    = c(0.01),
    padj    = c(0.05),
    log2err = c(0.5),
    ES      = c(0.5),
    NES     = c(1.5),
    size    = c(20),
    stringsAsFactors = FALSE
  )
  fgsea_res$leadingEdge <- list(c("TP53"))

  term2name <- data.frame(
    term = "GO:0006915",
    name = "apoptotic process",
    stringsAsFactors = FALSE
  )

  gene_list <- sort(setNames(rnorm(50), paste0("G", 1:50)), decreasing = TRUE)
  gene_sets <- list("GO:0006915" = c("TP53", "BCL2"))

  result <- fgsea_to_gseaResult(fgsea_res, gene_list, gene_sets, term2name = term2name)
  res_df <- result@result

  expect_equal(res_df$Description[1], "apoptotic process")
})

test_that("fgsea_to_gseaResult preserves gene list", {
  skip_if_not_installed("fgsea")
  skip_if_not_installed("clusterProfiler")

  fgsea_res <- data.frame(
    pathway = c("PathA"),
    pval    = c(0.01),
    padj    = c(0.05),
    log2err = c(0.5),
    ES      = c(0.5),
    NES     = c(1.5),
    size    = c(20),
    stringsAsFactors = FALSE
  )
  fgsea_res$leadingEdge <- list(c("TP53"))

  gene_list <- sort(setNames(rnorm(20), paste0("G", 1:20)), decreasing = TRUE)
  gene_sets <- list(PathA = c("TP53"))

  result <- fgsea_to_gseaResult(fgsea_res, gene_list, gene_sets)
  expect_equal(result@geneList, gene_list)
})
