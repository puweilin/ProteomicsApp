# =============================================================================
# Test Fixtures: Synthetic data generators for ProteomicsApp tests
# =============================================================================
# This file is automatically loaded by testthat before each test file.

# --- Load required packages that source files use without namespace ---
suppressPackageStartupMessages({
  library(R6)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(broom)
  library(readr)
  library(splines)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(limma)
  library(DEP)
  library(stringr)
  library(openxlsx)
})

# --- Source R6 class definitions ---
pkg_dir <- system.file("scripts", package = "ProteomicsApp")
if (pkg_dir == "") {
  pkg_dir <- file.path(testthat::test_path(), "..", "..", "inst", "scripts")
}

# Source into globalenv() so R6 methods can find search-path functions
# (colData, left_join, etc.) through globalenv's parent chain
source(file.path(pkg_dir, "GeneSetCacheManager.R"), local = globalenv())
source(file.path(pkg_dir, "ProteomicsAnalysis.R"), local = globalenv())
source(file.path(pkg_dir, "ProteomicsGSVA.R"), local = globalenv())
source(file.path(pkg_dir, "run_proteomics_pipeline.R"), local = globalenv())

# Copy class generators into the current (test) environment so helper functions
# and test files can find them. testthat's test env may not traverse globalenv().
ProteomicsDataManager <- get("ProteomicsDataManager", envir = globalenv())
DiffExpAnalyst <- get("DiffExpAnalyst", envir = globalenv())
EnrichmentAnalyst <- get("EnrichmentAnalyst", envir = globalenv())
GeneSetCacheManager <- get("GeneSetCacheManager", envir = globalenv())
ProteomicsGSVA <- get("ProteomicsGSVA", envir = globalenv())
run_proteomics_pipeline <- get("run_proteomics_pipeline", envir = globalenv())
if (exists("fgsea_to_gseaResult", envir = globalenv())) {
  fgsea_to_gseaResult <- get("fgsea_to_gseaResult", envir = globalenv())
}
if (exists("MfuzzClusterer", envir = globalenv())) {
  MfuzzClusterer <- get("MfuzzClusterer", envir = globalenv())
}

# --- Source Shiny modules ---
mod_dir <- system.file("app", "modules", package = "ProteomicsApp")
if (mod_dir == "") {
  mod_dir <- file.path(testthat::test_path(), "..", "..", "inst", "app", "modules")
}
if (dir.exists(mod_dir)) {
  for (f in list.files(mod_dir, pattern = "\\.R$", full.names = TRUE)) {
    tryCatch(source(f, local = globalenv()), error = function(e) NULL)
  }
}

# --- Gene symbol pool (realistic human gene symbols) ---
GENE_SYMBOLS <- c(
  "TP53", "BRCA1", "EGFR", "MYC", "KRAS", "PIK3CA", "AKT1", "PTEN",
  "RB1", "CDH1", "VEGFA", "MAPK1", "STAT3", "JUN", "FOS", "BCL2",
  "BAX", "CASP3", "TNF", "IL6"
)

# =============================================================================
# 1. create_test_matrix: Synthetic protein expression matrix
# =============================================================================
create_test_matrix <- function(n_proteins = 20, n_samples = 8) {
  set.seed(42)
  protein_ids <- paste0("P", sprintf("%05d", seq_len(n_proteins)))
  sample_names <- paste0("S", seq_len(n_samples))

  # Generate expression values (log2 scale, centered around 20)
  mat_values <- matrix(
    rnorm(n_proteins * n_samples, mean = 20, sd = 2),
    nrow = n_proteins, ncol = n_samples
  )

  # Inject ~10% NAs
  na_idx <- sample(seq_len(n_proteins * n_samples),
                   size = round(0.10 * n_proteins * n_samples))
  mat_values[na_idx] <- NA

  df <- data.frame(ProteinID = protein_ids, mat_values,
                   check.names = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("ProteinID", sample_names)
  df
}

# =============================================================================
# 2. create_test_metadata: Sample metadata
# =============================================================================
create_test_metadata <- function(n_samples = 8) {
  data.frame(
    label      = paste0("S", seq_len(n_samples)),
    condition  = rep(c("A", "B"), each = n_samples / 2),
    replicate  = rep(seq_len(n_samples / 2), times = 2),
    batch      = rep(c("1", "2"), times = n_samples / 2),
    time_point = seq_len(n_samples),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# 3. create_test_annotation: Protein annotation mapping
# =============================================================================
create_test_annotation <- function(n_proteins = 20) {
  gene_syms <- if (n_proteins <= length(GENE_SYMBOLS)) {
    GENE_SYMBOLS[seq_len(n_proteins)]
  } else {
    c(GENE_SYMBOLS, paste0("GENE", seq_len(n_proteins - length(GENE_SYMBOLS))))
  }

  data.frame(
    ProteinID   = paste0("P", sprintf("%05d", seq_len(n_proteins))),
    Gene        = gene_syms,
    Description = paste("Description for", gene_syms),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# 4. create_test_se: Build a SummarizedExperiment directly (bypass file I/O)
# =============================================================================
create_test_se <- function(n_proteins = 20, n_samples = 8) {
  set.seed(42)
  gene_syms <- if (n_proteins <= length(GENE_SYMBOLS)) {
    GENE_SYMBOLS[seq_len(n_proteins)]
  } else {
    c(GENE_SYMBOLS, paste0("GENE", seq_len(n_proteins - length(GENE_SYMBOLS))))
  }
  gene_syms <- make.unique(gene_syms)

  sample_names <- paste0("S", seq_len(n_samples))

  # Build expression matrix — half samples in condition A, half in B
  # Add a small condition effect to first 5 proteins
  mat <- matrix(rnorm(n_proteins * n_samples, mean = 20, sd = 2),
                nrow = n_proteins, ncol = n_samples,
                dimnames = list(gene_syms, sample_names))

  # Add differential signal for first 5 proteins in group B
  half <- n_samples / 2
  mat[1:min(5, n_proteins), (half + 1):n_samples] <-
    mat[1:min(5, n_proteins), (half + 1):n_samples] + 3

  col_data <- S4Vectors::DataFrame(
    label     = sample_names,
    condition = factor(rep(c("A", "B"), each = half)),
    batch     = factor(rep(c("1", "2"), times = half)),
    time_point = seq_len(n_samples),
    row.names = sample_names
  )

  row_data <- S4Vectors::DataFrame(
    name = gene_syms,
    ID   = paste0("P", sprintf("%05d", seq_len(n_proteins))),
    row.names = gene_syms
  )

  SummarizedExperiment::SummarizedExperiment(
    assays  = list(assay = mat),
    colData = col_data,
    rowData = row_data
  )
}

# =============================================================================
# 4b. create_test_data_manager: Build a ProteomicsDataManager with imputed data
#     Bypasses process_data() which requires DEP::make_se() + normalize_vsn()
# =============================================================================
create_test_data_manager <- function(n_proteins = 20, n_samples = 8) {
  mat   <- create_test_matrix(n_proteins, n_samples)
  meta  <- create_test_metadata(n_samples)
  annot <- create_test_annotation(n_proteins)

  mgr <- ProteomicsDataManager$new(mat, meta, annot, tolerate_missing_percent = 0.8)

  # Build SE object directly (bypass process_data which needs DEP)
  se <- create_test_se(n_proteins, n_samples)
  mgr$se_obj <- se
  mgr$meta_data <- as.data.frame(SummarizedExperiment::colData(se))
  mgr$valid_protein_numbers <- nrow(se)

  # Perform imputation manually (fill NAs with row means)
  imp_mat <- SummarizedExperiment::assay(se)
  # Our test SE has no NAs, but create a missing mask anyway
  mgr$missing_mask <- matrix(FALSE, nrow = nrow(imp_mat), ncol = ncol(imp_mat))

  mgr$imputed_se <- se
  mgr$imputed_se_backup <- se
  mgr$meta_data_backup <- mgr$meta_data
  mgr$imputation_method <- "test"

  mgr
}

# =============================================================================
# 5. create_test_excel: Write a 3-sheet Excel file for upload testing
# =============================================================================
create_test_excel <- function(dir = tempdir(), n_proteins = 20, n_samples = 8) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    skip("openxlsx not available")
  }

  mat   <- create_test_matrix(n_proteins, n_samples)
  meta  <- create_test_metadata(n_samples)
  annot <- create_test_annotation(n_proteins)

  file_path <- file.path(dir, "test_proteomics.xlsx")
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Matrix")
  openxlsx::writeData(wb, "Matrix", mat)
  openxlsx::addWorksheet(wb, "Metadata")
  openxlsx::writeData(wb, "Metadata", meta)
  openxlsx::addWorksheet(wb, "Annotation")
  openxlsx::writeData(wb, "Annotation", annot)
  openxlsx::saveWorkbook(wb, file_path, overwrite = TRUE)

  file_path
}

# =============================================================================
# 6. create_test_cache_manager: Mock GeneSetCacheManager with tiny gene sets
# =============================================================================
create_test_cache_manager <- function() {
  cache_dir <- file.path(tempdir(), paste0("test_cache_",
                         paste0(sample(letters, 10), collapse = "")))
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  # Build mock gene set caches — each term needs >= 10 genes for enricher(minGSSize=10)
  all_genes <- c(GENE_SYMBOLS,
                 paste0("GENE", seq_len(10)))  # Extra genes to reach 15 per term

  mock_go_bp <- list(
    TERM2GENE = data.frame(
      term = rep(c("GO:0006915", "GO:0007049"), each = 15),
      gene = c(all_genes[1:15], all_genes[c(6:20, 21:30)][1:15]),
      stringsAsFactors = FALSE
    ),
    TERM2NAME = data.frame(
      term = c("GO:0006915", "GO:0007049"),
      name = c("apoptotic process", "cell cycle"),
      stringsAsFactors = FALSE
    )
  )

  mock_kegg <- list(
    TERM2GENE = data.frame(
      term = rep(c("hsa04110", "hsa04151"), each = 15),
      gene = c(all_genes[1:15], all_genes[c(6:20, 21:30)][1:15]),
      stringsAsFactors = FALSE
    ),
    TERM2NAME = data.frame(
      term = c("hsa04110", "hsa04151"),
      name = c("Cell cycle", "PI3K-Akt signaling pathway"),
      stringsAsFactors = FALSE
    )
  )

  mock_generic <- list(
    TERM2GENE = data.frame(
      term = rep(c("TERM1", "TERM2"), each = 15),
      gene = c(all_genes[1:15], all_genes[c(6:20, 21:30)][1:15]),
      stringsAsFactors = FALSE
    ),
    TERM2NAME = data.frame(
      term = c("TERM1", "TERM2"),
      name = c("Test pathway 1", "Test pathway 2"),
      stringsAsFactors = FALSE
    )
  )

  saveRDS(mock_go_bp, file.path(cache_dir, "GO_BP.rds"))
  saveRDS(mock_generic, file.path(cache_dir, "GO_MF.rds"))
  saveRDS(mock_generic, file.path(cache_dir, "GO_CC.rds"))
  saveRDS(mock_kegg, file.path(cache_dir, "KEGG.rds"))
  saveRDS(mock_generic, file.path(cache_dir, "Reactome.rds"))
  saveRDS(mock_generic, file.path(cache_dir, "Wiki.rds"))

  # Save metadata so is_cache_exists() returns TRUE
  cache_meta <- list(
    date_created    = Sys.time(),
    org_hs_version  = "3.18.0",
    go_db_version   = "3.18.0",
    reactome_db_version = "1.86.0",
    msigdbr_version = "7.5.1",
    kegg_date       = Sys.Date(),
    db_counts       = list(GO_BP = 10, KEGG = 10)
  )
  saveRDS(cache_meta, file.path(cache_dir, "cache_meta.rds"))

  GeneSetCacheManager$new(cache_dir = cache_dir)
}
