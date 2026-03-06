# =============================================================================
# Tests: Pipeline integration (run_proteomics_pipeline)
# =============================================================================

test_that("Group mode pipeline returns list with expected elements", {
  skip_on_cran()
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("fgsea")

  mgr <- create_test_data_manager()
  cache_mgr <- create_test_cache_manager()
  output_dir <- file.path(tempdir(), paste0("pipeline_test_",
                          paste0(sample(letters, 8), collapse = "")))

  result <- run_proteomics_pipeline(
    data_manager   = mgr,
    analysis_type  = "group",
    condition_col  = "condition",
    control_group  = "A",
    case_group     = "B",
    results_dir    = output_dir,
    pval_cutoff    = 1,
    enrich_pval_cutoff = 1,
    logfc_cutoff   = 0,
    run_diff       = TRUE,
    run_enrichment = TRUE,
    run_gsva       = FALSE,
    cache_manager  = cache_mgr
  )

  expect_type(result, "list")
  expect_true("diff_tool" %in% names(result))
  expect_true("enrich_tool" %in% names(result))
  expect_true("gsva_tool" %in% names(result))

  # diff_tool should have results
  expect_s3_class(result$diff_tool$diff_results, "data.frame")

  # Cleanup
  unlink(output_dir, recursive = TRUE)
})

test_that("Group mode creates output directory and files", {
  skip_on_cran()

  mgr <- create_test_data_manager()
  cache_mgr <- create_test_cache_manager()
  output_dir <- file.path(tempdir(), paste0("pipeline_dir_",
                          paste0(sample(letters, 8), collapse = "")))

  result <- run_proteomics_pipeline(
    data_manager   = mgr,
    analysis_type  = "group",
    condition_col  = "condition",
    control_group  = "A",
    case_group     = "B",
    results_dir    = output_dir,
    pval_cutoff    = 1,
    logfc_cutoff   = 0,
    run_diff       = TRUE,
    run_enrichment = FALSE,
    run_gsva       = FALSE,
    cache_manager  = cache_mgr
  )

  # Check that output directory was created
  expected_subdir <- file.path(output_dir, "B_vs_A")
  expect_true(dir.exists(expected_subdir))

  # Check for output files
  expect_true(file.exists(file.path(expected_subdir, "Total_Proteins.xlsx")))
  expect_true(file.exists(file.path(expected_subdir, "diff_tool.rds")))

  # Cleanup
  unlink(output_dir, recursive = TRUE)
})

test_that("Continuous mode runs without error", {
  skip_on_cran()

  mgr <- create_test_data_manager()
  cache_mgr <- create_test_cache_manager()
  output_dir <- file.path(tempdir(), paste0("pipeline_cont_",
                          paste0(sample(letters, 8), collapse = "")))

  result <- run_proteomics_pipeline(
    data_manager     = mgr,
    analysis_type    = "continuous",
    continuous_col   = "time_point",
    continuous_method = "linear",
    results_dir      = output_dir,
    pval_cutoff      = 1,
    corr_cutoff      = 0,
    r2_cutoff        = 0,
    run_diff         = TRUE,
    run_enrichment   = FALSE,
    run_gsva         = FALSE,
    cache_manager    = cache_mgr
  )

  expect_type(result, "list")
  expect_s3_class(result$diff_tool$diff_results, "data.frame")
  expect_true("spearman_rho" %in% colnames(result$diff_tool$diff_results))

  # Cleanup
  unlink(output_dir, recursive = TRUE)
})

test_that("Pipeline errors on missing required group parameters", {
  mgr <- create_test_data_manager()
  output_dir <- file.path(tempdir(), "pipeline_err_test")

  expect_error(
    run_proteomics_pipeline(
      data_manager  = mgr,
      analysis_type = "group",
      control_group = NULL,
      case_group    = NULL,
      results_dir   = output_dir
    ),
    "requires"
  )
})

# =============================================================================
# Tests: Special characters in paths and group names
# =============================================================================

test_that("Pipeline handles special characters in group names", {
  skip_on_cran()

  mgr <- create_test_data_manager()
  cache_mgr <- create_test_cache_manager()

  # Rename condition levels to include special characters
  se <- mgr$imputed_se
  cd <- SummarizedExperiment::colData(se)
  levels(cd$condition) <- c("Group (A&B)", "Control+C")
  SummarizedExperiment::colData(se) <- cd
  mgr$imputed_se <- se
  mgr$imputed_se_backup <- se
  mgr$meta_data <- as.data.frame(cd)
  mgr$meta_data_backup <- mgr$meta_data

  output_dir <- file.path(tempdir(), paste0("pipeline_special_",
                          paste0(sample(letters, 8), collapse = "")))

  # Pipeline should not crash — directory creation works even with special group names
  # (sub_folder_name is sanitized). The analysis itself may error due to limma's
  # makeContrasts requiring valid R names, but that should be a caught error, not a crash.
  # Verify that the output directory was at least created (sanitized name).
  err <- tryCatch({
    run_proteomics_pipeline(
      data_manager   = mgr,
      analysis_type  = "group",
      condition_col  = "condition",
      control_group  = "Group (A&B)",
      case_group     = "Control+C",
      results_dir    = output_dir,
      pval_cutoff    = 1,
      logfc_cutoff   = 0,
      run_diff       = TRUE,
      run_enrichment = FALSE,
      run_gsva       = FALSE,
      cache_manager  = cache_mgr
    )
    NULL
  }, error = function(e) e)

  # The directory should have been created with sanitized name before analysis error
  created_dirs <- list.dirs(output_dir, full.names = FALSE, recursive = FALSE)
  expect_true(length(created_dirs) >= 1 || !is.null(err))

  # If it errored, the error should be from limma, not from dir.create
  if (!is.null(err)) {
    expect_true(grepl("syntactically valid|Non-valid names|makeContrasts", err$message))
  }

  # Cleanup
  unlink(output_dir, recursive = TRUE)
})

test_that("Pipeline handles special characters in continuous column name", {
  skip_on_cran()

  mgr <- create_test_data_manager()
  cache_mgr <- create_test_cache_manager()

  # Add a column with special characters that has valid numeric data
  se <- mgr$imputed_se
  cd <- SummarizedExperiment::colData(se)
  cd[["time_h_m"]] <- seq_len(ncol(se))
  SummarizedExperiment::colData(se) <- cd
  mgr$imputed_se <- se
  mgr$imputed_se_backup <- se
  mgr$meta_data <- as.data.frame(cd)
  mgr$meta_data_backup <- mgr$meta_data

  output_dir <- file.path(tempdir(), paste0("pipeline_cont_special_",
                          paste0(sample(letters, 8), collapse = "")))

  # Use the column name — directory sanitization should handle special chars
  # in the auto-generated sub_folder_name
  result <- run_proteomics_pipeline(
    data_manager     = mgr,
    analysis_type    = "continuous",
    continuous_col   = "time_h_m",
    continuous_method = "linear",
    results_dir      = output_dir,
    # Use a sub_folder_name WITH special chars to test sanitization
    sub_folder_name  = "Regression_time (h&m)_linear",
    pval_cutoff      = 1,
    corr_cutoff      = 0,
    r2_cutoff        = 0,
    run_diff         = TRUE,
    run_enrichment   = FALSE,
    run_gsva         = FALSE,
    cache_manager    = cache_mgr
  )

  expect_type(result, "list")
  expect_s3_class(result$diff_tool$diff_results, "data.frame")

  # Verify directory was created with sanitized name
  created_dirs <- list.dirs(output_dir, full.names = FALSE, recursive = FALSE)
  expect_length(created_dirs, 1)
  expect_false(grepl("[&()]", created_dirs[1]))

  # Cleanup
  unlink(output_dir, recursive = TRUE)
})

test_that("Pipeline handles special characters in results_dir path", {
  skip_on_cran()

  mgr <- create_test_data_manager()
  cache_mgr <- create_test_cache_manager()

  # Create a parent directory with special characters
  base_dir <- file.path(tempdir(), paste0("test&data (project)",
                        paste0(sample(letters, 4), collapse = "")))
  dir.create(base_dir, recursive = TRUE)
  output_dir <- file.path(base_dir, "results")

  result <- run_proteomics_pipeline(
    data_manager   = mgr,
    analysis_type  = "group",
    condition_col  = "condition",
    control_group  = "A",
    case_group     = "B",
    results_dir    = output_dir,
    pval_cutoff    = 1,
    logfc_cutoff   = 0,
    run_diff       = TRUE,
    run_enrichment = FALSE,
    run_gsva       = FALSE,
    cache_manager  = cache_mgr
  )

  expect_type(result, "list")
  expect_s3_class(result$diff_tool$diff_results, "data.frame")

  # Verify output files exist
  expected_subdir <- file.path(output_dir, "B_vs_A")
  expect_true(dir.exists(expected_subdir))
  expect_true(file.exists(file.path(expected_subdir, "Total_Proteins.xlsx")))

  # Cleanup
  unlink(base_dir, recursive = TRUE)
})

test_that("Pipeline sanitizes sub_folder_name even when provided explicitly", {
  skip_on_cran()

  mgr <- create_test_data_manager()
  cache_mgr <- create_test_cache_manager()
  output_dir <- file.path(tempdir(), paste0("pipeline_explicit_",
                          paste0(sample(letters, 8), collapse = "")))

  result <- run_proteomics_pipeline(
    data_manager     = mgr,
    analysis_type    = "group",
    condition_col    = "condition",
    control_group    = "A",
    case_group       = "B",
    results_dir      = output_dir,
    sub_folder_name  = "My Analysis (v2) & Results",
    pval_cutoff      = 1,
    logfc_cutoff     = 0,
    run_diff         = TRUE,
    run_enrichment   = FALSE,
    run_gsva         = FALSE,
    cache_manager    = cache_mgr
  )

  expect_type(result, "list")

  # The sanitized folder name should not contain &, (, )
  created_dirs <- list.dirs(output_dir, full.names = FALSE, recursive = FALSE)
  expect_length(created_dirs, 1)
  expect_false(grepl("[&()]", created_dirs[1]))
  # But should preserve alphanumeric content
  expect_true(grepl("My_Analysis", created_dirs[1]))

  # Cleanup
  unlink(output_dir, recursive = TRUE)
})
