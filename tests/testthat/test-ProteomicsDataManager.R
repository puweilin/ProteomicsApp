# =============================================================================
# Tests: ProteomicsDataManager
# =============================================================================

test_that("ProteomicsDataManager initialization stores inputs", {
  mat   <- create_test_matrix()
  meta  <- create_test_metadata()
  annot <- create_test_annotation()

  mgr <- ProteomicsDataManager$new(mat, meta, annot, tolerate_missing_percent = 0.5)

  expect_identical(mgr$raw_mat, mat)
  expect_identical(mgr$meta_data, meta)
  expect_identical(mgr$annot_data, annot)
  expect_equal(mgr$tolerate_missing_percent, 0.5)
})

test_that("process_data() creates SE object with correct dimensions", {
  skip_if_not_installed("DEP")
  skip_if_not_installed("SummarizedExperiment")

  # VSN normalization needs >= 42 rows (minDataPointsPerStratum)
  n_prot <- 200
  mat   <- create_test_matrix(n_proteins = n_prot)
  meta  <- create_test_metadata()
  annot <- create_test_annotation(n_proteins = n_prot)

  mgr <- ProteomicsDataManager$new(mat, meta, annot, tolerate_missing_percent = 0.8)
  mgr$process_data()

  expect_s4_class(mgr$se_obj, "SummarizedExperiment")
  expect_equal(ncol(mgr$se_obj), 8)
  expect_gt(nrow(mgr$se_obj), 0)
})

test_that("process_data() with filter reduces samples", {
  skip_if_not_installed("DEP")
  skip_if_not_installed("SummarizedExperiment")

  n_prot <- 200
  mat   <- create_test_matrix(n_proteins = n_prot)
  meta  <- create_test_metadata()
  annot <- create_test_annotation(n_proteins = n_prot)

  mgr <- ProteomicsDataManager$new(mat, meta, annot, tolerate_missing_percent = 0.8)
  mgr$process_data(filter_col = "condition", filter_value = "A")

  expect_s4_class(mgr$se_obj, "SummarizedExperiment")
  expect_equal(ncol(mgr$se_obj), 4)
})

test_that("process_data() with invalid filter column errors", {
  mat   <- create_test_matrix()
  meta  <- create_test_metadata()
  annot <- create_test_annotation()

  mgr <- ProteomicsDataManager$new(mat, meta, annot, tolerate_missing_percent = 0.8)
  expect_error(mgr$process_data(filter_col = "nonexistent", filter_value = "X"),
               "not found")
})

test_that("perform_imputation() stores method and missing mask", {
  mgr <- create_test_data_manager()

  expect_equal(mgr$imputation_method, "test")
  expect_true(is.logical(mgr$missing_mask))
  expect_equal(dim(mgr$missing_mask), dim(assay(mgr$se_obj)))
})

test_that("detect_outliers() returns list with outlier_samples", {
  mgr <- create_test_data_manager()

  result <- mgr$detect_outliers(method = "pca", do_plot = FALSE)
  expect_type(result, "list")
  expect_true("outliers" %in% names(result))
  expect_type(result$outliers, "character")
})

test_that("detect_outliers() correlation method works", {
  mgr <- create_test_data_manager()

  result <- mgr$detect_outliers(method = "correlation", do_plot = FALSE)
  expect_type(result, "list")
  expect_true("outliers" %in% names(result))
  expect_true("stats" %in% names(result))
})

test_that("remove_outliers() reduces samples", {
  mgr <- create_test_data_manager()

  original_ncol <- ncol(mgr$imputed_se)
  mgr$remove_outliers("S1")

  expect_equal(ncol(mgr$imputed_se), original_ncol - 1)
  expect_false("S1" %in% colnames(mgr$imputed_se))
})

test_that("reset_data() restores backup", {
  mgr <- create_test_data_manager()

  original_ncol <- ncol(mgr$imputed_se)
  mgr$remove_outliers("S1")
  expect_equal(ncol(mgr$imputed_se), original_ncol - 1)

  mgr$reset_data()
  expect_equal(ncol(mgr$imputed_se), original_ncol)
})

test_that("subset_samples() filters by metadata column", {
  mgr <- create_test_data_manager()

  original_ncol <- ncol(mgr$imputed_se)
  mgr$subset_samples("condition", "A")

  expect_equal(ncol(mgr$imputed_se), original_ncol / 2)
})

test_that("subset_samples() with invalid column errors", {
  mgr <- create_test_data_manager()

  expect_error(mgr$subset_samples("nonexistent", "A"), "not found")
})

test_that("plot_pca() returns ggplot", {
  mgr <- create_test_data_manager()

  p <- mgr$plot_pca(color_col = "condition")
  expect_s3_class(p, "ggplot")
})

test_that("plot_protein_expression() returns ggplot for categorical", {
  mgr <- create_test_data_manager()

  proteins <- rownames(mgr$imputed_se)[1:2]
  p <- mgr$plot_protein_expression(proteins, variable = "condition")
  expect_s3_class(p, "ggplot")
})
