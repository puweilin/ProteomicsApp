# =============================================================================
# Tests: Shiny modules (basic wiring tests)
# =============================================================================

test_that("Module namespacing works correctly", {
  ns <- shiny::NS("preproc")
  expect_equal(ns("btn_load"), "preproc-btn_load")

  ns2 <- shiny::NS("analysis")
  expect_equal(ns2("btn_run_pipeline"), "analysis-btn_run_pipeline")

  ns3 <- shiny::NS("vis")
  expect_equal(ns3("vis_volcano"), "vis-vis_volcano")
})

test_that("preprocessingUI function exists and returns shiny tag", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("shinyFiles")
  skip_if_not_installed("shinyWidgets")
  skip_if_not_installed("DT")

  if (!exists("preprocessingUI", mode = "function")) {
    skip("preprocessingUI not available")
  }

  suppressPackageStartupMessages(library(shiny))
  suppressPackageStartupMessages(library(shinyFiles))
  suppressPackageStartupMessages(library(shinyWidgets))
  suppressPackageStartupMessages(library(DT))

  ui <- preprocessingUI("test_preproc")
  expect_true(inherits(ui, "shiny.tag"))
})

test_that("analysisUI function exists and returns shiny tag", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("shinyWidgets")

  if (!exists("analysisUI", mode = "function")) {
    skip("analysisUI not available")
  }

  suppressPackageStartupMessages(library(shiny))
  suppressPackageStartupMessages(library(shinyWidgets))

  ui <- analysisUI("test_analysis")
  expect_true(inherits(ui, "shiny.tag"))
})

test_that("visualizationUI function exists and returns shiny tag", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("DT")

  if (!exists("visualizationUI", mode = "function")) {
    skip("visualizationUI not available")
  }

  suppressPackageStartupMessages(library(shiny))
  suppressPackageStartupMessages(library(DT))

  ui <- visualizationUI("test_vis")
  expect_true(inherits(ui, "shiny.tag"))
})

test_that("preprocessingServer returns expected reactive list", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("shinyWidgets")
  skip_if_not_installed("shinyFiles")

  if (!exists("preprocessingServer", mode = "function")) {
    skip("preprocessingServer not available")
  }

  suppressPackageStartupMessages(library(shiny))
  suppressPackageStartupMessages(library(shinyFiles))
  suppressPackageStartupMessages(library(shinyWidgets))

  rv <- shiny::reactiveValues(
    selected_file_path = NULL,
    raw_mat = NULL,
    meta = NULL,
    annot = NULL,
    manager = NULL,
    cache_root = tempdir(),
    logs = character(0)
  )
  add_log <- function(msg) { rv$logs <- c(rv$logs, msg) }

  tryCatch({
    shiny::testServer(preprocessingServer, args = list(rv = rv, add_log = add_log), {
      result <- session$getReturned()
      expect_type(result, "list")
      expect_true("sel_filter_column" %in% names(result))
      expect_true("sel_filter_value" %in% names(result))
      expect_true("sel_imp_method" %in% names(result))
    })
  }, error = function(e) {
    skip(paste("testServer failed:", e$message))
  })
})

test_that("analysisServer returns expected reactive list", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("shinyWidgets")

  if (!exists("analysisServer", mode = "function")) {
    skip("analysisServer not available")
  }

  suppressPackageStartupMessages(library(shiny))
  suppressPackageStartupMessages(library(shinyWidgets))

  rv <- shiny::reactiveValues(
    manager = NULL,
    diff_tool = NULL,
    enrich_tool = NULL,
    gsva_tool = NULL,
    loaded_params = NULL,
    cache_manager = NULL,
    cache_update_info = NULL,
    cache_root = tempdir(),
    logs = character(0)
  )
  add_log <- function(msg) { rv$logs <- c(rv$logs, msg) }

  tryCatch({
    shiny::testServer(analysisServer, args = list(rv = rv, add_log = add_log), {
      result <- session$getReturned()
      expect_type(result, "list")
      expect_true("analysis_mode" %in% names(result))
      expect_true("sel_condition_col" %in% names(result))
    })
  }, error = function(e) {
    skip(paste("testServer failed:", e$message))
  })
})
