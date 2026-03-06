# ==============================================================================
# 0. Auto-dependency management (Auto-Installation)
# ==============================================================================
setup_environment <- function() {
  cran_packages <- c(
    "shiny", "shinydashboard", "shinyWidgets", "shinyjs", "shinybusy", "shinythemes", "shinyFiles",
    "tidyverse", "readxl", "openxlsx", "DT", "here",
    "ggplot2", "ggrepel", "patchwork", "R6", "broom", "ggsci", "circlize",
    "missForest", "doParallel", "glmnet", "caret", "splines",
    "pheatmap", "viridis", "digest", "jsonlite", "ggplotify", "vegan", "cowplot",
    "rmarkdown"
  )

  bioc_packages <- c(
    "DEP", "SummarizedExperiment", "limma",
    "clusterProfiler", "enrichplot", "org.Hs.eg.db", "ReactomePA",
    "Mfuzz", "ComplexHeatmap", "GSVA", "fgsea", "BiocParallel", "GSEABase", "KEGGREST"
  )

  options(repos = c(CRAN = "https://cloud.r-project.org"))

  all_pkgs <- c(cran_packages, bioc_packages)
  suppressPackageStartupMessages({
    lapply(all_pkgs, function(x) tryCatch(library(x, character.only = TRUE), error = function(e) message(paste("Error loading", x))))
  })
}
setup_environment()

# --- Load core scripts (package paths) ---
pkg_scripts_dir <- system.file("scripts", package = "ProteomicsApp")
if (pkg_scripts_dir == "") {
  # Development mode: load from relative path
  pkg_scripts_dir <- file.path(dirname(getwd()), "scripts")
  if (!dir.exists(pkg_scripts_dir)) {
    pkg_scripts_dir <- here::here("scripts")
  }
}
source(file.path(pkg_scripts_dir, "GeneSetCacheManager.R"))
source(file.path(pkg_scripts_dir, "ProteomicsAnalysis.R"))
source(file.path(pkg_scripts_dir, "ProteomicsGSVA.R"))
source(file.path(pkg_scripts_dir, "run_proteomics_pipeline.R"))

# --- Load Shiny modules ---
pkg_modules_dir <- system.file("app", "modules", package = "ProteomicsApp")
if (pkg_modules_dir == "") {
  pkg_modules_dir <- file.path(getwd(), "modules")
  if (!dir.exists(pkg_modules_dir)) {
    pkg_modules_dir <- file.path(dirname(getwd()), "app", "modules")
  }
}
source(file.path(pkg_modules_dir, "mod_preprocessing.R"))
source(file.path(pkg_modules_dir, "mod_analysis.R"))
source(file.path(pkg_modules_dir, "mod_visualization.R"))

options(shiny.maxRequestSize = 100 * 1024^2)

# ==============================================================================
# UI
# ==============================================================================
ui <- navbarPage(
  title = div(
    img(src = "", height = "0px"),
    "Proteomics Analysis Platform"
  ),
  id = "main_nav",
  theme = shinythemes::shinytheme("flatly"),
  header = tagList(
    useShinyjs(),
    add_busy_spinner(spin = "fading-circle", color = "#3498db", position = "top-right"),
    tags$head(
      tags$style(HTML("
        /* Global styles */
        body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; background-color: #f5f7fa; }
        .navbar { box-shadow: 0 2px 10px rgba(0,0,0,0.1); }
        .navbar-brand { font-weight: 700; letter-spacing: 1px; }

        /* Panel and card styles */
        .well { background-color: #ffffff; border: 1px solid #e8ecf1; border-radius: 8px;
                box-shadow: 0 2px 8px rgba(0,0,0,0.04); }
        .tab-content { background-color: #ffffff; padding: 15px; border-radius: 0 0 8px 8px;
                       border: 1px solid #e8ecf1; border-top: none; }
        .nav-tabs > li.active > a { border-top: 3px solid #3498db !important; font-weight: 600; }

        /* Sidebar styles */
        .sidebar-panel, .well { transition: all 0.3s ease; }
        .sidebar-panel:hover, .well:hover { box-shadow: 0 4px 12px rgba(0,0,0,0.08); }

        /* Button style optimization */
        .btn { border-radius: 6px; font-weight: 500; transition: all 0.2s ease; }
        .btn:hover { transform: translateY(-1px); box-shadow: 0 4px 8px rgba(0,0,0,0.15); }
        .btn-primary { background: linear-gradient(135deg, #3498db, #2980b9); border: none; }
        .btn-danger { background: linear-gradient(135deg, #e74c3c, #c0392b); border: none; }
        .btn-warning { background: linear-gradient(135deg, #f39c12, #e67e22); border: none; }
        .btn-success { background: linear-gradient(135deg, #2ecc71, #27ae60); border: none; }
        .btn-info { background: linear-gradient(135deg, #00bcd4, #0097a7); border: none; }
        .btn-block { margin-bottom: 8px; }

        /* Title styles */
        h4 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 8px; margin-top: 15px; }
        h5 { color: #34495e; }

        /* Input field styles */
        .form-control { border-radius: 6px; border: 1px solid #dce4ec;
                         transition: border-color 0.3s ease; }
        .form-control:focus { border-color: #3498db; box-shadow: 0 0 5px rgba(52,152,219,0.3); }
        .selectize-input { border-radius: 6px !important; }

        /* Progress bar */
        .progress { border-radius: 10px; height: 8px; }
        .progress-bar { border-radius: 10px; }

        /* Data table styles */
        .dataTables_wrapper { font-size: 13px; }
        table.dataTable thead th { background-color: #f8f9fa; border-bottom: 2px solid #3498db; }

        /* Chart container */
        .shiny-plot-output { border-radius: 8px; overflow: hidden;
                             box-shadow: 0 1px 4px rgba(0,0,0,0.06); }

        /* Notification customization */
        .shiny-notification { border-radius: 8px; border-left: 4px solid #3498db; }

        /* Scrollbar customization */
        ::-webkit-scrollbar { width: 8px; }
        ::-webkit-scrollbar-track { background: #f1f1f1; border-radius: 4px; }
        ::-webkit-scrollbar-thumb { background: #c1c1c1; border-radius: 4px; }
        ::-webkit-scrollbar-thumb:hover { background: #a1a1a1; }

        /* Card-style panels */
        .card-panel { background: white; border-radius: 8px; padding: 20px;
                      margin-bottom: 15px; box-shadow: 0 2px 8px rgba(0,0,0,0.06); }

        /* Animation */
        .fade-in { animation: fadeIn 0.5s ease-in; }
        @keyframes fadeIn { from { opacity: 0; transform: translateY(10px); }
                           to { opacity: 1; transform: translateY(0); } }

        /* Sweet alert customization */
        .swal2-popup { border-radius: 12px !important; }

        /* Sidebar section dividers */
        hr { border-top: 1px solid #eee; margin: 15px 0; }

        /* Status badges */
        .badge-success { background-color: #2ecc71; }
        .badge-warning { background-color: #f39c12; }
        .badge-danger { background-color: #e74c3c; }

        /* Download button */
        .btn-success.shiny-download-link { display: block; text-align: center; }
      "))
    )
  ),

  # --- Tab 1: Preprocessing ---
  preprocessingUI("preproc"),

  # --- Tab 2: Analysis Config ---
  analysisUI("analysis"),

  # --- Tab 3: Visualization ---
  visualizationUI("viz")
)

# ==============================================================================
# Server
# ==============================================================================
server <- function(input, output, session) {

  # Shared reactive values across all modules
  rv <- reactiveValues(
    raw_mat = NULL, meta = NULL, annot = NULL,
    manager = NULL, diff_tool = NULL, enrich_tool = NULL, gsva_tool = NULL,
    logs = character(0),
    cache_root = NULL,
    selected_file_path = NULL,
    loaded_params = NULL,
    cache_manager = NULL,
    cache_update_info = NULL,
    # Values forwarded from preprocessing module to analysis module
    filter_col_val = NULL,
    filter_val_val = NULL,
    imp_method_val = NULL
  )

  add_log <- function(msg) {
    rv$logs <- c(paste0(format(Sys.time(), "[%H:%M:%S] "), msg), rv$logs)
  }

  # --- Initialize modules ---
  preproc_vals <- preprocessingServer("preproc", rv, add_log)
  analysis_vals <- analysisServer("analysis", rv, add_log)
  visualizationServer("viz", rv, analysis_vals)

  # --- Forward preprocessing inputs to rv for analysis module ---
  observe({
    rv$filter_col_val <- tryCatch(preproc_vals$sel_filter_column(), error = function(e) NULL)
  })
  observe({
    rv$filter_val_val <- tryCatch(preproc_vals$sel_filter_value(), error = function(e) NULL)
  })
  observe({
    rv$imp_method_val <- tryCatch(preproc_vals$sel_imp_method(), error = function(e) NULL)
  })
}

shinyApp(ui, server)
