# ==============================================================================
# 0. 自动依赖管理 (Auto-Installation)
# ==============================================================================
setup_environment <- function() {
  cran_packages <- c(
    "shiny", "shinydashboard", "shinyWidgets", "shinyjs", "shinybusy", "shinythemes", "shinyFiles",
    "tidyverse", "readxl", "openxlsx", "DT", "here",
    "ggplot2", "ggrepel", "patchwork", "R6", "broom", "ggsci", "circlize",
    "missForest", "doParallel", "glmnet", "caret", "splines",
    "pheatmap", "viridis", "digest", "jsonlite", "shinythemes", "ggplotify"
  )
  
  bioc_packages <- c(
    "DEP", "SummarizedExperiment", "limma",                            
    "clusterProfiler", "enrichplot", "org.Hs.eg.db", "ReactomePA",     
    "Mfuzz", "ComplexHeatmap", "GSVA", "GSEABase"                                                          
  )
  
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  
  all_pkgs <- c(cran_packages, bioc_packages)
  suppressPackageStartupMessages({
    lapply(all_pkgs, function(x) tryCatch(library(x, character.only = TRUE), error = function(e) message(paste("Error loading", x))))
  })
}
setup_environment()

# --- 加载核心脚本 (包内路径) ---
pkg_scripts_dir <- system.file("scripts", package = "ProteomicsApp")
if (pkg_scripts_dir == "") {
  # 开发模式: 直接从相对路径加载
  pkg_scripts_dir <- file.path(dirname(getwd()), "scripts")
  if (!dir.exists(pkg_scripts_dir)) {
    pkg_scripts_dir <- here::here("scripts")
  }
}
source(file.path(pkg_scripts_dir, "ProteomicsAnalysis.R"))
source(file.path(pkg_scripts_dir, "ProteomicsGSVA.R"))
source(file.path(pkg_scripts_dir, "run_proteomics_pipeline.R"))

options(shiny.maxRequestSize = 200 * 1024^2)

# ==============================================================================
# UI 设计 (Navbar 布局)
# ==============================================================================
ui <- navbarPage(
  title = "Proteomics Analysis Platform",
  id = "main_nav",
  theme = shinythemes::shinytheme("flatly"),  # [新增] 使用专业主题
  header = tagList(
    useShinyjs(),
    add_busy_spinner(spin = "fading-circle", position = "bottom-right"),
    # [新增] 专业学术风格CSS
    tags$head(
      tags$style(HTML("
        /* 全局样式 */
        body {
          font-family: 'Segoe UI', 'Helvetica Neue', Arial, sans-serif;
          font-size: 14px;
          color: #2c3e50;
          background-color: #f8f9fa;
        }

        /* Navbar 样式 */
        .navbar-default {
          background-color: #34495e;
          border-color: #2c3e50;
          box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }

        .navbar-default .navbar-brand,
        .navbar-default .navbar-nav > li > a {
          color: #ecf0f1;
          font-weight: 500;
        }

        .navbar-default .navbar-nav > .active > a {
          background-color: #2c3e50;
          color: #3498db;
        }

        /* 面板和卡片样式 */
        .well {
          background-color: white;
          border: 1px solid #e1e8ed;
          border-radius: 6px;
          box-shadow: 0 1px 3px rgba(0,0,0,0.06);
        }

        /* 侧边栏样式 */
        .col-sm-3 {
          background-color: white;
          border-radius: 6px;
          padding: 20px;
          box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        }

        /* 按钮样式优化 */
        .btn {
          border-radius: 4px;
          font-weight: 500;
          padding: 8px 16px;
          transition: all 0.3s ease;
          border: none;
        }

        .btn-primary {
          background-color: #3498db;
          box-shadow: 0 2px 4px rgba(52, 152, 219, 0.3);
        }

        .btn-primary:hover {
          background-color: #2980b9;
          transform: translateY(-2px);
          box-shadow: 0 4px 8px rgba(52, 152, 219, 0.4);
        }

        .btn-warning {
          background-color: #f39c12;
          color: white;
          box-shadow: 0 2px 4px rgba(243, 156, 18, 0.3);
        }

        .btn-warning:hover {
          background-color: #e67e22;
          transform: translateY(-2px);
        }

        .btn-danger {
          background-color: #e74c3c;
          box-shadow: 0 2px 4px rgba(231, 76, 60, 0.3);
        }

        .btn-danger:hover {
          background-color: #c0392b;
          transform: translateY(-2px);
        }

        .btn-info {
          background-color: #16a085;
          box-shadow: 0 2px 4px rgba(22, 160, 133, 0.3);
        }

        /* 标题样式 */
        h4 {
          color: #2c3e50;
          font-weight: 600;
          margin-bottom: 15px;
          padding-bottom: 10px;
          border-bottom: 2px solid #ecf0f1;
        }

        /* Tab 样式 - 主标签 (Preprocessing, Analysis Config, Visualization) */
        .nav-tabs > li > a {
          color: #7f8c8d;
          font-weight: 600;
          font-size: 16px;
          padding: 12px 20px;
        }

        .nav-tabs > li.active > a {
          color: #2c3e50;
          font-weight: 700;
          font-size: 16px;
          border-top: 3px solid #3498db;
          background-color: #f8f9fa;
        }

        /* 子标签样式 (Data Overview, Imputation QC, etc.) */
        .tab-content .nav-tabs > li > a {
          font-size: 14px;
          font-weight: 500;
          padding: 10px 15px;
        }

        .tab-content .nav-tabs > li.active > a {
          font-size: 14px;
          font-weight: 600;
          color: #3498db;
          border-top: 2px solid #3498db;
        }

        /* 输入框样式 */
        .form-control {
          border: 1px solid #dfe6e9;
          border-radius: 4px;
          transition: border-color 0.3s ease;
        }

        .form-control:focus {
          border-color: #3498db;
          box-shadow: 0 0 0 0.2rem rgba(52, 152, 219, 0.15);
        }

        /* Help text */
        .help-block {
          color: #95a5a6;
          font-size: 12px;
          font-style: italic;
          margin-top: 5px;
        }

        /* 进度条 */
        .progress {
          height: 25px;
          border-radius: 4px;
          box-shadow: inset 0 1px 2px rgba(0,0,0,0.1);
        }

        .progress-bar {
          background-color: #3498db;
          font-weight: 500;
        }

        /* 数据表格样式 */
        .dataTables_wrapper {
          font-size: 13px;
        }

        .dataTables_wrapper .dataTables_paginate .paginate_button.current {
          background: #3498db;
          border-color: #3498db;
        }

        /* Alert 框样式 */
        .alert {
          border-radius: 4px;
          border-left: 4px solid;
        }

        .alert-success {
          border-left-color: #27ae60;
        }

        .alert-warning {
          border-left-color: #f39c12;
        }

        .alert-danger {
          border-left-color: #e74c3c;
        }

        /* 图表容器 */
        .shiny-plot-output {
          border: 1px solid #ecf0f1;
          border-radius: 4px;
          background: white;
          padding: 10px;
          box-shadow: 0 1px 3px rgba(0,0,0,0.06);
        }

        /* Value Box (如果使用shinydashboard组件) */
        .small-box {
          border-radius: 4px;
          box-shadow: 0 2px 4px rgba(0,0,0,0.08);
        }

        /* Selectize样式 */
        .selectize-input {
          border: 1px solid #dfe6e9;
          border-radius: 4px;
        }

        .selectize-input.focus {
          border-color: #3498db;
        }

        /* File input */
        .btn-file {
          background-color: #ecf0f1;
          color: #2c3e50;
          border: 1px solid #bdc3c7;
        }

        /* 响应式表格 */
        .table {
          font-size: 13px;
        }

        .table thead th {
          background-color: #ecf0f1;
          color: #2c3e50;
          font-weight: 600;
          border-bottom: 2px solid #bdc3c7;
        }

        .table-striped tbody tr:nth-of-type(odd) {
          background-color: #f8f9fa;
        }

        /* pickerInput 样式修复 - 确保文字可见 */
        .bootstrap-select .dropdown-toggle {
          background-color: #fff !important;
          color: #2c3e50 !important;
          border: 1px solid #dfe6e9;
        }

        .bootstrap-select .dropdown-toggle .filter-option-inner-inner {
          color: #2c3e50 !important;
        }

        .bootstrap-select .dropdown-menu li a {
          color: #2c3e50 !important;
        }

        .bootstrap-select .dropdown-menu li.selected a {
          background-color: #3498db !important;
          color: #fff !important;
        }

        .bootstrap-select .dropdown-menu li a:hover {
          background-color: #ecf0f1 !important;
          color: #2c3e50 !important;
        }

        .bootstrap-select .dropdown-menu li.active a,
        .bootstrap-select .dropdown-menu li.active a:hover {
          background-color: #d6eaf8 !important;
          color: #2c3e50 !important;
        }

        .bootstrap-select .dropdown-menu .text mark,
        .bootstrap-select .dropdown-menu .text .text-muted {
          background-color: #f9e79f !important;
          color: #2c3e50 !important;
          padding: 0 2px;
        }

        .bootstrap-select .bs-searchbox input {
          color: #2c3e50 !important;
          background-color: #fff !important;
        }
      "))
    )
  ),
  
  # --- Tab 1: Preprocessing ---
  tabPanel("Preprocessing", icon = icon("vial"),
           sidebarLayout(
             sidebarPanel(width = 3,
                          h4(icon("upload"), "Data Import"),
                          shinyFilesButton("upload_file", "Select Excel File",
                                          title = "Please select an Excel file (.xlsx)",
                                          multiple = FALSE,
                                          buttonType = "primary",
                                          class = "btn-block"),
                          verbatimTextOutput("selected_file_path", placeholder = TRUE),
                          helpText("Sheet1: Matrix; Sheet2: Metadata; Sheet3: Annot"),

                          # [修改] 添加metadata列选择和值选择
                          uiOutput("ui_filter_column_selector"),
                          uiOutput("ui_filter_value_selector"),

                          hr(),
                          h4(icon("filter"), "Quality Control"),
                          sliderInput("missing_rate_threshold", "Max Missing Rate per Protein:",
                                     min = 0, max = 100, value = 70, step = 5,
                                     post = "%",
                                     ticks = TRUE),
                          helpText("Proteins with missing values above this threshold will be removed."),

                          uiOutput("ui_manager_cache_alert"),
                          checkboxInput("use_manager_cache", "Load Processed Data (if avail.)", TRUE),

                          actionButton("btn_load_process", "Load Data", icon = icon("play"), class = "btn-primary", width = "100%"),

                          hr(),
                          h4(icon("magic"), "Imputation"),
                          selectInput("sel_imp_method", "Method:", 
                                      choices = c("missForest (Accurate)" = "missForest", "MinProb" = "MinProb", "QRILC" = "QRILC", 
                                                  "MinDet" = "MinDet", "KNN" = "knn", "BPCA" = "bpca", "Zero" = "zero", "Minimum" = "min"),          
                                      selected = "missForest"),
                          actionButton("btn_impute", "Run Imputation", icon = icon("bolt"), class = "btn-warning", width = "100%"),
                          
                          hr(),
                          h4(icon("search-minus"), "Outlier Removal"),
                          helpText("Step 1: Auto-detect or manually select samples."),
                          actionButton("btn_detect_outlier", "Auto-Detect (PCA)", width = "100%"),
                          br(), br(),
                          pickerInput("sel_outliers", "Select Samples to Remove:",
                                      choices = NULL,
                                      selected = NULL,
                                      multiple = TRUE,
                                      options = list(`actions-box` = TRUE, `live-search` = TRUE,
                                                     `selected-text-format` = "count > 3",
                                                     `none-selected-text` = "No samples loaded yet")),
                          br(),
                          helpText("Step 2: Click Confirm to apply changes."),
                          actionButton("btn_confirm_remove", "Confirm Update", class = "btn-danger", width = "100%"),
                          hr(),
                          htmlOutput("txt_outlier_status")
             ),
             mainPanel(width = 9,
                       tabsetPanel(id = "subtab_wizard",
                                   tabPanel("Data Overview", icon = icon("table"), br(), verbatimTextOutput("data_summary"), hr(), h4("Initial Missing Values Pattern"), plotOutput("missing_pattern_plot", height = "500px")),
                                   tabPanel("Imputation QC", icon = icon("check-double"), br(), h4("Imputation Quality Assessment"), plotOutput("qc_imp_assess", height = "600px")),
                                   tabPanel("PCA (Outlier Check)", icon = icon("crosshairs"), br(),
                                           h4("PCA Overview"),
                                           fluidRow(
                                             column(6, plotOutput("qc_pca_plot", height = "450px")),
                                             column(6, plotOutput("qc_pca_outlier_plot", height = "450px"))
                                           ),
                                           hr(),
                                           h4("Sample Correlation (Outlier Detection)"),
                                           plotOutput("qc_cor_outlier_plot", height = "400px")
                                  ),
                                   tabPanel("Clean Matrix", icon = icon("th"), br(),
                                           fluidRow(
                                             column(12,
                                                   h4("Filtered & Imputed Protein Expression Matrix"),
                                                   helpText("Red values indicate imputed data. First column: Protein ID, Second column: Gene Name, Remaining columns: Sample expression values."),
                                                   downloadButton("download_clean_matrix", "Download Matrix (Excel)", class = "btn-success"),
                                                   hr(),
                                                   DTOutput("clean_matrix_table")
                                             )
                                           ))
                       )
             )
           )
  ),
  
  # --- Tab 2: Analysis Config ---
  tabPanel("Analysis Config", icon = icon("cogs"),
           sidebarLayout(
             sidebarPanel(width = 3,
                          h4("Analysis Settings"),
                          radioButtons("analysis_mode", "Mode:", 
                                       choices = c("Group Comparison" = "group", 
                                                   "ANOVA (Multi-group)" = "anova",
                                                   "Continuous Regression" = "continuous")),
                          uiOutput("ui_analysis_params"),
                          
                          hr(),
                          h4("Thresholds & Params"),
                          numericInput("param_pval", "P-val Cutoff:", 0.05, step = 0.01),
                          conditionalPanel("input.analysis_mode == 'group'",
                            numericInput("param_logfc", "LogFC Cutoff:", 0.263, step = 0.1)
                          ),
                          conditionalPanel("input.analysis_mode == 'continuous'",
                            numericInput("param_corr", "Correlation Cutoff (Rho):", 0.5, step = 0.1, min = 0, max = 1),
                            numericInput("param_r2", "R-squared Cutoff:", 0.5, step = 0.1, min = 0, max = 1)
                          ),
                          numericInput("param_enrich_pval", "Enrich P-val:", 0.05, step = 0.01),
                          checkboxInput("param_use_adj", "Use Adjusted P-val?", TRUE),
                          checkboxInput("run_gsva", "Run GSVA? (Time Consuming)", TRUE),
                          
                          hr(),
                          h4("Result Management"),
                          uiOutput("ui_analysis_history"),
                          splitLayout(
                            actionButton("btn_run_pipeline", "Run New Analysis", icon = icon("rocket"), class = "btn-danger", width = "100%"),
                            uiOutput("ui_btn_load_analysis")
                          )
             ),
             mainPanel(width = 9,
                       h4("Analysis Log Console"),
                       verbatimTextOutput("analysis_log")
             )
           )
  ),
  
  # --- Tab 3: Visualization ---
  tabPanel("Visualization", icon = icon("chart-bar"),
           sidebarLayout(
             sidebarPanel(width = 3,
                          h4("Visual Controls"),
                          htmlOutput("txt_current_params"),
                          hr(),
                          uiOutput("vis_slider_fc"),
                          fluidRow(
                            column(8, sliderInput("vis_slider_pval", "P-val Cutoff (View):", min=0, max=0.2, value=0.05, step=0.001)),
                            column(4, numericInput("vis_input_pval", "Value:", value=0.05, min=0, max=1, step=0.001))
                          ),
                          hr()
             ),
             mainPanel(width = 9,
                       tabsetPanel(
                         tabPanel("Volcano & Table", br(), fluidRow(column(7, plotOutput("vis_volcano", brush = "volcano_brush", height = "500px")), column(5, h5("Selected Points:"), tableOutput("vis_volcano_info"))), hr(), DTOutput("vis_dep_table")),
                         tabPanel("Heatmap", br(), fluidRow(column(4, numericInput("heatmap_n", "Top N Proteins:", 50, min = 10, max = 2000)), column(4, checkboxInput("heat_cluster_row", "Cluster Rows", TRUE)), column(4, checkboxInput("heat_cluster_col", "Cluster Cols", TRUE))), plotOutput("vis_heatmap", height = "700px")),
                         tabPanel("PCA", plotOutput("vis_pca", height = "600px")),
                         tabPanel("Deep Dive", br(),
                                  fluidRow(
                                    column(6, uiOutput("vis_selector_protein")),
                                    column(6, helpText("Select up to 9 proteins to display simultaneously"))
                                  ),
                                  plotOutput("vis_single_prot", height = "700px")
                         ),
                         tabPanel("Enrichment (ORA)", br(), fluidRow(column(6, selectInput("vis_ora_db", "Database:", choices = c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome", "Wiki"))), column(6, radioButtons("vis_ora_dir", "Direction:", choices = c("UP", "DOWN"), inline = TRUE))), plotOutput("vis_ora_plot", height = "600px"), DTOutput("vis_ora_table")),
                         tabPanel("Enrichment (GSEA)", br(), fluidRow(column(4, selectInput("vis_gsea_db", "Database:", choices = c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome"))), column(8, uiOutput("ui_vis_gsea_pathway"))), plotOutput("vis_gsea_dotplot", height = "500px"), hr(), plotOutput("vis_gsea_curve", height = "400px")),
                         tabPanel("GSVA Analysis", br(), fluidRow(column(4, selectInput("vis_gsva_db", "Database:", choices = c("GOBP", "GOMF", "GOCC", "KEGG", "Wiki", "Reactome"))), column(4, selectInput("vis_gsva_plot_type", "Plot Type:", choices = c("Volcano", "Heatmap", "Barplot")))), plotOutput("vis_gsva_plot", height = "700px"))
                       )
             )
           )
  )
)

# ==============================================================================
# Server 逻辑
# ==============================================================================
server <- function(input, output, session) {
  
  rv <- reactiveValues(
    raw_mat = NULL, meta = NULL, annot = NULL,
    manager = NULL, diff_tool = NULL, enrich_tool = NULL, gsva_tool = NULL,
    logs = character(0),
    cache_root = NULL,
    selected_file_path = NULL,  # 存储用户选择的文件真实路径
    loaded_params = NULL 
  )
  
  add_log <- function(msg) {
    rv$logs <- c(paste0(format(Sys.time(), "[%H:%M:%S] "), msg), rv$logs)
  }
  output$analysis_log <- renderText({ paste(rv$logs, collapse = "\n") })

  # --- shinyFiles: 文件选择器初始化 ---
  # 使用普通环境变量存储上次目录 (shinyFileChoose 的 roots 回调不在 reactive 上下文中)
  file_mem <- environment()
  file_mem$last_dir <- NULL

  get_volumes <- function() {
    base_vols <- c(Home = path.expand("~"), getVolumes()())
    ld <- file_mem$last_dir
    if (!is.null(ld) && dir.exists(ld) && ld != path.expand("~")) {
      c("Last Used" = ld, base_vols)
    } else {
      base_vols
    }
  }

  shinyFileChoose(input, "upload_file", roots = get_volumes, filetypes = c("xlsx"))

  # 显示选中的文件路径
  output$selected_file_path <- renderText({
    if (is.null(rv$selected_file_path)) {
      "No file selected"
    } else {
      basename(rv$selected_file_path)
    }
  })

  # [新增] 生成Clean Matrix的辅助函数
  generate_clean_matrix <- function(manager) {
    req(manager, manager$imputed_se)

    tryCatch({
      # 获取imputed数据矩阵
      imp_mat <- assay(manager$imputed_se)
      missing_mask <- manager$missing_mask

      # 确保missing_mask可用且维度匹配
      if (!is.null(missing_mask) && nrow(missing_mask) == nrow(imp_mat)) {
        # 如果删除了样本，需要调整missing_mask
        current_samples <- colnames(imp_mat)
        backup_samples <- if (!is.null(manager$imputed_se_backup)) colnames(manager$imputed_se_backup) else current_samples

        if (length(current_samples) < length(backup_samples)) {
          sample_indices <- match(current_samples, backup_samples)
          sample_indices <- sample_indices[!is.na(sample_indices)]
          if (length(sample_indices) == ncol(imp_mat)) {
            missing_mask <- missing_mask[, sample_indices, drop = FALSE]
          } else {
            missing_mask <- matrix(FALSE, nrow = nrow(imp_mat), ncol = ncol(imp_mat))
          }
        }
      } else {
        # 如果missing_mask不可用，创建一个全FALSE的mask
        missing_mask <- matrix(FALSE, nrow = nrow(imp_mat), ncol = ncol(imp_mat))
        add_log("Warning: Missing mask not available, imputed values will not be highlighted.")
      }

      # [修复] 获取正确的 Protein ID 和 Gene Name
      # rowData(manager$imputed_se)$ID 是 Protein ID
      # rowData(manager$imputed_se)$name 是 Gene Name
      protein_ids <- rowData(manager$imputed_se)$ID
      gene_names <- rowData(manager$imputed_se)$name
      if (is.null(protein_ids)) protein_ids <- rownames(imp_mat)
      if (is.null(gene_names)) gene_names <- rownames(imp_mat)

      # [修复] 将列名从 DEP 编号转为 metadata 中的 label
      col_meta <- as.data.frame(colData(manager$imputed_se))
      current_colnames <- colnames(imp_mat)

      # 获取 label 映射
      if ("label" %in% colnames(col_meta)) {
        # colData 的 rownames 是当前样本名，label 列是原始样本名
        label_mapping <- col_meta$label
        names(label_mapping) <- rownames(col_meta)

        # 替换列名
        new_colnames <- sapply(current_colnames, function(x) {
          if (x %in% names(label_mapping)) label_mapping[x] else x
        })
        colnames(imp_mat) <- new_colnames
        colnames(missing_mask) <- new_colnames
      }

      # 构建数据框
      clean_df <- data.frame(
        ProteinID = protein_ids,
        GeneName = gene_names,
        imp_mat,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )

      # 计算 imputed 值统计信息
      imputed_count <- sum(missing_mask)
      imputed_pct <- round(imputed_count / (nrow(missing_mask) * ncol(missing_mask)) * 100, 2)

      # 渲染DT表格（移除逐单元格高亮以提升性能，imputed值高亮仅在下载的Excel中显示）
      output$clean_matrix_table <- renderDT({
        datatable(
          clean_df,
          options = list(
            pageLength = 25,
            scrollX = TRUE,
            scrollY = "600px",
            dom = 'Bfrtip',
            buttons = c('copy', 'csv', 'excel'),
            fixedColumns = list(leftColumns = 2),
            columnDefs = list(
              list(className = 'dt-center', targets = 0:1)
            )
          ),
          rownames = FALSE,
          class = 'cell-border stripe hover',
          extensions = c('Buttons', 'FixedColumns'),
          caption = htmltools::tags$caption(
            style = 'caption-side: top; text-align: left; color: #666; font-size: 12px;',
            paste0('Imputed values: ', imputed_count, ' (', imputed_pct, '%). ',
                   'Download Excel to see imputed cells highlighted in red.')
          )
        )
      })

      # 下载按钮
      output$download_clean_matrix <- downloadHandler(
        filename = function() {
          paste0("clean_matrix_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
        },
        content = function(file) {
          # 创建workbook
          wb <- createWorkbook()
          addWorksheet(wb, "Clean Matrix")

          # 写入数据
          writeData(wb, "Clean Matrix", clean_df, startRow = 1, startCol = 1)

          # [优化] 使用向量化操作标记imputed值（红色填充）
          idx <- which(missing_mask, arr.ind = TRUE)
          if (nrow(idx) > 0) {
            style <- createStyle(fgFill = "#ffcccc", fontColour = "#cc0000", textDecoration = "bold")
            addStyle(wb, "Clean Matrix", style,
                     rows = idx[, 1] + 1,    # +1 因为标题行
                     cols = idx[, 2] + 2,    # +2 因为前两列是 ProteinID 和 GeneName
                     gridExpand = FALSE)
          }

          # 自动列宽
          setColWidths(wb, "Clean Matrix", cols = 1:ncol(clean_df), widths = "auto")

          # 冻结前两列
          freezePane(wb, "Clean Matrix", firstActiveRow = 2, firstActiveCol = 3)

          # 保存
          saveWorkbook(wb, file, overwrite = TRUE)
        }
      )

      add_log("Clean matrix generated successfully.")
      updateTabsetPanel(session, "subtab_wizard", selected = "Clean Matrix")

    }, error = function(e) {
      add_log(paste("Failed to generate clean matrix:", e$message))
      sendSweetAlert(session, "Error", paste("Failed to generate clean matrix:", e$message), "error")
    })
  }
  
  # --- 1. 文件选择与 Metadata 列选择 (shinyFiles) ---
  observeEvent(input$upload_file, {
    # shinyFiles 返回的是一个列表结构，需要用 parseFilePaths 解析
    if (is.integer(input$upload_file)) return()  # 初始状态为整数，忽略

    file_selected <- parseFilePaths(get_volumes(), input$upload_file)
    req(nrow(file_selected) > 0)

    path <- as.character(file_selected$datapath)
    rv$selected_file_path <- path  # 保存真实路径
    file_mem$last_dir <- dirname(path)  # 记忆上次目录

    tryCatch({
      add_log(paste("File selected:", path))

      temp_meta <- read_excel(path, sheet = 2)

      # [修改] 动态生成metadata列选择器
      meta_cols <- colnames(temp_meta)
      output$ui_filter_column_selector <- renderUI({
        selectInput("sel_filter_column", "Filter by Column:",
                    choices = c("None (All samples)" = "none", meta_cols),
                    selected = if("tissue" %in% meta_cols) "tissue" else "none")
      })

      rv$raw_mat <- read_excel(path, sheet = 1)
      rv$meta    <- temp_meta
      rv$annot   <- read_excel(path, sheet = 3)

    }, error = function(e) { sendSweetAlert(session, "File Error", e$message, "error") })
  })

  # [新增] 根据选择的列动态生成值选择器
  observeEvent(input$sel_filter_column, {
    req(rv$meta, input$sel_filter_column)

    if (input$sel_filter_column == "none") {
      output$ui_filter_value_selector <- renderUI({
        helpText("All samples will be used.", style = "color: #27ae60; font-weight: bold;")
      })
    } else {
      values <- unique(rv$meta[[input$sel_filter_column]])
      output$ui_filter_value_selector <- renderUI({
        pickerInput("sel_filter_value", paste0("Select ", input$sel_filter_column, ":"),
                    choices = values,
                    selected = values[1],
                    options = list(`live-search` = TRUE))
      })
    }
  })
  
  # --- 2. 动态缓存路径监听 ---
  observe({
    req(rv$selected_file_path, input$sel_filter_column)

    # 使用选中文件的上一层目录作为主目录
    parent_dir <- dirname(rv$selected_file_path)
    file_id <- tools::file_path_sans_ext(basename(rv$selected_file_path))

    # [修改] 根据选择的列和值生成缓存路径
    if (input$sel_filter_column == "none") {
      safe_filter <- "all_samples"
    } else {
      req(input$sel_filter_value)
      safe_filter <- paste0(gsub("[^a-zA-Z0-9]", "_", input$sel_filter_column),
                           "_",
                           gsub("[^a-zA-Z0-9]", "_", input$sel_filter_value))
    }

    # 在上传文件的同级目录下创建 results_web_session 文件夹
    rv$cache_root <- file.path(parent_dir, "results_web_session", file_id, safe_filter)

    if(!dir.exists(rv$cache_root)) dir.create(rv$cache_root, recursive = TRUE)

    manager_path <- file.path(rv$cache_root, "manager.rds")
    if (file.exists(manager_path)) {
      filter_display <- if (input$sel_filter_column == "none") {
        "All samples"
      } else {
        paste0(input$sel_filter_column, " = ", input$sel_filter_value)
      }
      output$ui_manager_cache_alert <- renderUI({
        div(style="color: #00a65a; font-weight: bold; background-color: #eafcf0; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
            icon("check-circle"), paste("Cache found for", filter_display))
      })
    } else {
      output$ui_manager_cache_alert <- renderUI({
        div(style="color: #8a6d3b; background-color: #fcf8e3; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
            icon("info-circle"), " No cache. Click Load to process.")
      })
    }
  })
  
  # --- 3. 加载数据 ---
  observeEvent(input$btn_load_process, {
    req(rv$raw_mat, input$sel_filter_column, rv$cache_root)
    manager_path <- file.path(rv$cache_root, "manager.rds")

    if (input$use_manager_cache && file.exists(manager_path)) {
      filter_display <- if (input$sel_filter_column == "none") "All samples" else paste0(input$sel_filter_column, " = ", input$sel_filter_value)
      add_log(paste("Loading cached data for", filter_display, "..."))
      withProgress(message = "Loading Pre-processed Data...", {
        rv$manager <- read_rds(manager_path)
      })
      sendSweetAlert(session, "Loaded", paste("Loaded data for", filter_display), "success")

      output$data_summary <- renderText({ paste0("Samples: ", ncol(rv$manager$se_obj), "\nProteins: ", nrow(rv$manager$se_obj)) })
      output$missing_pattern_plot <- renderPlot({ rv$manager$plot_missing_pattern() })

      # [修改] 恢复 Imputation Method 选择器显示实际使用的方法
      if (!is.null(rv$manager$imputed_se)) {
        # 从manager对象中读取实际使用的imputation方法
        actual_method <- if (!is.null(rv$manager$imputation_method)) {
          rv$manager$imputation_method
        } else {
          "missForest"  # 默认值（旧数据可能没有记录）
        }
        updateSelectInput(session, "sel_imp_method", selected = actual_method)

        # [修改] 恢复 Outlier Detection 状态并正确显示
        if (!is.null(rv$manager$imputed_se_backup)) {
          # 计算被移除的样本数
          removed_samples <- setdiff(colnames(rv$manager$imputed_se_backup),
                                     colnames(rv$manager$imputed_se))
          removed_count <- length(removed_samples)
          current_count <- ncol(rv$manager$imputed_se)

          if (removed_count > 0) {
            output$txt_outlier_status <- renderUI({
              HTML(paste0("<b style='color:green'>Active Samples: ", current_count, "</b> | ",
                         "<b style='color:red'>Removed: ", removed_count, "</b><br>",
                         "<small style='color:#7f8c8d;'>", paste(removed_samples, collapse = ", "), "</small>"))
            })
          } else {
            output$txt_outlier_status <- renderUI({
              HTML(paste0("<b style='color:green'>Active Samples: ", current_count, "</b> | ",
                         "<b>No samples removed</b>"))
            })
          }
        } else {
          # 没有backup，说明是第一次imputation或没有删除样本
          total_samples <- ncol(rv$manager$imputed_se)
          output$txt_outlier_status <- renderUI({
            HTML(paste0("<b style='color:green'>Active Samples: ", total_samples, "</b> | ",
                       "<b>No samples removed yet</b>"))
          })
        }

        output$qc_imp_assess <- renderPlot({ rv$manager$assess_imputation(plot_type = c("structure", "density"))$plots })
        output$qc_pca_plot <- renderPlot({ rv$manager$plot_pca(color_col = "condition") })
        output$qc_pca_outlier_plot <- renderPlot({
          res <- rv$manager$detect_outliers(method = "pca", sd_threshold = 3, remove = FALSE, do_plot = FALSE)
          res$plot
        })
        output$qc_cor_outlier_plot <- renderPlot({
          res <- rv$manager$detect_outliers(method = "correlation", sd_threshold = 3, remove = FALSE, do_plot = FALSE)
          res$plot
        })

        # [新增] 生成Clean Matrix
        generate_clean_matrix(rv$manager)

        # [新增] 更新 outlier picker
        update_outlier_picker()

        updateTabsetPanel(session, "subtab_wizard", selected = "Clean Matrix")
      }
      return()
    }

    tryCatch({
      filter_display <- if (input$sel_filter_column == "none") "All samples" else paste0(input$sel_filter_column, " = ", input$sel_filter_value)
      add_log(paste("Processing new data for", filter_display, "..."))
      add_log(paste("Missing rate threshold:", input$missing_rate_threshold, "%"))

      # [修改] 使用UI中设置的missing rate threshold
      mgr <- ProteomicsDataManager$new(rv$raw_mat, rv$meta, rv$annot, input$missing_rate_threshold / 100)

      # 记录过滤前的蛋白数量
      initial_protein_count <- nrow(rv$raw_mat) - 1  # 减去第一列（ID列）

      # [修改] 使用新的process_data参数
      if (input$sel_filter_column == "none") {
        mgr$process_data(filter_col = NULL, filter_value = NULL)
      } else {
        mgr$process_data(filter_col = input$sel_filter_column, filter_value = input$sel_filter_value)
      }

      rv$manager <- mgr

      # [修改] 增强data_summary信息，显示过滤统计
      filtered_count <- initial_protein_count - nrow(mgr$se_obj)
      output$data_summary <- renderText({
        paste0("=== Data Loading Summary ===\n",
              "Samples: ", ncol(mgr$se_obj), "\n",
              "Initial Proteins: ", initial_protein_count, "\n",
              "Filtered Out (>", input$missing_rate_threshold, "% missing): ", filtered_count, "\n",
              "Remaining Proteins: ", nrow(mgr$se_obj), "\n",
              "Retention Rate: ", round(nrow(mgr$se_obj) / initial_protein_count * 100, 1), "%")
      })
      output$missing_pattern_plot <- renderPlot({ mgr$plot_missing_pattern() })

      # [新增] 初始化outlier status
      output$txt_outlier_status <- renderUI({
        HTML("<b style='color:#95a5a6;'>Please run imputation first</b>")
      })

      updateTabsetPanel(session, "subtab_wizard", selected = "Data Overview")
      add_log(paste("Data loaded successfully:", ncol(mgr$se_obj), "samples,", nrow(mgr$se_obj), "proteins"))

    }, error = function(e) {
      sendSweetAlert(session, "Error", e$message, "error")
      add_log(paste("Data loading failed:", e$message))
    })
  })
  
  # --- 4. 插补 ---
  observeEvent(input$btn_impute, {
    req(rv$manager)
    disable("btn_impute")
    method_sel <- input$sel_imp_method

    tryCatch({
      withProgress(message = paste('Running', method_sel, '...'), value = 0.2, {
        rv$manager$perform_imputation(method = method_sel, cores = 4)
        write_rds(rv$manager, file.path(rv$cache_root, "manager.rds"))
        assess_res <- rv$manager$assess_imputation(plot_type = c("structure", "density"))
        output$qc_imp_assess <- renderPlot({ assess_res$plots })
        output$qc_pca_plot <- renderPlot({ rv$manager$plot_pca(color_col = "condition") })
        output$qc_pca_outlier_plot <- renderPlot({
          res <- rv$manager$detect_outliers(method = "pca", sd_threshold = 3, remove = FALSE, do_plot = FALSE)
          res$plot
        })
        output$qc_cor_outlier_plot <- renderPlot({
          res <- rv$manager$detect_outliers(method = "correlation", sd_threshold = 3, remove = FALSE, do_plot = FALSE)
          res$plot
        })

        # [新增] 初始化outlier status显示
        output$txt_outlier_status <- renderUI({
          total_samples <- ncol(rv$manager$imputed_se)
          HTML(paste0("<b style='color:green'>Active Samples: ", total_samples, "</b> | ",
                     "<b>No samples removed yet</b>"))
        })

        # [新增] 生成初始Clean Matrix
        generate_clean_matrix(rv$manager)

        # [新增] 更新 outlier picker
        update_outlier_picker()
      })
      sendSweetAlert(session, "Success", "Imputation done & Saved.\nClean Matrix generated.", "success")
      updateTabsetPanel(session, "subtab_wizard", selected = "Clean Matrix")
      enable("btn_impute")
      add_log(paste("Imputation completed using", method_sel))
    }, error = function(e) {
      add_log(paste("Impute Failed:", e$message))
      sendSweetAlert(session, "Imputation Error", e$message, "error")
      enable("btn_impute")
    })
  })
  
  # --- 5. Outlier Detection and Removal ---

  # 辅助函数：更新 outlier picker 的选项
  update_outlier_picker <- function() {
    if (is.null(rv$manager)) return()

    # 获取完整的样本列表
    se_obj <- if (!is.null(rv$manager$imputed_se_backup)) {
      rv$manager$imputed_se_backup
    } else if (!is.null(rv$manager$imputed_se)) {
      rv$manager$imputed_se
    } else {
      return()
    }

    all_samples <- colnames(se_obj)

    # 计算当前已被移除的样本
    removed_samples <- if (!is.null(rv$manager$imputed_se_backup)) {
      setdiff(colnames(rv$manager$imputed_se_backup), colnames(rv$manager$imputed_se))
    } else {
      character(0)
    }

    updatePickerInput(session, "sel_outliers",
                      choices = all_samples,
                      selected = removed_samples)
  }

  # Auto-Detect按钮：自动检测outliers
  observeEvent(input$btn_detect_outlier, {
    req(rv$manager$imputed_se)

    tryCatch({
      # 先重置数据到完整状态（如果有backup的话）
      if(!is.null(rv$manager$imputed_se_backup)) {
        rv$manager$reset_data()
        add_log("Reset to full dataset before detection.")
      }

      # [修复] 先更新 picker 的 choices，确保下拉菜单显示所有样本
      update_outlier_picker()

      # 运行outlier检测
      res <- rv$manager$detect_outliers(method = "pca", sd_threshold = 2.5, do_plot = TRUE)

      if(length(res$outliers) > 0) {
        # 更新picker，选中检测到的outliers
        updatePickerInput(session, "sel_outliers", selected = res$outliers)
        sendSweetAlert(session, "Detected",
                      paste("Found", length(res$outliers), "potential outlier(s):\n",
                            paste(res$outliers, collapse = ", "),
                            "\n\nClick 'Confirm Update' to remove them."),
                      "warning")
        add_log(paste("Detected outliers:", paste(res$outliers, collapse = ", ")))
      } else {
        updatePickerInput(session, "sel_outliers", selected = character(0))
        sendSweetAlert(session, "Clean", "No outliers detected. Data looks good!", "success")
        add_log("No outliers detected.")
      }

      # 更新PCA plot
      output$qc_pca_plot <- renderPlot({
        rv$manager$plot_pca(color_col = "condition") +
          ggtitle("PCA - After Detection (No Changes Yet)")
      })
      output$qc_pca_outlier_plot <- renderPlot({
        res <- rv$manager$detect_outliers(method = "pca", sd_threshold = 3, remove = FALSE, do_plot = FALSE)
        res$plot
      })
      output$qc_cor_outlier_plot <- renderPlot({
        res <- rv$manager$detect_outliers(method = "correlation", sd_threshold = 3, remove = FALSE, do_plot = FALSE)
        res$plot
      })

    }, error = function(e) {
      sendSweetAlert(session, "Detection Error", paste("Error:", e$message), "error")
      add_log(paste("Outlier detection failed:", e$message))
    })
  })

  # Confirm Update按钮：应用选择的outliers删除
  observeEvent(input$btn_confirm_remove, {
    req(rv$manager, rv$cache_root)
    disable("btn_confirm_remove")

    tryCatch({
      # 先重置到完整数据集
      if(!is.null(rv$manager$imputed_se_backup)) {
        rv$manager$reset_data()
        add_log("Reset to full dataset.")
      } else {
        # 第一次操作，创建backup
        rv$manager$imputed_se_backup <- rv$manager$imputed_se
        rv$manager$meta_data_backup <- rv$manager$meta_data
        if(!is.null(rv$manager$se_obj)) {
          # 也备份se_obj（用于assess_imputation）
          # 注意：se_obj在remove_outliers时会自动同步更新
        }
        add_log("Created backup of original data.")
      }

      target_samples <- input$sel_outliers
      total_samples <- ncol(rv$manager$imputed_se)

      # 安全检查
      if (length(target_samples) > 0 && (total_samples - length(target_samples)) < 3) {
        stop("Too many samples selected! At least 3 samples are required for analysis.")
      }

      # 执行删除
      if (length(target_samples) > 0) {
        rv$manager$remove_outliers(target_samples)
        add_log(paste("Removed", length(target_samples), "sample(s):", paste(target_samples, collapse = ", ")))

        # 保存更新后的manager
        write_rds(rv$manager, file.path(rv$cache_root, "manager.rds"))
        add_log("Changes saved to cache.")

        # 更新状态显示
        output$txt_outlier_status <- renderUI({
          curr <- ncol(rv$manager$imputed_se)
          rem <- length(target_samples)
          HTML(paste0("<b style='color:green'>Active Samples: ", curr, "</b> | ",
                     "<b style='color:red'>Removed: ", rem, "</b><br>",
                     "<small style='color:#7f8c8d;'>", paste(target_samples, collapse = ", "), "</small>"))
        })

        # 更新PCA plot
        output$qc_pca_plot <- renderPlot({
          rv$manager$plot_pca(color_col = "condition") +
            ggtitle(paste("PCA - After Removing", length(target_samples), "Sample(s)"))
        })
        output$qc_pca_outlier_plot <- renderPlot({
          res <- rv$manager$detect_outliers(method = "pca", sd_threshold = 3, remove = FALSE, do_plot = FALSE)
          res$plot
        })
        output$qc_cor_outlier_plot <- renderPlot({
          res <- rv$manager$detect_outliers(method = "correlation", sd_threshold = 3, remove = FALSE, do_plot = FALSE)
          res$plot
        })

        # 可选：更新imputation QC plot
        output$qc_imp_assess <- renderPlot({
          tryCatch({
            rv$manager$assess_imputation(plot_type = c("structure", "density"))$plots
          }, error = function(e) {
            # 如果维度不匹配，只显示density plot
            rv$manager$assess_imputation(plot_type = "density")$plots
          })
        })

        # [新增] 生成Clean Matrix表格
        generate_clean_matrix(rv$manager)

        # [新增] 更新 outlier picker
        update_outlier_picker()

        sendSweetAlert(session, "Success",
                      paste("Successfully removed", length(target_samples), "sample(s).\nClean Matrix generated."),
                      "success")

      } else {
        # 没有选择任何样本，恢复到完整数据集
        add_log("No samples selected. Data restored to full set.")

        output$txt_outlier_status <- renderUI({
          HTML(paste0("<b style='color:green'>Active Samples: ", total_samples, "</b> | ",
                     "<b>No samples removed</b>"))
        })

        output$qc_pca_plot <- renderPlot({
          rv$manager$plot_pca(color_col = "condition") + ggtitle("PCA - Full Dataset")
        })
        output$qc_pca_outlier_plot <- renderPlot({
          res <- rv$manager$detect_outliers(method = "pca", sd_threshold = 3, remove = FALSE, do_plot = FALSE)
          res$plot
        })
        output$qc_cor_outlier_plot <- renderPlot({
          res <- rv$manager$detect_outliers(method = "correlation", sd_threshold = 3, remove = FALSE, do_plot = FALSE)
          res$plot
        })

        # 保存恢复后的状态
        write_rds(rv$manager, file.path(rv$cache_root, "manager.rds"))

        sendSweetAlert(session, "Restored", "Data restored to full dataset.", "info")
      }

    }, error = function(e) {
      sendSweetAlert(session, "Update Error", paste("Error:", e$message), "error")
      add_log(paste("Update failed:", e$message))
    })

    enable("btn_confirm_remove")
  })
  
  # --- 6. 分析参数 UI ---
  output$ui_analysis_params <- renderUI({
    req(rv$manager$imputed_se)
    cols <- colnames(colData(rv$manager$imputed_se))
    tagList(
      selectInput("sel_condition_col", "Condition Column:", choices = cols, selected = "condition"),
      conditionalPanel("input.analysis_mode == 'group'",
                       selectInput("sel_control", "Control:", choices = NULL),
                       selectInput("sel_case", "Case:", choices = NULL),
                       selectInput("sel_paired_col", "Paired Column (Repeated Measures):", choices = c("None", cols), selected = "None")),
      pickerInput("sel_covariates", "Covariates:", choices = cols, multiple = TRUE)
    )
  })
  
  observeEvent(c(input$sel_condition_col, rv$manager), {
    req(rv$manager$imputed_se, input$sel_condition_col)
    meta <- colData(rv$manager$imputed_se) %>% as.data.frame()
    if(input$sel_condition_col %in% colnames(meta)) {
      levels_vec <- unique(meta[[input$sel_condition_col]])
      updateSelectInput(session, "sel_control", choices = levels_vec)
      updateSelectInput(session, "sel_case", choices = levels_vec, selected = levels_vec[2])
    }
  })
  
  # --- 7. 历史记录检测与加载 ---
  get_base_name <- reactive({
    req(input$analysis_mode, input$sel_condition_col)
    if (input$analysis_mode == "group") {
      req(input$sel_case, input$sel_control)
      paste0(input$sel_case, "_vs_", input$sel_control)
    } else if (input$analysis_mode == "continuous") {
      paste0("Regression_", input$sel_condition_col)
    } else {
      "ANOVA_Multigroup"
    }
  })
  
  observe({
    req(rv$cache_root, get_base_name())
    base_name <- get_base_name()
    all_dirs <- list.dirs(rv$cache_root, full.names = FALSE, recursive = FALSE)
    matched_dirs <- all_dirs[grep(paste0("^", base_name), all_dirs)]
    
    choices_list <- list()
    if (length(matched_dirs) > 0) {
      for (d in matched_dirs) {
        param_file <- file.path(rv$cache_root, d, "analysis_params.txt")
        if (file.exists(param_file)) {
          lines <- readLines(param_file)
          time_line <- grep("Time:", lines, value=T)[1]
          pval_line <- grep("P-val:", lines, value=T)[1]
          fc_line <- grep("LogFC:", lines, value=T)[1]
          label <- paste0(d, " [", time_line, " | ", pval_line, " | ", fc_line, "]")
          choices_list[[label]] <- d
        }
      }
    }
    
    if (length(choices_list) > 0) {
      output$ui_analysis_history <- renderUI({ selectInput("sel_existing_result", "Found History:", choices = choices_list, width = "100%") })
      output$ui_btn_load_analysis <- renderUI({ actionButton("btn_load_analysis", "Load Selected", icon = icon("folder-open"), class = "btn-info", width = "100%") })
    } else {
      output$ui_analysis_history <- renderUI({ helpText("No previous results found.") })
      output$ui_btn_load_analysis <- renderUI({ NULL })
    }
  })
  
  observeEvent(input$btn_load_analysis, {
    req(input$sel_existing_result)
    target_dir <- file.path(rv$cache_root, input$sel_existing_result)
    withProgress(message = "Loading Result...", {
      rv$diff_tool <- read_rds(file.path(target_dir, "diff_tool.rds"))
      if(file.exists(file.path(target_dir, "enrich_tool.rds"))) rv$enrich_tool <- read_rds(file.path(target_dir, "enrich_tool.rds"))
      if(file.exists(file.path(target_dir, "gsva_tool.rds"))) rv$gsva_tool <- read_rds(file.path(target_dir, "gsva_tool.rds"))
      rv$loaded_params <- paste(readLines(file.path(target_dir, "analysis_params.txt")), collapse = "<br>")
    })
    sendSweetAlert(session, "Loaded", "Result loaded.", "success")
    updateNavbarPage(session, "main_nav", selected = "Visualization")
  })
  
  # --- 8. 运行新分析 (Hash) ---
  observeEvent(input$btn_run_pipeline, {
    req(rv$manager$imputed_se)
    disable("btn_run_pipeline")
    add_log("Starting Analysis Pipeline...")

    # Build comprehensive params list for hash (all parameters)
    params_list <- list(
      filter_col = input$sel_filter_column,
      filter_val = if(input$sel_filter_column != "none") input$sel_filter_value else "all",
      imp_method = input$sel_imp_method,
      mode = input$analysis_mode,
      condition_col = input$sel_condition_col,
      case = if(input$analysis_mode == "group") input$sel_case else NA,
      control = if(input$analysis_mode == "group") input$sel_control else NA,
      covariates = input$sel_covariates,
      pval_cutoff = input$param_pval,
      logfc_cutoff = if(input$analysis_mode == "group") input$param_logfc else NA,
      corr_cutoff = if(input$analysis_mode == "continuous") input$param_corr else NA,
      r2_cutoff = if(input$analysis_mode == "continuous") input$param_r2 else NA,
      enrich_pval = input$param_enrich_pval,
      use_adj_pval = input$param_use_adj,
      run_gsva = input$run_gsva
    )
    short_hash <- substr(digest::digest(params_list, algo = "md5"), 1, 6)
    final_sub_folder <- paste0(get_base_name(), "_", short_hash)
    full_res_dir <- file.path(rv$cache_root, final_sub_folder)
    if(!dir.exists(full_res_dir)) dir.create(full_res_dir, recursive = TRUE)

    # Save all parameters to analysis_params.txt
    filter_info <- if(input$sel_filter_column == "none") {
      "All samples"
    } else {
      paste0(input$sel_filter_column, " = ", input$sel_filter_value)
    }

    param_text <- c(
      paste("Analysis Time:", Sys.time()),
      paste("Hash ID:", short_hash),
      paste("Data Filter:", filter_info),
      paste("Imputation Method:", input$sel_imp_method),
      paste("Mode:", input$analysis_mode),
      paste("Condition Column:", input$sel_condition_col),
      paste("Case Group:", if(input$analysis_mode == "group") input$sel_case else "N/A"),
      paste("Control Group:", if(input$analysis_mode == "group") input$sel_control else "N/A"),
      paste("Covariates:", if(length(input$sel_covariates) > 0) paste(input$sel_covariates, collapse = ", ") else "None"),
      paste("Paired Column:", if(!is.null(input$sel_paired_col) && input$sel_paired_col != "None") input$sel_paired_col else "None"),
      paste("P-value Cutoff:", input$param_pval),
      paste("LogFC Cutoff:", if(input$analysis_mode == "group") input$param_logfc else "N/A"),
      paste("Correlation Cutoff:", if(input$analysis_mode == "continuous") input$param_corr else "N/A"),
      paste("R-squared Cutoff:", if(input$analysis_mode == "continuous") input$param_r2 else "N/A"),
      paste("Enrichment P-value Cutoff:", input$param_enrich_pval),
      paste("Use Adjusted P-value:", input$param_use_adj),
      paste("Run GSVA:", input$run_gsva)
    )
    writeLines(param_text, file.path(full_res_dir, "analysis_params.txt"))
    
    tryCatch({
      withProgress(message = 'Analyzing...', value = 0.1, {
        if (input$analysis_mode == "anova") {
          dt <- DiffExpAnalyst$new(rv$manager$imputed_se)
          dt$run_anova_analysis(condition_col = input$sel_condition_col)
          rv$diff_tool <- dt; write_rds(dt, file.path(full_res_dir, "diff_tool.rds"))

          # [修复] 先调用 get_sig_proteins 确保 sig_results 被填充
          dt$get_sig_proteins(pval_cutoff = input$param_pval, use_adjusted_pval = input$param_use_adj)

          # [修复] 安全检查 sig_results 是否存在且有数据
          if(!is.null(dt$sig_results) && nrow(dt$sig_results) > 0) {
            et <- EnrichmentAnalyst$new(); sig_genes <- dt$sig_results$Protein[dt$sig_results$adj.P.Val < input$param_pval]
            et$ora_up <- et$run_comprehensive_ora(sig_genes, rownames(rv$manager$imputed_se), input$param_enrich_pval)
            rv$enrich_tool <- et; write_rds(et, file.path(full_res_dir, "enrich_tool.rds"))
          }
          if(input$run_gsva) { gt <- ProteomicsGSVA$new(rv$manager); gt$run_gsva(dbs = c("GOBP", "GOMF", "GOCC", "KEGG", "Wiki")); rv$gsva_tool <- gt; write_rds(gt, file.path(full_res_dir, "gsva_tool.rds")) }
        } else {
          paired_col_val <- if(!is.null(input$sel_paired_col) && input$sel_paired_col != "None") input$sel_paired_col else NULL
          res_objs <- run_proteomics_pipeline(data_manager = rv$manager, analysis_type = input$analysis_mode, condition_col = input$sel_condition_col, control_group = input$sel_control, case_group = input$sel_case, continuous_col = if(input$analysis_mode == "continuous") input$sel_condition_col else NULL, covariates = if(length(input$sel_covariates)==0) NULL else input$sel_covariates, paired_col = paired_col_val, results_dir = rv$cache_root, sub_folder_name = final_sub_folder, pval_cutoff = input$param_pval, enrich_pval_cutoff = input$param_enrich_pval, logfc_cutoff = input$param_logfc, corr_cutoff = input$param_corr, r2_cutoff = input$param_r2, top_n_labels = 10, use_adj_pval_sig = input$param_use_adj, run_gsva = input$run_gsva, gsva_dbs = c("GOBP", "GOMF", "GOCC", "KEGG", "Wiki", "Reactome"))
          rv$diff_tool <- res_objs$diff_tool; rv$enrich_tool <- res_objs$enrich_tool; rv$gsva_tool <- res_objs$gsva_tool
        }
      })
      rv$loaded_params <- paste(param_text, collapse = "<br>")
      sendSweetAlert(session, "Success", "Analysis complete!", "success")
      updateNavbarPage(session, "main_nav", selected = "Visualization")
      enable("btn_run_pipeline")
    }, error = function(e) { add_log(paste("Error:", e$message)); enable("btn_run_pipeline") })
  })
  
  # --- Visualization ---
  output$txt_current_params <- renderUI({ if(is.null(rv$loaded_params)) return(HTML("<i>No analysis loaded.</i>")); wellPanel(style = "background-color: #f9fafe; font-size: 11px; padding: 5px;", HTML(paste0("<b>Current Analysis:</b><br>", rv$loaded_params))) })

  output$vis_slider_fc <- renderUI({
    req(rv$diff_tool)
    if("logFC" %in% colnames(rv$diff_tool$diff_results)) {
      fluidRow(
        column(8, sliderInput("vis_fc", "LogFC Cutoff:", min = 0, max = 2, value = input$param_logfc, step = 0.01)),
        column(4, numericInput("vis_input_fc", "Value:", value = input$param_logfc, min = 0, max = 10, step = 0.01))
      )
    } else {
      NULL
    }
  })

  # --- 同步 LogFC 滑动条和输入框 ---
  observeEvent(input$vis_fc, {
    if(!is.null(input$vis_input_fc) && !isTRUE(all.equal(input$vis_fc, input$vis_input_fc))) {
      updateNumericInput(session, "vis_input_fc", value = input$vis_fc)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$vis_input_fc, {
    if(!is.null(input$vis_fc) && !isTRUE(all.equal(input$vis_input_fc, input$vis_fc))) {
      updateSliderInput(session, "vis_fc", value = input$vis_input_fc)
    }
  }, ignoreInit = TRUE)

  # --- 同步 P-val 滑动条和输入框 ---
  observeEvent(input$vis_slider_pval, {
    if(!is.null(input$vis_input_pval) && !isTRUE(all.equal(input$vis_slider_pval, input$vis_input_pval))) {
      updateNumericInput(session, "vis_input_pval", value = input$vis_slider_pval)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$vis_input_pval, {
    if(!is.null(input$vis_slider_pval) && !isTRUE(all.equal(input$vis_input_pval, input$vis_slider_pval))) {
      updateSliderInput(session, "vis_slider_pval", value = input$vis_input_pval)
    }
  }, ignoreInit = TRUE)

  output$vis_volcano <- renderPlot({
    req(rv$diff_tool)
    if(input$analysis_mode == "anova") {
      plot(1, type="n", axes=F, xlab="", ylab=""); text(1, 1, "N/A for ANOVA", cex=1.5)
    } else if(input$analysis_mode == "continuous") {
      rv$diff_tool$plot_volcano(
        corr_cutoff = input$param_corr,
        r2_cutoff = input$param_r2,
        pval_cutoff = input$vis_slider_pval,
        use_adjusted_pval = input$param_use_adj
      )
    } else {
      rv$diff_tool$plot_volcano(
        logfc_cutoff = if(is.null(input$vis_fc)) 0.58 else input$vis_fc,
        pval_cutoff = input$vis_slider_pval,
        use_adjusted_pval = input$param_use_adj
      )
    }
  })
  output$vis_volcano_info <- renderTable({ req(rv$diff_tool, input$volcano_brush); if(input$analysis_mode == "anova") return(NULL); df <- rv$diff_tool$diff_results; is_reg <- "spearman_rho" %in% colnames(df); p_col <- if(input$param_use_adj) { if(is_reg) "adj_pval" else "adj.P.Val" } else { if(is_reg) "p_value" else "P.Value" }; plot_df <- df %>% mutate(x_val = if(is_reg) spearman_rho else logFC, log_pval = -log10(!!sym(p_col))); brushedPoints(plot_df, input$volcano_brush, xvar = "x_val", yvar = "log_pval") })
  output$vis_dep_table <- renderDT({ req(rv$diff_tool); datatable(rv$diff_tool$diff_results, options=list(scrollX=T)) })
  output$vis_heatmap <- renderPlot({
    req(rv$diff_tool, rv$manager$imputed_se)
    df <- rv$diff_tool$diff_results
    n_each <- ceiling(input$heatmap_n / 2)  # 每个方向取一半

    if(input$analysis_mode == "group") {
      # Group Comparison: 选择 top N 上调和下调 DEGs
      p_col <- if(input$param_use_adj) "adj.P.Val" else "P.Value"
      fc_cutoff <- if(is.null(input$vis_fc)) 0 else input$vis_fc
      pval_cutoff <- input$vis_slider_pval

      # 上调蛋白 (logFC > 0, 显著)
      up_prots <- df %>%
        filter(logFC > fc_cutoff & !!sym(p_col) < pval_cutoff) %>%
        arrange(!!sym(p_col)) %>%
        head(n_each) %>%
        pull(Protein)

      # 下调蛋白 (logFC < 0, 显著)
      down_prots <- df %>%
        filter(logFC < -fc_cutoff & !!sym(p_col) < pval_cutoff) %>%
        arrange(!!sym(p_col)) %>%
        head(n_each) %>%
        pull(Protein)

      sig_prots <- c(up_prots, down_prots)
      heatmap_title <- paste0("Top ", length(sig_prots), " DEGs (Up: ", length(up_prots), ", Down: ", length(down_prots), ")")

    } else if(input$analysis_mode == "continuous") {
      # Continuous Regression: 选择 top N 正相关和负相关 DEGs
      p_col <- if(input$param_use_adj) "adj_pval" else "p_value"
      corr_cutoff <- input$param_corr
      pval_cutoff <- input$vis_slider_pval

      # 正相关蛋白
      pos_prots <- df %>%
        filter(spearman_rho > corr_cutoff & !!sym(p_col) < pval_cutoff) %>%
        arrange(!!sym(p_col)) %>%
        head(n_each) %>%
        pull(Protein)

      # 负相关蛋白
      neg_prots <- df %>%
        filter(spearman_rho < -corr_cutoff & !!sym(p_col) < pval_cutoff) %>%
        arrange(!!sym(p_col)) %>%
        head(n_each) %>%
        pull(Protein)

      sig_prots <- c(pos_prots, neg_prots)
      heatmap_title <- paste0("Top ", length(sig_prots), " Correlated Proteins (Pos: ", length(pos_prots), ", Neg: ", length(neg_prots), ")")

    } else {
      # ANOVA: 使用原有方法获取显著蛋白
      sig_prots <- rv$diff_tool$get_sig_proteins(
        pval_cutoff = input$vis_slider_pval,
        use_adjusted_pval = input$param_use_adj
      )
      sig_prots <- head(sig_prots, input$heatmap_n)
      heatmap_title <- paste0("ANOVA Significant Proteins (n=", length(sig_prots), ")")
    }

    if(length(sig_prots) == 0) {
      plot(1, type="n", axes=F, xlab="", ylab="")
      text(1, 1, "No significant proteins found", cex=1.5)
      return()
    }

    # 自定义热图绘制，使用正确的标题
    se <- rv$manager$imputed_se
    mat <- assay(se)[sig_prots, , drop = FALSE]
    mat_scaled <- t(scale(t(mat)))

    # 获取annotation
    anno_df <- as.data.frame(colData(se))

    # 使用 ComplexHeatmap 替代 pheatmap，确保在 Shiny 网页中正确显示
    ha <- NULL
    if(input$sel_condition_col %in% colnames(anno_df)) {
      annot_vals <- anno_df[[input$sel_condition_col]]
      # 对 groups 排序以确保颜色分配一致
      groups <- sort(unique(as.character(annot_vals)))
      # 使用 ggsci 配色
      if(length(groups) <= 10) {
        colors <- ggsci::pal_npg()(length(groups))
      } else {
        colors <- ggsci::pal_d3("category20")(length(groups))
      }
      group_colors <- setNames(colors, groups)

      # 创建注释数据框，列名必须与 col_list 的键名一致
      annot_df_for_heatmap <- data.frame(annot_vals, row.names = colnames(mat))
      colnames(annot_df_for_heatmap) <- input$sel_condition_col

      col_list <- list()
      col_list[[input$sel_condition_col]] <- group_colors

      ha <- ComplexHeatmap::HeatmapAnnotation(
        df = annot_df_for_heatmap,
        col = col_list,
        show_annotation_name = TRUE
      )
    }

    ht <- ComplexHeatmap::Heatmap(
      mat_scaled,
      name = "Z-score",
      column_title = heatmap_title,
      top_annotation = ha,
      cluster_rows = input$heat_cluster_row,
      cluster_columns = input$heat_cluster_col,
      show_row_names = (length(sig_prots) <= 50),
      show_column_names = TRUE,
      row_names_gp = grid::gpar(fontsize = 8),
      column_names_gp = grid::gpar(fontsize = 8),
      col = circlize::colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))
    )

    grid::grid.newpage()
    ComplexHeatmap::draw(ht)
  })
  output$vis_pca <- renderPlot({ req(rv$manager); rv$manager$plot_pca(color_col = input$sel_condition_col) })
  output$vis_selector_protein <- renderUI({
    req(rv$manager)
    pickerInput("vis_sel_prot", "Select Proteins (max 9):",
                choices = rownames(rv$manager$imputed_se),
                multiple = TRUE,
                options = list(
                  `max-options` = 9,
                  `live-search` = TRUE,
                  `actions-box` = TRUE,
                  `selected-text-format` = "count > 3"
                ))
  })
  output$vis_single_prot <- renderPlot({
    req(rv$manager, input$vis_sel_prot)
    rv$manager$plot_protein_expression(input$vis_sel_prot, variable = input$sel_condition_col)
  })
  
  output$vis_ora_plot <- renderPlot({
    req(rv$diff_tool)
    if(input$analysis_mode == "anova") {
      plot(1, type="n", axes=F, xlab="", ylab=""); text(1, 1, "N/A for ANOVA", cex=1.5)
      return()
    }
    req(rv$enrich_tool)
    tryCatch({
      et <- rv$enrich_tool
      target <- if(input$vis_ora_dir=="UP") et$ora_up[[input$vis_ora_db]] else et$ora_down[[input$vis_ora_db]]
      if(!is.null(target) && inherits(target, "enrichResult") && nrow(target@result) > 0) {
        dotplot(target, showCategory=15)
      } else {
        plot(1, type="n", axes=F, xlab="", ylab=""); text(1, 1, "No enriched pathways found", cex=1.5)
      }
    }, error = function(e) {
      plot(1, type="n", axes=F, xlab="", ylab=""); text(1, 1, paste("Error:", e$message), cex=1.2, col="red")
    })
  })
  output$vis_ora_table <- renderDT({
    req(rv$diff_tool)
    if(input$analysis_mode == "anova") return(NULL)
    req(rv$enrich_tool)
    et <- rv$enrich_tool
    target <- if(input$vis_ora_dir=="UP") et$ora_up[[input$vis_ora_db]] else et$ora_down[[input$vis_ora_db]]
    if(!is.null(target)) datatable(as.data.frame(target), options=list(scrollX=T))
  })
  output$ui_vis_gsea_pathway <- renderUI({
    req(rv$diff_tool)
    if(input$analysis_mode == "anova") return(NULL)
    req(rv$enrich_tool)
    res <- rv$enrich_tool$gsea_res[[input$vis_gsea_db]]
    if(!is.null(res) && nrow(as.data.frame(res)) > 0) selectInput("vis_sel_gsea_path", "Pathway:", choices = res$ID) else NULL
  })
  output$vis_gsea_dotplot <- renderPlot({
    req(rv$diff_tool)
    if(input$analysis_mode == "anova") {
      plot(1, type="n", axes=F, xlab="", ylab=""); text(1, 1, "N/A for ANOVA", cex=1.5)
      return()
    }
    req(rv$enrich_tool)
    tryCatch({
      res <- rv$enrich_tool$gsea_res[[input$vis_gsea_db]]
      if(!is.null(res) && inherits(res, "gseaResult") && nrow(res@result) > 0) {
        dotplot(res, split=".sign") + facet_grid(.~.sign)
      } else {
        plot(1, type="n", axes=F, xlab="", ylab=""); text(1, 1, "No GSEA results", cex=1.5)
      }
    }, error = function(e) {
      plot(1, type="n", axes=F, xlab="", ylab=""); text(1, 1, paste("Error:", e$message), cex=1.2, col="red")
    })
  })
  output$vis_gsea_curve <- renderPlot({
    req(rv$diff_tool)
    if(input$analysis_mode == "anova") {
      plot(1, type="n", axes=F, xlab="", ylab=""); text(1, 1, "N/A for ANOVA", cex=1.5)
      return()
    }
    req(rv$enrich_tool, input$vis_sel_gsea_path)
    res <- rv$enrich_tool$gsea_res[[input$vis_gsea_db]]
    if(!is.null(res)) enrichplot::gseaplot2(res, input$vis_sel_gsea_path)
  })
  output$vis_gsva_plot <- renderPlot({
    req(rv$diff_tool)
    if(input$analysis_mode == "anova") {
      plot(1, type="n", axes=F, xlab="", ylab=""); text(1, 1, "N/A for ANOVA", cex=1.5)
      return()
    }
    req(rv$gsva_tool)
    gt <- rv$gsva_tool
    db <- input$vis_gsva_db
    if(input$vis_gsva_plot_type == "Volcano") {
      gt$plot_pathway_volcano(db, pval_cutoff = input$vis_slider_pval)
    } else if (input$vis_gsva_plot_type == "Heatmap") {
      # ComplexHeatmap 需要用 draw() 并通过 grid.grabExpr 捕获
      ht <- gt$plot_pathway_heatmap(db, group_col = input$sel_condition_col)
      grid::grid.newpage()
      ComplexHeatmap::draw(ht)
    } else {
      gt$plot_pathway_gsea_bar(db)
    }
  })
}

shinyApp(ui, server)