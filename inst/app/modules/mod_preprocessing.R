# =============================================================================
# Module: Preprocessing (Data Import, QC, Imputation, Outlier Removal)
# =============================================================================

preprocessingUI <- function(id) {
  ns <- NS(id)
  tabPanel("Preprocessing", icon = icon("vial"),
    sidebarLayout(
      sidebarPanel(width = 3,
        h4(icon("upload"), "Data Import"),
        shinyFilesButton(ns("upload_file"), "Select Excel File",
                        title = "Please select an Excel file (.xlsx)",
                        multiple = FALSE,
                        buttonType = "primary",
                        class = "btn-block"),
        verbatimTextOutput(ns("selected_file_path"), placeholder = TRUE),
        helpText("Sheet1: Matrix; Sheet2: Metadata; Sheet3: Annot"),

        uiOutput(ns("ui_filter_column_selector")),
        uiOutput(ns("ui_filter_value_selector")),

        hr(),
        h4(icon("filter"), "Quality Control"),
        sliderInput(ns("missing_rate_threshold"), "Max Missing Rate per Protein:",
                   min = 0, max = 100, value = 70, step = 5,
                   post = "%", ticks = TRUE),
        helpText("Proteins with missing values above this threshold will be removed."),

        uiOutput(ns("ui_manager_cache_alert")),
        checkboxInput(ns("use_manager_cache"), "Load Processed Data (if avail.)", TRUE),

        actionButton(ns("btn_load_process"), "Load Data", icon = icon("play"),
                    class = "btn-primary", width = "100%"),

        hr(),
        h4(icon("magic"), "Imputation"),
        selectInput(ns("sel_imp_method"), "Method:",
                    choices = c("missForest (Accurate)" = "missForest",
                                "MinProb" = "MinProb", "QRILC" = "QRILC",
                                "MinDet" = "MinDet", "KNN" = "knn",
                                "BPCA" = "bpca", "Zero" = "zero", "Minimum" = "min"),
                    selected = "missForest"),
        conditionalPanel(
          condition = paste0("input['", ns("sel_imp_method"), "'] == 'missForest'"),
          numericInput(ns("imp_ntree"), "Number of Trees:", value = 100, min = 10, max = 500, step = 10),
          numericInput(ns("imp_maxiter"), "Max Iterations:", value = 10, min = 1, max = 20, step = 1),
          helpText("Fewer trees/iterations = faster but less accurate.")
        ),
        actionButton(ns("btn_impute"), "Run Imputation", icon = icon("bolt"),
                    class = "btn-warning", width = "100%"),

        hr(),
        h4(icon("search-minus"), "Outlier Removal"),
        helpText("Step 1: Auto-detect or manually select samples."),
        actionButton(ns("btn_detect_outlier"), "Auto-Detect (PCA)", width = "100%"),
        br(), br(),
        pickerInput(ns("sel_outliers"), "Select Samples to Remove:",
                    choices = NULL, selected = NULL, multiple = TRUE,
                    options = list(`actions-box` = TRUE, `live-search` = TRUE,
                                   `selected-text-format` = "count > 3",
                                   `none-selected-text` = "No samples loaded yet")),
        br(),
        helpText("Step 2: Click Confirm to apply changes."),
        actionButton(ns("btn_confirm_remove"), "Confirm Update",
                    class = "btn-danger", width = "100%"),
        hr(),
        htmlOutput(ns("txt_outlier_status"))
      ),
      mainPanel(width = 9,
        tabsetPanel(id = ns("subtab_wizard"),
          tabPanel("Data Overview", icon = icon("table"), br(),
                  verbatimTextOutput(ns("data_summary")), hr(),
                  h4("Initial Missing Values Pattern"),
                  plotOutput(ns("missing_pattern_plot"), height = "500px")),
          tabPanel("Imputation QC", icon = icon("check-double"), br(),
                  h4("Imputation Quality Assessment"),
                  plotOutput(ns("qc_imp_assess"), height = "600px")),
          tabPanel("PCA (Outlier Check)", icon = icon("crosshairs"), br(),
                  h4("PCA Overview"),
                  fluidRow(
                    column(6, plotOutput(ns("qc_pca_plot"), height = "450px")),
                    column(6, plotOutput(ns("qc_pca_outlier_plot"), height = "450px"))
                  ),
                  hr(),
                  h4("Sample Correlation (Outlier Detection)"),
                  plotOutput(ns("qc_cor_outlier_plot"), height = "400px")),
          tabPanel("Clean Matrix", icon = icon("th"), br(),
                  fluidRow(
                    column(12,
                          h4("Filtered & Imputed Protein Expression Matrix"),
                          helpText("Red values indicate imputed data."),
                          downloadButton(ns("download_clean_matrix"), "Download Matrix (Excel)",
                                        class = "btn-success"),
                          hr(),
                          DTOutput(ns("clean_matrix_table"))
                    )
                  ))
        )
      )
    )
  )
}

preprocessingServer <- function(id, rv, add_log) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # --- shinyFiles: File selector initialization ---
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

    shinyFileChoose(input, "upload_file", roots = get_volumes, session = session,
                    filetypes = c("xlsx"))

    output$selected_file_path <- renderText({
      if (is.null(rv$selected_file_path)) "No file selected"
      else basename(rv$selected_file_path)
    })

    # --- Input Validation Helper ---
    validate_excel_input <- function(file_path) {
      errors <- character(0)
      if (!file.exists(file_path)) return("File not found.")
      sheets <- readxl::excel_sheets(file_path)
      if (length(sheets) < 3) {
        return(paste0("Excel file must have at least 3 sheets (found ", length(sheets), ")."))
      }
      mat <- tryCatch(readxl::read_excel(file_path, sheet = 1), error = function(e) NULL)
      meta <- tryCatch(readxl::read_excel(file_path, sheet = 2), error = function(e) NULL)
      annot <- tryCatch(readxl::read_excel(file_path, sheet = 3), error = function(e) NULL)
      if (is.null(mat)) errors <- c(errors, "Cannot read Sheet 1 (Matrix).")
      if (is.null(meta)) errors <- c(errors, "Cannot read Sheet 2 (Metadata).")
      if (is.null(annot)) errors <- c(errors, "Cannot read Sheet 3 (Annotation).")
      if (length(errors) > 0) return(paste(errors, collapse = "\n"))
      if (ncol(mat) < 2) errors <- c(errors, "Matrix sheet must have at least 2 columns.")
      if (!"label" %in% colnames(meta)) errors <- c(errors, "Metadata must contain 'label' column.")
      if (ncol(annot) < 2) errors <- c(errors, "Annotation must have at least 2 columns.")
      if (length(errors) > 0) return(paste(errors, collapse = "\n"))
      return(NULL)
    }

    # --- Clean Matrix helper ---
    generate_clean_matrix <- function(manager) {
      req(manager, manager$imputed_se)
      tryCatch({
        imp_mat <- SummarizedExperiment::assay(manager$imputed_se)
        missing_mask <- manager$missing_mask
        if (!is.null(missing_mask) && nrow(missing_mask) == nrow(imp_mat)) {
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
          missing_mask <- matrix(FALSE, nrow = nrow(imp_mat), ncol = ncol(imp_mat))
        }

        protein_ids <- SummarizedExperiment::rowData(manager$imputed_se)$ID
        gene_names <- SummarizedExperiment::rowData(manager$imputed_se)$name
        if (is.null(protein_ids)) protein_ids <- rownames(imp_mat)
        if (is.null(gene_names)) gene_names <- rownames(imp_mat)

        col_meta <- as.data.frame(SummarizedExperiment::colData(manager$imputed_se))
        if ("label" %in% colnames(col_meta)) {
          label_mapping <- col_meta$label
          names(label_mapping) <- rownames(col_meta)
          new_colnames <- sapply(colnames(imp_mat), function(x) {
            if (x %in% names(label_mapping)) label_mapping[x] else x
          })
          colnames(imp_mat) <- new_colnames
          colnames(missing_mask) <- new_colnames
        }

        clean_df <- data.frame(
          ProteinID = protein_ids, GeneName = gene_names,
          imp_mat, check.names = FALSE, stringsAsFactors = FALSE
        )

        imputed_count <- sum(missing_mask)
        imputed_pct <- round(imputed_count / (nrow(missing_mask) * ncol(missing_mask)) * 100, 2)

        output$clean_matrix_table <- renderDT({
          datatable(clean_df,
            options = list(pageLength = 25, scrollX = TRUE, scrollY = "600px",
                          dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel'),
                          fixedColumns = list(leftColumns = 2)),
            rownames = FALSE, class = 'cell-border stripe hover',
            extensions = c('Buttons', 'FixedColumns'),
            caption = htmltools::tags$caption(
              style = 'caption-side: top; text-align: left; color: #666; font-size: 12px;',
              paste0('Imputed values: ', imputed_count, ' (', imputed_pct, '%).')
            )
          )
        })

        output$download_clean_matrix <- downloadHandler(
          filename = function() paste0("clean_matrix_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"),
          content = function(file) {
            wb <- openxlsx::createWorkbook()
            openxlsx::addWorksheet(wb, "Clean Matrix")
            openxlsx::writeData(wb, "Clean Matrix", clean_df, startRow = 1, startCol = 1)
            idx <- which(missing_mask, arr.ind = TRUE)
            if (nrow(idx) > 0) {
              style <- openxlsx::createStyle(fgFill = "#ffcccc", fontColour = "#cc0000", textDecoration = "bold")
              openxlsx::addStyle(wb, "Clean Matrix", style, rows = idx[, 1] + 1, cols = idx[, 2] + 2, gridExpand = FALSE)
            }
            openxlsx::setColWidths(wb, "Clean Matrix", cols = 1:ncol(clean_df), widths = "auto")
            openxlsx::freezePane(wb, "Clean Matrix", firstActiveRow = 2, firstActiveCol = 3)
            openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
          }
        )

        add_log("Clean matrix generated successfully.")
        updateTabsetPanel(session, "subtab_wizard", selected = "Clean Matrix")
      }, error = function(e) {
        add_log(paste("Failed to generate clean matrix:", e$message))
      })
    }

    # --- File upload observer ---
    observeEvent(input$upload_file, {
      if (is.integer(input$upload_file)) return()
      file_selected <- parseFilePaths(get_volumes(), input$upload_file)
      req(nrow(file_selected) > 0)
      path <- as.character(file_selected$datapath)
      rv$selected_file_path <- path
      file_mem$last_dir <- dirname(path)

      tryCatch({
        add_log(paste("File selected:", path))
        validation_error <- validate_excel_input(path)
        if (!is.null(validation_error)) {
          showNotification(validation_error, type = "error", duration = 10)
          add_log(paste("Validation failed:", validation_error))
          return()
        }
        temp_meta <- readxl::read_excel(path, sheet = 2)
        meta_cols <- colnames(temp_meta)
        output$ui_filter_column_selector <- renderUI({
          selectInput(ns("sel_filter_column"), "Filter by Column:",
                      choices = c("None (All samples)" = "none", meta_cols),
                      selected = if("tissue" %in% meta_cols) "tissue" else "none")
        })
        rv$raw_mat <- readxl::read_excel(path, sheet = 1)
        rv$meta <- temp_meta
        rv$annot <- readxl::read_excel(path, sheet = 3)
      }, error = function(e) {
        sendSweetAlert(session, "File Error", e$message, "error")
      })
    })

    # Dynamic filter value selector
    observeEvent(input$sel_filter_column, {
      req(rv$meta, input$sel_filter_column)
      if (input$sel_filter_column == "none") {
        output$ui_filter_value_selector <- renderUI({
          helpText("All samples will be used.", style = "color: #27ae60; font-weight: bold;")
        })
      } else {
        values <- unique(rv$meta[[input$sel_filter_column]])
        output$ui_filter_value_selector <- renderUI({
          pickerInput(ns("sel_filter_value"), paste0("Select ", input$sel_filter_column, ":"),
                      choices = values, selected = values[1],
                      options = list(`live-search` = TRUE))
        })
      }
    })

    # Dynamic cache path observer
    observe({
      req(rv$selected_file_path, input$sel_filter_column)
      tryCatch({
        parent_dir <- dirname(rv$selected_file_path)
        file_id <- tools::file_path_sans_ext(basename(rv$selected_file_path))
        # Sanitize file_id to remove special characters that break directory creation
        file_id <- gsub("[^a-zA-Z0-9._-]", "_", file_id)
        if (input$sel_filter_column == "none") {
          safe_filter <- "all_samples"
        } else {
          req(input$sel_filter_value)
          safe_filter <- paste0(gsub("[^a-zA-Z0-9]", "_", input$sel_filter_column), "_",
                               gsub("[^a-zA-Z0-9]", "_", input$sel_filter_value))
        }
        rv$cache_root <- file.path(parent_dir, "results_web_session", file_id, safe_filter)
        if(!dir.exists(rv$cache_root)) {
          dir.create(rv$cache_root, recursive = TRUE, showWarnings = TRUE)
          if (!dir.exists(rv$cache_root)) {
            add_log(paste("Warning: Could not create results directory:", rv$cache_root))
          }
        }
      }, error = function(e) {
        add_log(paste("Error setting up results directory:", e$message))
      })

      manager_path <- file.path(rv$cache_root, "manager.rds")
      if (file.exists(manager_path)) {
        filter_display <- if (input$sel_filter_column == "none") "All samples"
                          else paste0(input$sel_filter_column, " = ", input$sel_filter_value)
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

    # --- Load Data ---
    observeEvent(input$btn_load_process, {
      req(rv$raw_mat, input$sel_filter_column, rv$cache_root)
      manager_path <- file.path(rv$cache_root, "manager.rds")

      if (input$use_manager_cache && file.exists(manager_path)) {
        filter_display <- if (input$sel_filter_column == "none") "All samples" else paste0(input$sel_filter_column, " = ", input$sel_filter_value)
        add_log(paste("Loading cached data for", filter_display, "..."))
        withProgress(message = "Loading Pre-processed Data...", {
          rv$manager <- readr::read_rds(manager_path)
        })
        sendSweetAlert(session, "Loaded", paste("Loaded data for", filter_display), "success")
        output$data_summary <- renderText({ paste0("Samples: ", ncol(rv$manager$se_obj), "\nProteins: ", nrow(rv$manager$se_obj)) })
        output$missing_pattern_plot <- renderPlot({ rv$manager$plot_missing_pattern() })

        if (!is.null(rv$manager$imputed_se)) {
          actual_method <- if (!is.null(rv$manager$imputation_method)) rv$manager$imputation_method else "missForest"
          updateSelectInput(session, "sel_imp_method", selected = actual_method)

          if (!is.null(rv$manager$imputed_se_backup)) {
            removed_samples <- setdiff(colnames(rv$manager$imputed_se_backup), colnames(rv$manager$imputed_se))
            removed_count <- length(removed_samples)
            current_count <- ncol(rv$manager$imputed_se)
            output$txt_outlier_status <- renderUI({
              if (removed_count > 0) {
                HTML(paste0("<b style='color:green'>Active Samples: ", current_count, "</b> | ",
                           "<b style='color:red'>Removed: ", removed_count, "</b><br>",
                           "<small style='color:#7f8c8d;'>", paste(removed_samples, collapse = ", "), "</small>"))
              } else {
                HTML(paste0("<b style='color:green'>Active Samples: ", current_count, "</b> | <b>No samples removed</b>"))
              }
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
          generate_clean_matrix(rv$manager)
          update_outlier_picker()
          updateTabsetPanel(session, "subtab_wizard", selected = "Clean Matrix")
        }
        return()
      }

      tryCatch({
        filter_display <- if (input$sel_filter_column == "none") "All samples" else paste0(input$sel_filter_column, " = ", input$sel_filter_value)
        add_log(paste("Processing new data for", filter_display, "..."))
        mgr <- ProteomicsDataManager$new(rv$raw_mat, rv$meta, rv$annot, input$missing_rate_threshold / 100)
        initial_protein_count <- nrow(rv$raw_mat) - 1
        if (input$sel_filter_column == "none") {
          mgr$process_data(filter_col = NULL, filter_value = NULL)
        } else {
          mgr$process_data(filter_col = input$sel_filter_column, filter_value = input$sel_filter_value)
        }
        rv$manager <- mgr
        filtered_count <- initial_protein_count - nrow(mgr$se_obj)
        output$data_summary <- renderText({
          paste0("=== Data Loading Summary ===\n",
                "Samples: ", ncol(mgr$se_obj), "\nInitial Proteins: ", initial_protein_count,
                "\nFiltered Out (>", input$missing_rate_threshold, "% missing): ", filtered_count,
                "\nRemaining Proteins: ", nrow(mgr$se_obj))
        })
        output$missing_pattern_plot <- renderPlot({ mgr$plot_missing_pattern() })
        output$txt_outlier_status <- renderUI({ HTML("<b style='color:#95a5a6;'>Please run imputation first</b>") })
        updateTabsetPanel(session, "subtab_wizard", selected = "Data Overview")
        add_log(paste("Data loaded:", ncol(mgr$se_obj), "samples,", nrow(mgr$se_obj), "proteins"))
      }, error = function(e) {
        sendSweetAlert(session, "Error", e$message, "error")
        add_log(paste("Data loading failed:", e$message))
      })
    })

    # --- Imputation ---
    observeEvent(input$btn_impute, {
      req(rv$manager)
      disable("btn_impute")
      method_sel <- input$sel_imp_method
      tryCatch({
        withProgress(message = paste('Running', method_sel, '...'), value = 0.2, {
          rv$manager$perform_imputation(
            method = method_sel, cores = 4,
            ntree = if (!is.null(input$imp_ntree)) input$imp_ntree else 100,
            maxiter = if (!is.null(input$imp_maxiter)) input$imp_maxiter else 10
          )
          readr::write_rds(rv$manager, file.path(rv$cache_root, "manager.rds"))
          assess_res <- rv$manager$assess_imputation(plot_type = c("structure", "density"))
          output$qc_imp_assess <- renderPlot({ assess_res$plots })
          output$qc_pca_plot <- renderPlot({ rv$manager$plot_pca(color_col = "condition") })
          output$qc_pca_outlier_plot <- renderPlot({
            res <- rv$manager$detect_outliers(method = "pca", sd_threshold = 3, remove = FALSE, do_plot = FALSE); res$plot
          })
          output$qc_cor_outlier_plot <- renderPlot({
            res <- rv$manager$detect_outliers(method = "correlation", sd_threshold = 3, remove = FALSE, do_plot = FALSE); res$plot
          })
          output$txt_outlier_status <- renderUI({
            HTML(paste0("<b style='color:green'>Active Samples: ", ncol(rv$manager$imputed_se), "</b> | <b>No samples removed yet</b>"))
          })
          generate_clean_matrix(rv$manager)
          update_outlier_picker()
        })
        sendSweetAlert(session, "Success", "Imputation done & Saved.", "success")
        updateTabsetPanel(session, "subtab_wizard", selected = "Clean Matrix")
        enable("btn_impute")
        add_log(paste("Imputation completed using", method_sel))
      }, error = function(e) {
        add_log(paste("Impute Failed:", e$message))
        sendSweetAlert(session, "Imputation Error", e$message, "error")
        enable("btn_impute")
      })
    })

    # --- Outlier Detection ---
    update_outlier_picker <- function() {
      if (is.null(rv$manager)) return()
      se_obj <- if (!is.null(rv$manager$imputed_se_backup)) rv$manager$imputed_se_backup
                else if (!is.null(rv$manager$imputed_se)) rv$manager$imputed_se
                else return()
      all_samples <- colnames(se_obj)
      removed_samples <- if (!is.null(rv$manager$imputed_se_backup))
        setdiff(colnames(rv$manager$imputed_se_backup), colnames(rv$manager$imputed_se))
      else character(0)
      updatePickerInput(session, "sel_outliers", choices = all_samples, selected = removed_samples)
    }

    observeEvent(input$btn_detect_outlier, {
      req(rv$manager$imputed_se)
      tryCatch({
        if(!is.null(rv$manager$imputed_se_backup)) {
          rv$manager$reset_data()
          add_log("Reset to full dataset before detection.")
        }
        update_outlier_picker()
        res <- rv$manager$detect_outliers(method = "pca", sd_threshold = 2.5, do_plot = TRUE)
        if(length(res$outliers) > 0) {
          updatePickerInput(session, "sel_outliers", selected = res$outliers)
          sendSweetAlert(session, "Detected",
                        paste("Found", length(res$outliers), "potential outlier(s):", paste(res$outliers, collapse = ", ")),
                        "warning")
          add_log(paste("Detected outliers:", paste(res$outliers, collapse = ", ")))
        } else {
          updatePickerInput(session, "sel_outliers", selected = character(0))
          sendSweetAlert(session, "Clean", "No outliers detected.", "success")
        }
        output$qc_pca_plot <- renderPlot({ rv$manager$plot_pca(color_col = "condition") })
        output$qc_pca_outlier_plot <- renderPlot({
          res <- rv$manager$detect_outliers(method = "pca", sd_threshold = 3, remove = FALSE, do_plot = FALSE); res$plot
        })
        output$qc_cor_outlier_plot <- renderPlot({
          res <- rv$manager$detect_outliers(method = "correlation", sd_threshold = 3, remove = FALSE, do_plot = FALSE); res$plot
        })
      }, error = function(e) {
        sendSweetAlert(session, "Detection Error", paste("Error:", e$message), "error")
      })
    })

    observeEvent(input$btn_confirm_remove, {
      req(rv$manager, rv$cache_root)
      disable("btn_confirm_remove")
      tryCatch({
        if(!is.null(rv$manager$imputed_se_backup)) {
          rv$manager$reset_data()
        } else {
          rv$manager$imputed_se_backup <- rv$manager$imputed_se
          rv$manager$meta_data_backup <- rv$manager$meta_data
        }
        target_samples <- input$sel_outliers
        total_samples <- ncol(rv$manager$imputed_se)
        if (length(target_samples) > 0 && (total_samples - length(target_samples)) < 3) {
          stop("Too many samples selected! At least 3 required.")
        }
        if (length(target_samples) > 0) {
          rv$manager$remove_outliers(target_samples)
          add_log(paste("Removed", length(target_samples), "sample(s)"))
          readr::write_rds(rv$manager, file.path(rv$cache_root, "manager.rds"))
          output$txt_outlier_status <- renderUI({
            HTML(paste0("<b style='color:green'>Active: ", ncol(rv$manager$imputed_se), "</b> | ",
                       "<b style='color:red'>Removed: ", length(target_samples), "</b>"))
          })
        } else {
          output$txt_outlier_status <- renderUI({
            HTML(paste0("<b style='color:green'>Active: ", total_samples, "</b> | <b>No samples removed</b>"))
          })
          readr::write_rds(rv$manager, file.path(rv$cache_root, "manager.rds"))
        }
        output$qc_pca_plot <- renderPlot({ rv$manager$plot_pca(color_col = "condition") })
        output$qc_imp_assess <- renderPlot({
          tryCatch(rv$manager$assess_imputation(plot_type = c("structure", "density"))$plots,
                   error = function(e) rv$manager$assess_imputation(plot_type = "density")$plots)
        })
        generate_clean_matrix(rv$manager)
        update_outlier_picker()
        sendSweetAlert(session, "Success", "Update applied.", "success")
      }, error = function(e) {
        sendSweetAlert(session, "Error", e$message, "error")
      })
      enable("btn_confirm_remove")
    })

    # Return reactive inputs that other modules need
    return(list(
      sel_filter_column = reactive(input$sel_filter_column),
      sel_filter_value = reactive(input$sel_filter_value),
      sel_imp_method = reactive(input$sel_imp_method),
      missing_rate_threshold = reactive(input$missing_rate_threshold)
    ))
  })
}
