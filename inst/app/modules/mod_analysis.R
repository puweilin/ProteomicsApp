# =============================================================================
# Module: Analysis Config (Pipeline Execution, Result History)
# =============================================================================

analysisUI <- function(id) {
  ns <- NS(id)
  tabPanel("Analysis Config", icon = icon("cogs"),
    sidebarLayout(
      sidebarPanel(width = 3,
        h4("Analysis Settings"),
        radioButtons(ns("analysis_mode"), "Mode:",
                     choices = c("Group Comparison" = "group",
                                 "ANOVA (Multi-group)" = "anova",
                                 "Continuous Regression" = "continuous")),
        uiOutput(ns("ui_analysis_params")),

        hr(),
        h4("Thresholds & Params"),
        numericInput(ns("param_pval"), "P-val Cutoff:", 0.05, step = 0.01),
        conditionalPanel(
          condition = paste0("input['", ns("analysis_mode"), "'] == 'group'"),
          numericInput(ns("param_logfc"), "LogFC Cutoff:", 0.263, step = 0.1)
        ),
        conditionalPanel(
          condition = paste0("input['", ns("analysis_mode"), "'] == 'continuous'"),
          numericInput(ns("param_corr"), "Correlation Cutoff (Rho):", 0.5, step = 0.1, min = 0, max = 1),
          numericInput(ns("param_r2"), "R-squared Cutoff:", 0.5, step = 0.1, min = 0, max = 1)
        ),
        numericInput(ns("param_enrich_pval"), "Enrich P-val:", 0.05, step = 0.01),
        checkboxInput(ns("param_use_adj"), "Use Adjusted P-val?", TRUE),
        checkboxInput(ns("run_gsva"), "Run GSVA? (Time Consuming)", TRUE),

        hr(),
        h4("Gene Set Database"),
        uiOutput(ns("ui_cache_status")),
        uiOutput(ns("ui_cache_update_banner")),
        actionButton(ns("btn_rebuild_cache"), "Rebuild Cache",
                     icon = icon("refresh"), class = "btn-sm btn-default",
                     width = "100%"),

        hr(),
        actionButton(ns("btn_run_pipeline"), "Run New Analysis",
                     icon = icon("rocket"), class = "btn-danger", width = "100%"),

        hr(),
        h4("Result Management"),
        uiOutput(ns("ui_analysis_history")),
        uiOutput(ns("ui_btn_load_analysis")),

        hr(),
        h4(icon("save"), "Session Management"),
        downloadButton(ns("btn_save_session"), "Save Session (.rds)",
                       class = "btn-info", width = "100%"),
        br(), br(),
        fileInput(ns("upload_session"), "Restore Session:", accept = ".rds"),
        helpText("Save/restore complete analysis state including all results.")
      ),
      mainPanel(width = 9,
        h4("Analysis Log Console"),
        verbatimTextOutput(ns("analysis_log"))
      )
    )
  )
}

analysisServer <- function(id, rv, add_log) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # --- Log display ---
    output$analysis_log <- renderText({ paste(rv$logs, collapse = "\n") })

    # --- Gene Set Cache: auto-build on first launch ---
    observe({
      if (!is.null(rv$cache_manager)) return()
      rv$cache_manager <- GeneSetCacheManager$new()
      if (!rv$cache_manager$is_cache_exists()) {
        showModal(modalDialog(
          title = "Building Gene Set Cache",
          "Building gene set cache for first use. KEGG requires internet access. This may take a few minutes...",
          footer = NULL
        ))
        tryCatch({
          rv$cache_manager$build_cache()
          removeModal()
          showNotification("Gene set cache ready!", type = "message", duration = 5)
          add_log("Gene set cache built successfully.")
        }, error = function(e) {
          removeModal()
          showNotification(
            paste("Cache build failed:", e$message, "- Enrichment will use online APIs as fallback."),
            type = "warning", duration = 10
          )
          add_log(paste("Gene set cache build failed:", e$message))
        })
      } else {
        rv$cache_update_info <- rv$cache_manager$check_for_updates()
        add_log("Gene set cache loaded.")
      }
    })

    # --- Cache status display ---
    output$ui_cache_status <- renderUI({
      req(rv$cache_manager)
      info <- rv$cache_manager$get_cache_info()
      tags$small(style = "color:grey;", HTML(info))
    })

    # --- Cache update banner ---
    output$ui_cache_update_banner <- renderUI({
      req(rv$cache_update_info)
      if (rv$cache_update_info$needs_refresh) {
        reason <- if (rv$cache_update_info$has_update) {
          "Annotation database update available"
        } else {
          paste0("Cache expired (", rv$cache_update_info$cache_age_days, " days old)")
        }
        div(class = "alert alert-warning", style = "padding: 8px; margin-top: 5px;",
          icon("exclamation-triangle"), " ", reason,
          actionButton(ns("btn_update_cache"), "Update Now",
                       class = "btn-sm btn-warning pull-right")
        )
      }
    })

    # --- Cache update button ---
    observeEvent(input$btn_update_cache, {
      req(rv$cache_manager)
      withProgress(message = "Updating gene set cache...", value = 0, {
        rv$cache_manager$build_cache(progress_callback = function(db, i, total) {
          incProgress(1 / total, detail = db)
        })
      })
      rv$cache_update_info <- rv$cache_manager$check_for_updates()
      showNotification("Gene set cache updated!", type = "message", duration = 5)
      add_log("Gene set cache updated.")
    })

    # --- Cache rebuild button ---
    observeEvent(input$btn_rebuild_cache, {
      req(rv$cache_manager)
      showModal(modalDialog(
        title = "Rebuild Gene Set Cache",
        "This will delete the existing cache and rebuild from scratch. Continue?",
        footer = tagList(
          modalButton("Cancel"),
          actionButton(ns("btn_rebuild_confirm"), "Rebuild", class = "btn-danger")
        )
      ))
    })

    observeEvent(input$btn_rebuild_confirm, {
      removeModal()
      req(rv$cache_manager)
      if (dir.exists(rv$cache_manager$cache_dir)) {
        unlink(rv$cache_manager$cache_dir, recursive = TRUE)
        dir.create(rv$cache_manager$cache_dir, recursive = TRUE)
        rv$cache_manager$cache_metadata <- NULL
        add_log("Gene set cache cleared.")
      }
      withProgress(message = "Rebuilding gene set cache...", value = 0, {
        rv$cache_manager$build_cache(progress_callback = function(db, i, total) {
          incProgress(1 / total, detail = db)
        })
      })
      rv$cache_update_info <- rv$cache_manager$check_for_updates()
      showNotification("Gene set cache rebuilt!", type = "message", duration = 5)
      add_log("Gene set cache rebuilt from scratch.")
    })

    # --- Analysis parameters UI ---
    output$ui_analysis_params <- renderUI({
      req(rv$manager$imputed_se)
      cols <- colnames(SummarizedExperiment::colData(rv$manager$imputed_se))
      tagList(
        selectInput(ns("sel_condition_col"), "Condition Column:",
                    choices = cols, selected = "condition"),
        conditionalPanel(
          condition = paste0("input['", ns("analysis_mode"), "'] == 'group'"),
          selectInput(ns("sel_control"), "Control:", choices = NULL),
          selectInput(ns("sel_case"), "Case:", choices = NULL),
          selectInput(ns("sel_paired_col"), "Paired Column (Repeated Measures):",
                      choices = c("None", cols), selected = "None")
        ),
        pickerInput(ns("sel_covariates"), "Covariates:", choices = cols, multiple = TRUE)
      )
    })

    # Update case/control choices when condition column changes
    observeEvent(c(input$sel_condition_col, rv$manager), {
      req(rv$manager$imputed_se, input$sel_condition_col)
      meta <- as.data.frame(SummarizedExperiment::colData(rv$manager$imputed_se))
      if (input$sel_condition_col %in% colnames(meta)) {
        levels_vec <- unique(meta[[input$sel_condition_col]])
        updateSelectInput(session, "sel_control", choices = levels_vec)
        updateSelectInput(session, "sel_case", choices = levels_vec, selected = levels_vec[2])
      }
    })

    # --- Result history detection and loading ---
    get_base_name <- reactive({
      req(input$analysis_mode, input$sel_condition_col)
      raw_name <- if (input$analysis_mode == "group") {
        req(input$sel_case, input$sel_control)
        paste0(input$sel_case, "_vs_", input$sel_control)
      } else if (input$analysis_mode == "continuous") {
        paste0("Regression_", input$sel_condition_col)
      } else {
        "ANOVA_Multigroup"
      }
      # Sanitize to avoid special characters in directory names
      gsub("[^a-zA-Z0-9._-]", "_", raw_name)
    })

    observe({
      req(rv$cache_root, get_base_name())
      base_name <- get_base_name()
      all_dirs <- list.dirs(rv$cache_root, full.names = FALSE, recursive = FALSE)
      matched_dirs <- all_dirs[startsWith(all_dirs, base_name)]

      choices_list <- list()
      if (length(matched_dirs) > 0) {
        for (d in matched_dirs) {
          param_file <- file.path(rv$cache_root, d, "analysis_params.txt")
          if (file.exists(param_file)) {
            lines <- readLines(param_file)
            time_line <- grep("Time:", lines, value = TRUE)[1]
            pval_line <- grep("P-val:", lines, value = TRUE)[1]
            fc_line <- grep("LogFC:", lines, value = TRUE)[1]
            label <- paste0(d, " [", time_line, " | ", pval_line, " | ", fc_line, "]")
            choices_list[[label]] <- d
          }
        }
      }

      if (length(choices_list) > 0) {
        output$ui_analysis_history <- renderUI({
          selectInput(ns("sel_existing_result"), "Found History:",
                      choices = choices_list, width = "100%")
        })
        output$ui_btn_load_analysis <- renderUI({
          actionButton(ns("btn_load_analysis"), "Load Selected",
                       icon = icon("folder-open"), class = "btn-info", width = "100%")
        })
      } else {
        output$ui_analysis_history <- renderUI({ helpText("No previous results found.") })
        output$ui_btn_load_analysis <- renderUI({ NULL })
      }
    })

    observeEvent(input$btn_load_analysis, {
      req(input$sel_existing_result)
      target_dir <- file.path(rv$cache_root, input$sel_existing_result)
      withProgress(message = "Loading Result...", {
        rv$diff_tool <- readr::read_rds(file.path(target_dir, "diff_tool.rds"))
        if (file.exists(file.path(target_dir, "enrich_tool.rds")))
          rv$enrich_tool <- readr::read_rds(file.path(target_dir, "enrich_tool.rds"))
        if (file.exists(file.path(target_dir, "gsva_tool.rds")))
          rv$gsva_tool <- readr::read_rds(file.path(target_dir, "gsva_tool.rds"))
        rv$loaded_params <- paste(readLines(file.path(target_dir, "analysis_params.txt")),
                                  collapse = "<br>")
      })
      sendSweetAlert(session, "Loaded", "Result loaded.", "success")
      updateNavbarPage(session, "main_nav", selected = "Visualization")
    })

    # --- Run new analysis ---
    observeEvent(input$btn_run_pipeline, {
      req(rv$manager$imputed_se)
      disable("btn_run_pipeline")
      add_log("Starting Analysis Pipeline...")

      params_list <- list(
        filter_col = rv$filter_col_val,
        filter_val = rv$filter_val_val,
        imp_method = rv$imp_method_val,
        mode = input$analysis_mode,
        condition_col = input$sel_condition_col,
        case = if (input$analysis_mode == "group") input$sel_case else NA,
        control = if (input$analysis_mode == "group") input$sel_control else NA,
        covariates = input$sel_covariates,
        pval_cutoff = input$param_pval,
        logfc_cutoff = if (input$analysis_mode == "group") input$param_logfc else NA,
        corr_cutoff = if (input$analysis_mode == "continuous") input$param_corr else NA,
        r2_cutoff = if (input$analysis_mode == "continuous") input$param_r2 else NA,
        enrich_pval = input$param_enrich_pval,
        use_adj_pval = input$param_use_adj,
        run_gsva = input$run_gsva
      )
      short_hash <- substr(digest::digest(params_list, algo = "md5"), 1, 6)
      final_sub_folder <- paste0(get_base_name(), "_", short_hash)
      full_res_dir <- file.path(rv$cache_root, final_sub_folder)

      tryCatch({
        if (!dir.exists(full_res_dir)) dir.create(full_res_dir, recursive = TRUE)
        if (!dir.exists(full_res_dir)) stop("Could not create output directory: ", full_res_dir)
      }, error = function(e) {
        add_log(paste("Error:", e$message))
        sendSweetAlert(session, "Directory Error",
                       paste("Failed to create results directory. Check for special characters in file path.\n", e$message),
                       "error")
        enable("btn_run_pipeline")
        return()
      })

      filter_info <- if (is.null(rv$filter_col_val) || rv$filter_col_val == "none") {
        "All samples"
      } else {
        paste0(rv$filter_col_val, " = ", rv$filter_val_val)
      }

      param_text <- c(
        paste("Analysis Time:", Sys.time()),
        paste("Hash ID:", short_hash),
        paste("Data Filter:", filter_info),
        paste("Imputation Method:", if (!is.null(rv$imp_method_val)) rv$imp_method_val else "N/A"),
        paste("Mode:", input$analysis_mode),
        paste("Condition Column:", input$sel_condition_col),
        paste("Case Group:", if (input$analysis_mode == "group") input$sel_case else "N/A"),
        paste("Control Group:", if (input$analysis_mode == "group") input$sel_control else "N/A"),
        paste("Covariates:", if (length(input$sel_covariates) > 0) paste(input$sel_covariates, collapse = ", ") else "None"),
        paste("Paired Column:", if (!is.null(input$sel_paired_col) && input$sel_paired_col != "None") input$sel_paired_col else "None"),
        paste("P-value Cutoff:", input$param_pval),
        paste("LogFC Cutoff:", if (input$analysis_mode == "group") input$param_logfc else "N/A"),
        paste("Correlation Cutoff:", if (input$analysis_mode == "continuous") input$param_corr else "N/A"),
        paste("R-squared Cutoff:", if (input$analysis_mode == "continuous") input$param_r2 else "N/A"),
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
            rv$diff_tool <- dt
            readr::write_rds(dt, file.path(full_res_dir, "diff_tool.rds"))

            dt$get_sig_proteins(pval_cutoff = input$param_pval,
                                use_adjusted_pval = input$param_use_adj)

            if (!is.null(dt$sig_results) && nrow(dt$sig_results) > 0) {
              et <- EnrichmentAnalyst$new(rv$cache_manager)
              sig_genes <- dt$sig_results$Protein
              universe_genes <- rownames(rv$manager$imputed_se)
              et$ora_up <- et$run_comprehensive_ora(sig_genes, universe_genes,
                                                     input$param_enrich_pval)
              rv$enrich_tool <- et
              readr::write_rds(et, file.path(full_res_dir, "enrich_tool.rds"))
            }

            if (input$run_gsva) {
              gt <- ProteomicsGSVA$new(rv$manager, cache_manager = rv$cache_manager)
              gt$run_gsva(dbs = c("GOBP", "GOMF", "GOCC", "KEGG", "Wiki"))
              rv$gsva_tool <- gt
              readr::write_rds(gt, file.path(full_res_dir, "gsva_tool.rds"))
            }
          } else {
            paired_col_val <- if (!is.null(input$sel_paired_col) && input$sel_paired_col != "None") input$sel_paired_col else NULL
            res_objs <- run_proteomics_pipeline(
              data_manager = rv$manager,
              analysis_type = input$analysis_mode,
              condition_col = input$sel_condition_col,
              control_group = input$sel_control,
              case_group = input$sel_case,
              continuous_col = if (input$analysis_mode == "continuous") input$sel_condition_col else NULL,
              covariates = if (length(input$sel_covariates) == 0) NULL else input$sel_covariates,
              paired_col = paired_col_val,
              results_dir = rv$cache_root,
              sub_folder_name = final_sub_folder,
              pval_cutoff = input$param_pval,
              enrich_pval_cutoff = input$param_enrich_pval,
              logfc_cutoff = input$param_logfc,
              corr_cutoff = input$param_corr,
              r2_cutoff = input$param_r2,
              top_n_labels = 10,
              use_adj_pval_sig = input$param_use_adj,
              run_gsva = input$run_gsva,
              gsva_dbs = c("GOBP", "GOMF", "GOCC", "KEGG", "Wiki", "Reactome"),
              cache_manager = rv$cache_manager
            )
            rv$diff_tool <- res_objs$diff_tool
            rv$enrich_tool <- res_objs$enrich_tool
            rv$gsva_tool <- res_objs$gsva_tool
          }
        })
        rv$loaded_params <- paste(param_text, collapse = "<br>")
        sendSweetAlert(session, "Success", "Analysis complete!", "success")
        updateNavbarPage(session, "main_nav", selected = "Visualization")
        enable("btn_run_pipeline")
      }, error = function(e) {
        add_log(paste("Error:", e$message))
        sendSweetAlert(session, "Pipeline Error", e$message, "error")
        enable("btn_run_pipeline")
      })
    })

    # --- Session Save/Restore ---
    output$btn_save_session <- downloadHandler(
      filename = function() {
        paste0("ProteomicsApp_Session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
      },
      content = function(file) {
        session_data <- list(
          manager      = rv$manager,
          diff_tool    = rv$diff_tool,
          enrich_tool  = rv$enrich_tool,
          gsva_tool    = rv$gsva_tool,
          loaded_params = rv$loaded_params,
          save_time    = Sys.time(),
          app_version  = tryCatch(as.character(packageVersion("ProteomicsApp")),
                                  error = function(e) "dev")
        )
        saveRDS(session_data, file)
        showNotification("Session saved successfully!", type = "message")
      }
    )

    observeEvent(input$upload_session, {
      req(input$upload_session)
      tryCatch({
        session_data <- readRDS(input$upload_session$datapath)
        required_fields <- c("manager", "diff_tool", "loaded_params")
        if (!all(required_fields %in% names(session_data))) {
          showNotification("Invalid session file: missing required fields.", type = "error")
          return()
        }
        rv$manager       <- session_data$manager
        rv$diff_tool     <- session_data$diff_tool
        rv$enrich_tool   <- session_data$enrich_tool
        rv$gsva_tool     <- session_data$gsva_tool
        rv$loaded_params <- session_data$loaded_params

        save_info <- if (!is.null(session_data$save_time)) {
          format(session_data$save_time, "%Y-%m-%d %H:%M")
        } else "unknown"

        showNotification(paste("Session restored from", save_info),
                         type = "message", duration = 5)
        add_log(paste("Session restored from file saved at", save_info))
        updateNavbarPage(session, "main_nav", selected = "Visualization")
      }, error = function(e) {
        showNotification(paste("Failed to restore session:", e$message), type = "error")
      })
    })

    # Return reactive inputs that other modules need
    return(list(
      analysis_mode = reactive(input$analysis_mode),
      sel_condition_col = reactive(input$sel_condition_col),
      sel_control = reactive(input$sel_control),
      sel_case = reactive(input$sel_case),
      param_pval = reactive(input$param_pval),
      param_logfc = reactive(input$param_logfc),
      param_corr = reactive(input$param_corr),
      param_r2 = reactive(input$param_r2),
      param_enrich_pval = reactive(input$param_enrich_pval),
      param_use_adj = reactive(input$param_use_adj),
      run_gsva = reactive(input$run_gsva)
    ))
  })
}
