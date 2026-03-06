# =============================================================================
# Module: Visualization (All visualization tabs)
# =============================================================================

visualizationUI <- function(id) {
  ns <- NS(id)
  tabPanel("Visualization", icon = icon("chart-bar"),
    sidebarLayout(
      sidebarPanel(width = 3,
        h4("Visual Controls"),
        htmlOutput(ns("txt_current_params")),
        hr(),
        uiOutput(ns("vis_slider_fc")),
        fluidRow(
          column(8, sliderInput(ns("vis_slider_pval"), "P-val Cutoff (View):",
                                min = 0, max = 0.2, value = 0.05, step = 0.001)),
          column(4, numericInput(ns("vis_input_pval"), "Value:",
                                 value = 0.05, min = 0, max = 1, step = 0.001))
        ),
        hr(),
        downloadButton(ns("btn_download_report"), "Download Report (PDF)",
                        class = "btn-success", width = "100%")
      ),
      mainPanel(width = 9,
        tabsetPanel(
          tabPanel("Volcano & Table", br(),
            fluidRow(
              column(7, plotOutput(ns("vis_volcano"), brush = ns("volcano_brush"), height = "500px")),
              column(5, h5("Selected Points:"), tableOutput(ns("vis_volcano_info")))
            ),
            hr(),
            DTOutput(ns("vis_dep_table"))
          ),
          tabPanel("Heatmap", br(),
            fluidRow(
              column(4, numericInput(ns("heatmap_n"), "Top N Proteins:", 50, min = 10, max = 2000)),
              column(4, checkboxInput(ns("heat_cluster_row"), "Cluster Rows", TRUE)),
              column(4, checkboxInput(ns("heat_cluster_col"), "Cluster Cols", TRUE))
            ),
            plotOutput(ns("vis_heatmap"), height = "700px")
          ),
          tabPanel("PCA",
            plotOutput(ns("vis_pca"), height = "600px")
          ),
          tabPanel("Deep Dive", br(),
            fluidRow(
              column(6, uiOutput(ns("vis_selector_protein"))),
              column(6, helpText("Select up to 9 proteins to display simultaneously"))
            ),
            plotOutput(ns("vis_single_prot"), height = "700px")
          ),
          tabPanel("Enrichment (ORA)", br(),
            fluidRow(
              column(6, selectInput(ns("vis_ora_db"), "Database:",
                                    choices = c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome", "Wiki"))),
              column(6, radioButtons(ns("vis_ora_dir"), "Direction:",
                                     choices = c("UP", "DOWN"), inline = TRUE))
            ),
            plotOutput(ns("vis_ora_plot"), height = "600px"),
            DTOutput(ns("vis_ora_table"))
          ),
          tabPanel("Enrichment (GSEA)", br(),
            fluidRow(
              column(4, selectInput(ns("vis_gsea_db"), "Database:",
                                    choices = c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome"))),
              column(8, uiOutput(ns("ui_vis_gsea_pathway")))
            ),
            plotOutput(ns("vis_gsea_dotplot"), height = "500px"),
            hr(),
            plotOutput(ns("vis_gsea_curve"), height = "400px")
          ),
          tabPanel("GSVA Analysis", br(),
            fluidRow(
              column(4, selectInput(ns("vis_gsva_db"), "Database:",
                                    choices = c("GOBP", "GOMF", "GOCC", "KEGG", "Wiki", "Reactome"))),
              column(4, selectInput(ns("vis_gsva_plot_type"), "Plot Type:",
                                    choices = c("Volcano", "Heatmap", "Barplot")))
            ),
            plotOutput(ns("vis_gsva_plot"), height = "700px")
          )
        )
      )
    )
  )
}

visualizationServer <- function(id, rv, analysis) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # --- Current params display ---
    output$txt_current_params <- renderUI({
      if (is.null(rv$loaded_params)) return(HTML("<i>No analysis loaded.</i>"))
      wellPanel(
        style = "background-color: #f9fafe; font-size: 11px; padding: 5px;",
        HTML(paste0("<b>Current Analysis:</b><br>", rv$loaded_params))
      )
    })

    # --- Dynamic FC slider ---
    output$vis_slider_fc <- renderUI({
      req(rv$diff_tool)
      if ("logFC" %in% colnames(rv$diff_tool$diff_results)) {
        logfc_val <- tryCatch(analysis$param_logfc(), error = function(e) 0.263)
        fluidRow(
          column(8, sliderInput(ns("vis_fc"), "LogFC Cutoff:",
                                min = 0, max = 2, value = logfc_val, step = 0.01)),
          column(4, numericInput(ns("vis_input_fc"), "Value:",
                                 value = logfc_val, min = 0, max = 10, step = 0.01))
        )
      } else {
        NULL
      }
    })

    # --- Sync LogFC slider and input ---
    observeEvent(input$vis_fc, {
      if (!is.null(input$vis_input_fc) && !isTRUE(all.equal(input$vis_fc, input$vis_input_fc))) {
        updateNumericInput(session, "vis_input_fc", value = input$vis_fc)
      }
    }, ignoreInit = TRUE)

    observeEvent(input$vis_input_fc, {
      if (!is.null(input$vis_fc) && !isTRUE(all.equal(input$vis_input_fc, input$vis_fc))) {
        updateSliderInput(session, "vis_fc", value = input$vis_input_fc)
      }
    }, ignoreInit = TRUE)

    # --- Sync P-val slider and input ---
    observeEvent(input$vis_slider_pval, {
      if (!is.null(input$vis_input_pval) && !isTRUE(all.equal(input$vis_slider_pval, input$vis_input_pval))) {
        updateNumericInput(session, "vis_input_pval", value = input$vis_slider_pval)
      }
    }, ignoreInit = TRUE)

    observeEvent(input$vis_input_pval, {
      if (!is.null(input$vis_slider_pval) && !isTRUE(all.equal(input$vis_input_pval, input$vis_slider_pval))) {
        updateSliderInput(session, "vis_slider_pval", value = input$vis_input_pval)
      }
    }, ignoreInit = TRUE)

    # --- Volcano plot ---
    output$vis_volcano <- renderPlot({
      req(rv$diff_tool)
      mode <- analysis$analysis_mode()
      if (mode == "anova") {
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, "N/A for ANOVA", cex = 1.5)
      } else if (mode == "continuous") {
        rv$diff_tool$plot_volcano(
          corr_cutoff = analysis$param_corr(),
          r2_cutoff = analysis$param_r2(),
          pval_cutoff = input$vis_slider_pval,
          use_adjusted_pval = analysis$param_use_adj()
        )
      } else {
        rv$diff_tool$plot_volcano(
          logfc_cutoff = if (is.null(input$vis_fc)) 0.58 else input$vis_fc,
          pval_cutoff = input$vis_slider_pval,
          use_adjusted_pval = analysis$param_use_adj()
        )
      }
    })

    # --- Volcano brush info ---
    output$vis_volcano_info <- renderTable({
      req(rv$diff_tool, input$volcano_brush)
      mode <- analysis$analysis_mode()
      if (mode == "anova") return(NULL)
      df <- rv$diff_tool$diff_results
      is_reg <- "spearman_rho" %in% colnames(df)
      p_col <- if (analysis$param_use_adj()) {
        if (is_reg) "adj_pval" else "adj.P.Val"
      } else {
        if (is_reg) "p_value" else "P.Value"
      }
      plot_df <- df %>%
        dplyr::mutate(x_val = if (is_reg) spearman_rho else logFC,
                      log_pval = -log10(!!rlang::sym(p_col)))
      brushedPoints(plot_df, input$volcano_brush, xvar = "x_val", yvar = "log_pval")
    })

    # --- DEP table ---
    output$vis_dep_table <- renderDT({
      req(rv$diff_tool)
      datatable(rv$diff_tool$diff_results, options = list(scrollX = TRUE))
    })

    # --- Heatmap ---
    output$vis_heatmap <- renderPlot({
      req(rv$diff_tool, rv$manager$imputed_se)
      df <- rv$diff_tool$diff_results
      n_each <- ceiling(input$heatmap_n / 2)
      mode <- analysis$analysis_mode()
      condition_col <- analysis$sel_condition_col()

      if (mode == "group") {
        p_col <- if (analysis$param_use_adj()) "adj.P.Val" else "P.Value"
        fc_cutoff <- if (is.null(input$vis_fc)) 0 else input$vis_fc
        pval_cutoff <- input$vis_slider_pval

        up_prots <- df %>%
          dplyr::filter(logFC > fc_cutoff & !!rlang::sym(p_col) < pval_cutoff) %>%
          dplyr::arrange(!!rlang::sym(p_col)) %>%
          head(n_each) %>%
          dplyr::pull(Protein)
        down_prots <- df %>%
          dplyr::filter(logFC < -fc_cutoff & !!rlang::sym(p_col) < pval_cutoff) %>%
          dplyr::arrange(!!rlang::sym(p_col)) %>%
          head(n_each) %>%
          dplyr::pull(Protein)
        sig_prots <- c(up_prots, down_prots)
        heatmap_title <- paste0("Top ", length(sig_prots), " DEGs (Up: ", length(up_prots),
                                ", Down: ", length(down_prots), ")")

      } else if (mode == "continuous") {
        p_col <- if (analysis$param_use_adj()) "adj_pval" else "p_value"
        corr_cutoff <- analysis$param_corr()
        pval_cutoff <- input$vis_slider_pval

        pos_prots <- df %>%
          dplyr::filter(spearman_rho > corr_cutoff & !!rlang::sym(p_col) < pval_cutoff) %>%
          dplyr::arrange(!!rlang::sym(p_col)) %>%
          head(n_each) %>%
          dplyr::pull(Protein)
        neg_prots <- df %>%
          dplyr::filter(spearman_rho < -corr_cutoff & !!rlang::sym(p_col) < pval_cutoff) %>%
          dplyr::arrange(!!rlang::sym(p_col)) %>%
          head(n_each) %>%
          dplyr::pull(Protein)
        sig_prots <- c(pos_prots, neg_prots)
        heatmap_title <- paste0("Top ", length(sig_prots), " Correlated Proteins (Pos: ",
                                length(pos_prots), ", Neg: ", length(neg_prots), ")")

      } else {
        sig_prots <- rv$diff_tool$get_sig_proteins(
          pval_cutoff = input$vis_slider_pval,
          use_adjusted_pval = analysis$param_use_adj()
        )
        sig_prots <- head(sig_prots, input$heatmap_n)
        heatmap_title <- paste0("ANOVA Significant Proteins (n=", length(sig_prots), ")")
      }

      if (length(sig_prots) == 0) {
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, "No significant proteins found", cex = 1.5)
        return()
      }

      se <- rv$manager$imputed_se
      mat <- SummarizedExperiment::assay(se)[sig_prots, , drop = FALSE]
      mat_scaled <- t(scale(t(mat)))

      anno_df <- as.data.frame(SummarizedExperiment::colData(se))
      ha <- NULL
      if (condition_col %in% colnames(anno_df)) {
        annot_vals <- anno_df[[condition_col]]
        groups <- sort(unique(as.character(annot_vals)))
        if (length(groups) <= 10) {
          colors <- ggsci::pal_npg()(length(groups))
        } else {
          colors <- ggsci::pal_d3("category20")(length(groups))
        }
        group_colors <- setNames(colors, groups)

        annot_df_for_heatmap <- data.frame(annot_vals, row.names = colnames(mat))
        colnames(annot_df_for_heatmap) <- condition_col
        col_list <- list()
        col_list[[condition_col]] <- group_colors

        ha <- ComplexHeatmap::HeatmapAnnotation(
          df = annot_df_for_heatmap, col = col_list, show_annotation_name = TRUE
        )
      }

      ht <- ComplexHeatmap::Heatmap(
        mat_scaled, name = "Z-score", column_title = heatmap_title,
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

    # --- PCA ---
    output$vis_pca <- renderPlot({
      req(rv$manager)
      rv$manager$plot_pca(color_col = analysis$sel_condition_col())
    })

    # --- Deep Dive protein selector ---
    output$vis_selector_protein <- renderUI({
      req(rv$manager)
      pickerInput(ns("vis_sel_prot"), "Select Proteins (max 9):",
                  choices = rownames(rv$manager$imputed_se),
                  multiple = TRUE,
                  options = list(`max-options` = 9, `live-search` = TRUE,
                                 `actions-box` = TRUE,
                                 `selected-text-format` = "count > 3"))
    })

    output$vis_single_prot <- renderPlot({
      req(rv$manager, input$vis_sel_prot)
      rv$manager$plot_protein_expression(input$vis_sel_prot,
                                          variable = analysis$sel_condition_col())
    })

    # --- ORA ---
    output$vis_ora_plot <- renderPlot({
      req(rv$diff_tool)
      mode <- analysis$analysis_mode()
      if (mode == "anova") {
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, "N/A for ANOVA", cex = 1.5)
        return()
      }
      req(rv$enrich_tool)
      tryCatch({
        et <- rv$enrich_tool
        target <- if (input$vis_ora_dir == "UP") et$ora_up[[input$vis_ora_db]]
                  else et$ora_down[[input$vis_ora_db]]
        if (!is.null(target) && inherits(target, "enrichResult") && nrow(target@result) > 0) {
          enrichplot::dotplot(target, showCategory = 15)
        } else {
          plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
          text(1, 1, "No enriched pathways found", cex = 1.5)
        }
      }, error = function(e) {
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, paste("Error:", e$message), cex = 1.2, col = "red")
      })
    })

    output$vis_ora_table <- renderDT({
      req(rv$diff_tool)
      mode <- analysis$analysis_mode()
      if (mode == "anova") return(NULL)
      req(rv$enrich_tool)
      et <- rv$enrich_tool
      target <- if (input$vis_ora_dir == "UP") et$ora_up[[input$vis_ora_db]]
                else et$ora_down[[input$vis_ora_db]]
      if (!is.null(target)) datatable(as.data.frame(target), options = list(scrollX = TRUE))
    })

    # --- GSEA ---
    output$ui_vis_gsea_pathway <- renderUI({
      req(rv$diff_tool)
      mode <- analysis$analysis_mode()
      if (mode == "anova") return(NULL)
      req(rv$enrich_tool)
      res <- rv$enrich_tool$gsea_res[[input$vis_gsea_db]]
      if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
        selectInput(ns("vis_sel_gsea_path"), "Pathway:", choices = res$ID)
      } else {
        NULL
      }
    })

    output$vis_gsea_dotplot <- renderPlot({
      req(rv$diff_tool)
      mode <- analysis$analysis_mode()
      if (mode == "anova") {
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, "N/A for ANOVA", cex = 1.5)
        return()
      }
      req(rv$enrich_tool)
      tryCatch({
        res <- rv$enrich_tool$gsea_res[[input$vis_gsea_db]]
        if (!is.null(res) && inherits(res, "gseaResult") && nrow(res@result) > 0) {
          enrichplot::dotplot(res, split = ".sign") + ggplot2::facet_grid(. ~ .sign)
        } else {
          plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
          text(1, 1, "No GSEA results", cex = 1.5)
        }
      }, error = function(e) {
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, paste("Error:", e$message), cex = 1.2, col = "red")
      })
    })

    output$vis_gsea_curve <- renderPlot({
      req(rv$diff_tool)
      mode <- analysis$analysis_mode()
      if (mode == "anova") {
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, "N/A for ANOVA", cex = 1.5)
        return()
      }
      req(rv$enrich_tool, input$vis_sel_gsea_path)
      res <- rv$enrich_tool$gsea_res[[input$vis_gsea_db]]
      if (!is.null(res)) enrichplot::gseaplot2(res, input$vis_sel_gsea_path)
    })

    # --- GSVA ---
    output$vis_gsva_plot <- renderPlot({
      req(rv$diff_tool)
      mode <- analysis$analysis_mode()
      if (mode == "anova") {
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, "N/A for ANOVA", cex = 1.5)
        return()
      }
      req(rv$gsva_tool)
      gt <- rv$gsva_tool
      db <- input$vis_gsva_db
      condition_col <- analysis$sel_condition_col()
      if (input$vis_gsva_plot_type == "Volcano") {
        gt$plot_pathway_volcano(db, pval_cutoff = input$vis_slider_pval)
      } else if (input$vis_gsva_plot_type == "Heatmap") {
        ht <- gt$plot_pathway_heatmap(db, group_col = condition_col)
        grid::grid.newpage()
        ComplexHeatmap::draw(ht)
      } else {
        gt$plot_pathway_gsea_bar(db)
      }
    })

    # --- PDF Report download ---
    output$btn_download_report <- downloadHandler(
      filename = function() {
        paste0("Proteomics_Report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
      },
      content = function(file) {
        template <- system.file("app", "report_template.Rmd", package = "ProteomicsApp")
        if (template == "") {
          # Development mode
          template <- file.path(dirname(getwd()), "app", "report_template.Rmd")
          if (!file.exists(template)) {
            template <- file.path(getwd(), "report_template.Rmd")
          }
        }
        if (!file.exists(template)) {
          showNotification("Report template not found.", type = "error")
          return()
        }

        temp_rmd <- file.path(tempdir(), "report_template.Rmd")
        file.copy(template, temp_rmd, overwrite = TRUE)

        mode <- tryCatch(analysis$analysis_mode(), error = function(e) "group")

        withProgress(message = "Generating report...", value = 0.5, {
          rmarkdown::render(temp_rmd,
            output_file = file,
            params = list(
              diff_tool       = rv$diff_tool,
              enrich_tool     = rv$enrich_tool,
              gsva_tool       = rv$gsva_tool,
              analysis_params = rv$loaded_params %||% "",
              analysis_mode   = mode
            ),
            envir = new.env(parent = globalenv())
          )
        })
      }
    )
  })
}
