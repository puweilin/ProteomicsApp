#' Run the complete proteomics differential and enrichment analysis pipeline (modular execution)
#'
#' @param data_manager Pre-initialized ProteomicsDataManager object (must contain imputed_se)
#' @param analysis_type Analysis type: "group" (default) or "continuous"
#' @param condition_col (Group mode) Group column name
#' @param control_group (Group mode) Control group name
#' @param case_group (Group mode) Case/treatment group name
#' @param continuous_col (Continuous mode) Continuous variable column name (e.g., "age", "time_num")
#' @param continuous_method (Continuous mode) Regression method: "spline" (default) or "linear"
#' @param covariates Covariate vector (e.g., "age"), default NULL
#' @param results_dir Root output directory for results
#' @param sub_folder_name (Optional) Subfolder name. Auto-generated if not specified
#' @param pval_cutoff P-value threshold
#' @param logfc_cutoff (Group mode) LogFC threshold
#' @param corr_cutoff (Continuous mode) Correlation (Spearman Rho) threshold
#' @param r2_cutoff (Continuous mode) R-squared threshold
#' @param use_adj_pval_sig (Logical) Whether to use adjusted p-value for significant protein filtering (default TRUE)
#' @param top_n_labels Number of volcano plot labels
#'
#' @param run_diff Whether to run differential expression analysis (default TRUE)
#' @param run_enrichment Whether to run enrichment analysis (default TRUE)
#' @param run_gsva Whether to run GSVA analysis (default TRUE)
#'
#' @param gsva_dbs Database list for GSVA (default: all 6 databases)
#' @param gsva_min_size GSVA gene set minimum gene count
#' @param gsva_max_size GSVA gene set maximum gene count
#'
#' @return Returns a list containing diff_tool, enrich_tool, and gsva_tool objects
run_proteomics_pipeline <- function(data_manager,
                                    # --- Core mode selection ---
                                    analysis_type = c("group", "continuous"),

                                    # --- Group mode parameters ---
                                    condition_col = "condition",
                                    control_group = NULL,
                                    case_group = NULL,

                                    # --- Continuous mode parameters ---
                                    continuous_col = NULL,
                                    continuous_method = "spline",
                                    corr_cutoff = 0.5,
                                    r2_cutoff = 0.5,

                                    # --- General parameters ---
                                    covariates = NULL,
                                    paired_col = NULL,
                                    results_dir = "Results",
                                    sub_folder_name = NULL,
                                    pval_cutoff = 0.05,
                                    enrich_pval_cutoff = 1,
                                    logfc_cutoff = 0.263,
                                    top_n_labels = 10,
                                    use_adj_pval_sig = TRUE,

                                    # --- Module control switches ---
                                    run_diff = TRUE,       # Whether to run differential analysis
                                    run_enrichment = TRUE, # Whether to run enrichment analysis
                                    run_gsva = TRUE,       # Whether to run GSVA analysis

                                    # --- GSVA parameters ---
                                    gsva_dbs = c("GOBP", "GOMF", "GOCC", "KEGG", "Wiki", "Reactome"),
                                    gsva_cont_method = "spearman",
                                    gsva_min_size = 10,
                                    gsva_max_size = 500,
                                    gsva_adjusted = TRUE,

                                    # --- Cache manager ---
                                    cache_manager = NULL
                                    ) {

  # Parameter validation
  analysis_type <- match.arg(analysis_type)

  # --- 1. Directory setup and initialization ---

  # Auto-generate directory name (needed even when skipping diff, to locate data)
  if (is.null(sub_folder_name)) {
    if (analysis_type == "group") {
      if(is.null(control_group) || is.null(case_group)) stop("Group analysis requires 'control_group' and 'case_group'.")
      cov_str <- if(is.null(covariates)) "" else paste0("_adj_", paste(covariates, collapse="_"))
      sub_folder_name <- paste0(case_group, "_vs_", control_group, cov_str)
    } else {
      if(is.null(continuous_col)) stop("Continuous analysis requires 'continuous_col'.")
      sub_folder_name <- paste0("Regression_", continuous_col, "_", continuous_method)
    }
  }

  # Sanitize sub_folder_name to remove characters that can break directory creation
  sub_folder_name <- gsub("[^a-zA-Z0-9._-]", "_", sub_folder_name)

  subresult_dir <- file.path(results_dir, sub_folder_name)
  enrich_dir <- file.path(subresult_dir, "Enrich_ORA")
  gsea_dir <- file.path(subresult_dir, "Enrich_GSEA")
  gsva_dir <- file.path(subresult_dir, "GSVA_Results")
  gsva_plots_dir <- file.path(gsva_dir, "Plots")

  if (!dir.exists(subresult_dir)) dir.create(subresult_dir, recursive = TRUE)
  if (!dir.exists(subresult_dir)) stop("Failed to create output directory: ", subresult_dir,
                                       ". Check that the path does not contain unsupported characters.")
  message(paste("  Output directory:", subresult_dir))

  # Initialize diff_tool (compute if run_diff, otherwise try loading)
  diff_tool <- DiffExpAnalyst$new(data_manager$imputed_se)
  sig_proteins <- c() # Initialize as empty

  # --- 2. Module A: Differential Expression Analysis ---
  if (run_diff) {
    message(">>> Starting Module: Differential Analysis <<<")
    
    if (analysis_type == "group") {
      message(paste0("  Running Group Analysis: ", case_group, " vs ", control_group))
      diff_tool$run_dep_analysis(
        condition_col = condition_col,
        control_group = control_group,
        case_group = case_group,
        covariates = covariates,
        paired_col = paired_col
      )
    } else {
      message(paste0("  Running Continuous Analysis: ", continuous_col))
      diff_tool$run_continuous_analysis(
        time_col = continuous_col,
        method = continuous_method
      )
    }

    # Plot and save
    message("  Generating Volcano Plot...")
    pdf(file.path(subresult_dir, "volcano.pdf"), width = 10, height = 8)
    p <- diff_tool$plot_volcano(
      logfc_cutoff = logfc_cutoff,
      r2_cutoff = r2_cutoff,
      pval_cutoff = pval_cutoff,
      top_n_labels = top_n_labels,
      use_adjusted_pval = use_adj_pval_sig
    )
    print(p)
    dev.off()

    # Extract significant protein list
    sig_proteins <- diff_tool$get_sig_proteins(
      pval_cutoff = pval_cutoff,
      logfc_cutoff = logfc_cutoff,
      corr_cutoff = corr_cutoff,
      r2_cutoff = r2_cutoff,
      use_adjusted_pval = use_adj_pval_sig
    )
    message(paste0("  Significant proteins found: ", length(sig_proteins)))

    # Save results
    openxlsx::write.xlsx(diff_tool$diff_results, file = file.path(subresult_dir, "Total_Proteins.xlsx"))
    openxlsx::write.xlsx(diff_tool$sig_results, file = file.path(subresult_dir, "Sig_Proteins.xlsx"))
    write_rds(diff_tool, file.path(subresult_dir, "diff_tool.rds"))
    
  } else {
    message("--- Skipping Differential Analysis (run_diff = FALSE) ---")
    
    # If skipping diff but running enrichment, must try loading previous results
    if (run_enrichment) {
      rds_path <- file.path(subresult_dir, "diff_tool.rds")
      if (file.exists(rds_path)) {
        message("  [Auto-Load] Found existing diff_tool.rds, loading for enrichment analysis...")
        diff_tool <- read_rds(rds_path)
        
        # Re-extract significant proteins (ensure parameter consistency)
        sig_proteins <- diff_tool$get_sig_proteins(
          pval_cutoff = pval_cutoff,
          logfc_cutoff = logfc_cutoff,
          corr_cutoff = corr_cutoff,
          r2_cutoff = r2_cutoff,
          use_adjusted_pval = use_adj_pval_sig
        )
        message(paste0("  Loaded significant proteins: ", length(sig_proteins)))
      } else {
        warning("  [Error] Cannot run Enrichment: 'diff_tool.rds' not found and run_diff is FALSE.")
        run_enrichment <- FALSE # Force-disable enrichment analysis
      }
    }
  }

  # --- 3. Module B: Enrichment Analysis ---
  enrich_tool <- NULL
  
  if (run_enrichment) {
    message(">>> Starting Module: Enrichment Analysis <<<")
    
    if (!dir.exists(enrich_dir)) dir.create(enrich_dir, recursive = TRUE)
    if (!dir.exists(gsea_dir)) dir.create(gsea_dir, recursive = TRUE)

    if (length(sig_proteins) > 0) {
      enrich_tool <- EnrichmentAnalyst$new(cache_manager)

      # Run DiffExpAnalyst object analysis
      enrich_results <- enrich_tool$analyze_diff_obj(diff_tool, pval_cutoff = enrich_pval_cutoff)

      # A. ORA (Up/Down) export
      enrich_tool$enrich_to_excel(direction = "UP", output_prefix = "Enrich", target_dir = enrich_dir)
      enrich_tool$enrich_to_excel(direction = "DOWN", output_prefix = "Enrich", target_dir = enrich_dir)

      # B. GSEA export
      enrich_tool$gsea_to_excel(output_prefix = "Enrich_GSEA", target_dir = gsea_dir)

      # C. GSEA Plots (dotplot + enrichment curves for top pathways)
      gsea_plots_dir <- file.path(gsea_dir, "Plots")
      if (!dir.exists(gsea_plots_dir)) dir.create(gsea_plots_dir, recursive = TRUE)

      for (db in names(enrich_tool$gsea_res)) {
        res <- enrich_tool$gsea_res[[db]]
        if (is.null(res) || !inherits(res, "gseaResult") || nrow(res@result) == 0) next

        # Dotplot with activated/suppressed facets
        tryCatch({
          p_dot <- enrichplot::dotplot(res, showCategory = 15, split = ".sign",
                                        label_format = 50, color = "pvalue") +
            ggplot2::facet_grid(.~.sign) +
            ggplot2::ggtitle(paste("GSEA:", db))
          ggplot2::ggsave(plot = p_dot,
                           filename = file.path(gsea_plots_dir, paste0("GSEA_Dotplot_", db, ".pdf")),
                           width = 14, height = 8)
        }, error = function(e) {
          message("  Error plotting GSEA dotplot for ", db, ": ", e$message)
        })

        # Enrichment curves for top 5 pathways by NES
        tryCatch({
          top_pathways <- head(res@result[order(abs(res@result$NES), decreasing = TRUE), "ID"], 5)
          for (pw in top_pathways) {
            safe_name <- gsub("[^a-zA-Z0-9_]", "_", substr(pw, 1, 50))
            p_curve <- enrichplot::gseaplot2(res, geneSetID = pw, title = pw)
            ggplot2::ggsave(plot = p_curve,
                             filename = file.path(gsea_plots_dir,
                                                 paste0("GSEA_Curve_", db, "_", safe_name, ".pdf")),
                             width = 10, height = 6)
          }
        }, error = function(e) {
          message("  Error plotting GSEA curves for ", db, ": ", e$message)
        })
      }
      message("  GSEA plots saved.")

      # Save enrich_tool object
      write_rds(enrich_tool, file.path(subresult_dir, "enrich_tool.rds"))
      message("  Enrichment Analysis Completed.")

    } else {
      warning("  No significant proteins found (or loaded). Skipping Enrichment Analysis.")
    }
  } else {
    message("--- Skipping Enrichment Analysis (run_enrichment = FALSE) ---")
  }

  # --- 4. Module C: GSVA Pathway Analysis ---
  gsva_tool <- NULL

  if (run_gsva) {
    message(">>> Starting Module: GSVA Pathway Analysis <<<")
    
    if (!dir.exists(gsva_dir)) dir.create(gsva_dir, recursive = TRUE)
    if (!dir.exists(gsva_plots_dir)) dir.create(gsva_plots_dir, recursive = TRUE)

    # Initialize GSVA tool
    gsva_tool <- ProteomicsGSVA$new(data_manager, cache_manager = cache_manager)

    # Run GSVA for all specified databases
    message(paste("  Running GSVA for databases:", paste(gsva_dbs, collapse = ", ")))
    gsva_tool$run_gsva(
      dbs = gsva_dbs,
      min_size = gsva_min_size,
      max_size = gsva_max_size,
      save_dir = gsva_dir
    )

    # Run differential pathway analysis based on analysis type
    if (analysis_type == "group") {
      message(paste("  Running differential pathway analysis (Group):", case_group, "vs", control_group))

      gsva_tool$run_pathway_diff_analysis(
        dbs = gsva_dbs,
        group_col = condition_col,
        control_group = control_group,
        case_group = case_group,
        covariates = covariates,
        pval_cutoff = pval_cutoff,
        fc_cutoff = logfc_cutoff,
        adjusted = gsva_adjusted,
        save_dir = gsva_dir
      )
    } else {
      message(paste("  Running pathway correlation analysis (Continuous):", continuous_col))

      gsva_tool$run_pathway_continuous_analysis(
        dbs = gsva_dbs,
        continuous_col = continuous_col,
        covariates = covariates,
        method = gsva_cont_method,
        pval_cutoff = pval_cutoff,
        adjusted = gsva_adjusted,
        save_dir = gsva_dir
      )
    }

    # Generate and save plots for each database
    message("  Generating GSVA visualization plots...")

    for (db in gsva_dbs) {
      db_plots_dir <- file.path(gsva_plots_dir, db)
      if (!dir.exists(db_plots_dir)) dir.create(db_plots_dir, recursive = TRUE)

      if (is.null(gsva_tool$gsva_results[[db]])) {
        next
      }

      # A. Volcano plot
      tryCatch({
        # Smart threshold: don't use LogFC cutoff directly in continuous mode
        gsva_cor_cutoff <- if(analysis_type == "continuous") 0.3 else logfc_cutoff
        
        p_volcano <- gsva_tool$plot_pathway_volcano(
          db = db,
          pval_cutoff = pval_cutoff,
          fc_cutoff = gsva_cor_cutoff, 
          top_n = top_n_labels,
          adjusted = gsva_adjusted
        )
        ggsave(plot = p_volcano, filename = file.path(db_plots_dir, paste0("Volcano_", db, ".pdf")),
               width = 10, height = 8)
      }, error = function(e) {
        message(paste("  Error plotting volcano for", db, ":", e$message))
      })

      # B. Heatmap (Top N differential pathways)
      tryCatch({
        if (!is.null(gsva_tool$diff_pathways[[db]]) && nrow(gsva_tool$diff_pathways[[db]]) > 0) {
          ht_heatmap <- gsva_tool$plot_pathway_heatmap(
            db = db,
            top_n = min(30, nrow(gsva_tool$diff_pathways[[db]])),
            group_col = condition_col,
            show_sample_annot = TRUE
          )
          pdf(file.path(db_plots_dir, paste0("Heatmap_", db, ".pdf")), width = 12, height = 10)
          print(ht_heatmap)
          dev.off()
        }
      }, error = function(e) {
        message(paste("  Error plotting heatmap for", db, ":", e$message))
      })

      # C. GSEA-style waterfall barplot
      tryCatch({
        p_bar <- gsva_tool$plot_pathway_gsea_bar(
          db = db,
          top_n = 20,
          sort_by = ifelse(analysis_type == "group", "logFC", "Correlation")
        )
        ggsave(plot = p_bar, filename = file.path(db_plots_dir, paste0("GSEA_Bar_", db, ".pdf")),
               width = 12, height = 8)
      }, error = function(e) {
        message(paste("  Error plotting GSEA bar for", db, ":", e$message))
      })

      # D. Pathway correlation heatmap
      tryCatch({
        if (!is.null(gsva_tool$diff_pathways[[db]]) && nrow(gsva_tool$diff_pathways[[db]]) >= 2) {
          ht_corr <- gsva_tool$plot_pathway_correlation(
            db = db,
            top_n = min(20, nrow(gsva_tool$diff_pathways[[db]]))
          )
          pdf(file.path(db_plots_dir, paste0("Correlation_", db, ".pdf")), width = 12, height = 10)
          print(ht_corr)
          dev.off()
        }
      }, error = function(e) {
        message(paste("  Error plotting correlation for", db, ":", e$message))
      })
    }

    # Save all differential pathway summary
    gsva_tool$save_all_diff_results(file.path(gsva_dir, "All_DiffPathways_Summary.xlsx"))

    # Save gsva_tool object
    write_rds(gsva_tool, file.path(subresult_dir, "gsva_tool.rds"))
    message("  GSVA Analysis Completed.")
    
  } else {
    message("--- Skipping GSVA Analysis (run_gsva = FALSE) ---")
  }

  message(">>> Pipeline Run Finished <<<")

  return(list(diff_tool = diff_tool, enrich_tool = enrich_tool, gsva_tool = gsva_tool))
}