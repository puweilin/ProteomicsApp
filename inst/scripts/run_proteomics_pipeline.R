#' 运行完整的蛋白质组学差异与富集分析流程 (支持模块化运行)
#'
#' @param data_manager 已经初始化的 ProteomicsDataManager 对象 (需包含 imputed_se)
#' @param analysis_type 分析类型: "group" (默认) 或 "continuous"
#' @param condition_col (Group模式) 分组列名
#' @param control_group (Group模式) 对照组名称
#' @param case_group (Group模式) 实验组名称
#' @param continuous_col (Continuous模式) 连续变量列名 (如 "age", "time_num")
#' @param continuous_method (Continuous模式) 回归方法: "spline" (默认) 或 "linear"
#' @param covariates 协变量向量 (例如 "age")，默认为 NULL
#' @param results_dir 结果输出的总根目录
#' @param sub_folder_name (可选) 子文件夹名称。如果不填，会自动生成
#' @param pval_cutoff P值阈值
#' @param logfc_cutoff (Group模式) LogFC 阈值
#' @param corr_cutoff (Continuous模式) Correlation (Spearman Rho) 阈值
#' @param r2_cutoff (Continuous模式) R-squared 阈值
#' @param use_adj_pval_sig (逻辑值) 筛选显著蛋白表时，是否使用校正后 P 值 (默认 TRUE)
#' @param top_n_labels 火山图标注数量
#' 
#' @param run_diff [新增] 是否运行差异表达分析 (默认 TRUE)
#' @param run_enrichment [新增] 是否运行富集分析 (默认 TRUE)
#' @param run_gsva 是否运行GSVA分析 (默认 TRUE)
#' 
#' @param gsva_dbs 运行GSVA的数据库列表 (默认所有6种)
#' @param gsva_min_size GSVA基因集最小基因数
#' @param gsva_max_size GSVA基因集最大基因数
#'
#' @return 返回一个列表，包含 diff_tool, enrich_tool 和 gsva_tool 对象
run_proteomics_pipeline <- function(data_manager,
                                    # --- 核心模式选择 ---
                                    analysis_type = c("group", "continuous"),

                                    # --- Group 模式参数 ---
                                    condition_col = "condition",
                                    control_group = NULL,
                                    case_group = NULL,

                                    # --- Continuous 模式参数 ---
                                    continuous_col = NULL,
                                    continuous_method = "spline",
                                    corr_cutoff = 0.5,
                                    r2_cutoff = 0.5,

                                    # --- 通用参数 ---
                                    covariates = NULL,
                                    paired_col = NULL,
                                    results_dir = "Results",
                                    sub_folder_name = NULL,
                                    pval_cutoff = 0.05,
                                    enrich_pval_cutoff = 1,
                                    logfc_cutoff = 0.263,
                                    top_n_labels = 10,
                                    use_adj_pval_sig = TRUE,

                                    # --- [新增] 模块控制开关 ---
                                    run_diff = TRUE,       # 是否跑差异分析
                                    run_enrichment = TRUE, # 是否跑富集分析
                                    run_gsva = TRUE,       # 是否跑GSVA分析

                                    # --- GSVA 参数 ---
                                    gsva_dbs = c("GOBP", "GOMF", "GOCC", "KEGG", "Wiki", "Reactome"),
                                    gsva_cont_method = "spearman",
                                    gsva_min_size = 10,
                                    gsva_max_size = 500,
                                    gsva_adjusted = TRUE
                                    ) {

  require(here)
  require(openxlsx)
  require(readr)

  # 参数校验
  analysis_type <- match.arg(analysis_type)

  # --- 1. 目录准备与初始化 ---
  
  # 自动生成目录名 (即使不跑 diff，也需要知道目录在哪里以加载数据)
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

  subresult_dir <- file.path(results_dir, sub_folder_name)
  enrich_dir <- file.path(subresult_dir, "Enrich_ORA")
  gsea_dir <- file.path(subresult_dir, "Enrich_GSEA")
  gsva_dir <- file.path(subresult_dir, "GSVA_Results")
  gsva_plots_dir <- file.path(gsva_dir, "Plots")

  if (!dir.exists(subresult_dir)) dir.create(subresult_dir, recursive = TRUE)
  message(paste("  Output directory:", subresult_dir))

  # 初始化 diff_tool (如果跑 diff 则计算，不跑则尝试加载)
  diff_tool <- DiffExpAnalyst$new(data_manager$imputed_se)
  sig_proteins <- c() # 初始化为空

  # --- 2. 模块 A: 差异表达分析 (Differential Analysis) ---
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

    # 绘图与保存
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

    # 提取显著列表
    sig_proteins <- diff_tool$get_sig_proteins(
      pval_cutoff = pval_cutoff,
      logfc_cutoff = logfc_cutoff,
      corr_cutoff = corr_cutoff,
      r2_cutoff = r2_cutoff,
      use_adjusted_pval = use_adj_pval_sig
    )
    message(paste0("  Significant proteins found: ", length(sig_proteins)))

    # 保存结果
    openxlsx::write.xlsx(diff_tool$diff_results, file = file.path(subresult_dir, "Total_Proteins.xlsx"))
    openxlsx::write.xlsx(diff_tool$sig_results, file = file.path(subresult_dir, "Sig_Proteins.xlsx"))
    write_rds(diff_tool, file.path(subresult_dir, "diff_tool.rds"))
    
  } else {
    message("--- Skipping Differential Analysis (run_diff = FALSE) ---")
    
    # 如果不跑 Diff 但要跑 Enrich，必须尝试加载旧结果
    if (run_enrichment) {
      rds_path <- file.path(subresult_dir, "diff_tool.rds")
      if (file.exists(rds_path)) {
        message("  [Auto-Load] Found existing diff_tool.rds, loading for enrichment analysis...")
        diff_tool <- read_rds(rds_path)
        
        # 重新提取显著蛋白 (确保参数一致)
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
        run_enrichment <- FALSE # 强制关闭富集分析
      }
    }
  }

  # --- 3. 模块 B: 富集分析 (Enrichment) ---
  enrich_tool <- NULL
  
  if (run_enrichment) {
    message(">>> Starting Module: Enrichment Analysis <<<")
    
    if (!dir.exists(enrich_dir)) dir.create(enrich_dir, recursive = TRUE)
    if (!dir.exists(gsea_dir)) dir.create(gsea_dir, recursive = TRUE)

    if (length(sig_proteins) > 0) {
      enrich_tool <- EnrichmentAnalyst$new()

      # 运行 DiffExpAnalyst 对象分析
      enrich_results <- enrich_tool$analyze_diff_obj(diff_tool, pval_cutoff = enrich_pval_cutoff)

      # A. ORA (Up/Down) 导出
      enrich_tool$enrich_to_excel(direction = "UP", output_prefix = "Enrich", target_dir = enrich_dir)
      enrich_tool$enrich_to_excel(direction = "DOWN", output_prefix = "Enrich", target_dir = enrich_dir)

      # B. GSEA 导出
      enrich_tool$gsea_to_excel(output_prefix = "Enrich_GSEA", target_dir = gsea_dir)

      # 保存 enrich_tool 对象
      write_rds(enrich_tool, file.path(subresult_dir, "enrich_tool.rds"))
      message("  Enrichment Analysis Completed.")

    } else {
      warning("  No significant proteins found (or loaded). Skipping Enrichment Analysis.")
    }
  } else {
    message("--- Skipping Enrichment Analysis (run_enrichment = FALSE) ---")
  }

  # --- 4. 模块 C: GSVA Pathway 分析 ---
  gsva_tool <- NULL

  if (run_gsva) {
    message(">>> Starting Module: GSVA Pathway Analysis <<<")
    
    if (!dir.exists(gsva_dir)) dir.create(gsva_dir, recursive = TRUE)
    if (!dir.exists(gsva_plots_dir)) dir.create(gsva_plots_dir, recursive = TRUE)

    # 初始化 GSVA 工具
    gsva_tool <- ProteomicsGSVA$new(data_manager)

    # 6.1 运行 GSVA (所有指定数据库)
    message(paste("  Running GSVA for databases:", paste(gsva_dbs, collapse = ", ")))
    gsva_tool$run_gsva(
      dbs = gsva_dbs,
      min_size = gsva_min_size,
      max_size = gsva_max_size,
      save_dir = gsva_dir
    )

    # 6.2 根据分析类型运行差异通路分析
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

    # 6.3 绘制并保存各数据库的图表
    message("  Generating GSVA visualization plots...")

    for (db in gsva_dbs) {
      db_plots_dir <- file.path(gsva_plots_dir, db)
      if (!dir.exists(db_plots_dir)) dir.create(db_plots_dir, recursive = TRUE)

      if (is.null(gsva_tool$gsva_results[[db]])) {
        next
      }

      # A. 火山图
      tryCatch({
        # 智能阈值：Continuous 模式下不要直接用 LogFC cutoff
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

      # B. 热图 (Top N 差异通路)
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

      # C. GSEA 瀑布图
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

      # D. 通路相关性热图
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

    # 6.4 保存所有差异通路汇总
    gsva_tool$save_all_diff_results(file.path(gsva_dir, "All_DiffPathways_Summary.xlsx"))

    # 保存 gsva_tool 对象
    write_rds(gsva_tool, file.path(subresult_dir, "gsva_tool.rds"))
    message("  GSVA Analysis Completed.")
    
  } else {
    message("--- Skipping GSVA Analysis (run_gsva = FALSE) ---")
  }

  message(">>> Pipeline Run Finished <<<")

  return(list(diff_tool = diff_tool, enrich_tool = enrich_tool, gsva_tool = gsva_tool))
}