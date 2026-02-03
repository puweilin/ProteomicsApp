# =============================================================================
# ProteomicsGSVA: 蛋白组学GSVA分析类
# =============================================================================

library(R6)
library(tidyverse)
library(GSVA)
library(limma)
library(openxlsx)
library(ggrepel)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(reshape2)
library(msigdbr)

ProteomicsGSVA <- R6Class("ProteomicsGSVA",
  public = list(
    # --- 核心数据存储 ---
    data_manager = NULL,          # ProteomicsDataManager对象
    gsva_matrices = NULL,         # GSVA打分矩阵列表 (每个DB一个matrix)
    gsva_results = NULL,          # 差异分析结果列表 (每个DB一个data.frame)
    diff_pathways = NULL,         # 显著差异通路列表 (每个DB一个data.frame)
    pathway_genesets = NULL,      # 基因集列表 (每个DB一个list)

    # --- 初始化 ---
    initialize = function(data_manager, organism = "hsa") {
      # 检查 data_manager 是否为 R6 类 (简单的 check)
      if (!"ProteomicsDataManager" %in% class(data_manager)) {
        stop("data_manager must be a ProteomicsDataManager object")
      }
      self$data_manager <- data_manager
      self$gsva_matrices <- list()
      self$gsva_results <- list()
      self$diff_pathways <- list()
      self$pathway_genesets <- list()
      private$organism <- organism
      message("ProteomicsGSVA initialized successfully.")
    },

# --- 功能1: 构建单个数据库的基因集 (KEGG=Online, Others=MSigDB) ---
    prepare_geneset = function(db = "GOBP", min_size = 10, max_size = 500) {
      message(paste("--- Preparing gene sets for:", db, "---"))
      gs_list <- list()

      tryCatch({
        
        # >>> 1. KEGG: 强制使用在线 API (KEGGREST) <<<
        if (db == "KEGG") {
          message("  [Online] Fetching latest KEGG pathways from API (KEGGREST)...")
          
          # 检查必要依赖
          if (!requireNamespace("KEGGREST", quietly = TRUE)) {
             stop("Package 'KEGGREST' is required for online KEGG analysis. Please install it using BiocManager::install('KEGGREST').")
          }
          if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
             stop("Package 'org.Hs.eg.db' is required for ID conversion.")
          }
          
          # 1.1 下载数据 (网络请求)
          # 获取 通路ID -> 基因ID (Entrez)
          link_data <- KEGGREST::keggLink("pathway", "hsa") 
          # 获取 通路ID -> 通路名称
          name_data <- KEGGREST::keggList("pathway", "hsa") 
          
          if(length(link_data) == 0 || length(name_data) == 0) {
            stop("Failed to fetch data from KEGG API. Please check your internet connection.")
          }

          # 1.2 数据清洗
          # 去除前缀 (hsa:123 -> 123, path:hsa001 -> hsa001)
          gene_ids <- gsub("hsa:", "", names(link_data))
          path_ids <- gsub("path:", "", unname(link_data))
          
          # 1.3 ID 转换: Entrez ID -> Gene Symbol
          message("  Mapping Entrez IDs to Gene Symbols...")
          gene_map <- clusterProfiler::bitr(unique(gene_ids), 
                                            fromType = "ENTREZID", 
                                            toType = "SYMBOL", 
                                            OrgDb = "org.Hs.eg.db")
          
          # 1.4 构建列表
          df_map <- data.frame(ENTREZID = gene_ids, PathwayID = path_ids) %>%
            inner_join(gene_map, by = "ENTREZID") %>% # 仅保留能转为 Symbol 的基因
            mutate(PathwayName = name_data[PathwayID]) 
          
          # 清洗通路名称 (去除 " - Homo sapiens (human)")
          df_map$PathwayName <- gsub(" - Homo sapiens \\(human\\)", "", df_map$PathwayName)
          df_map$PathwayName <- trimws(df_map$PathwayName) # 去除可能残留的空格
          
          # 生成列表
          gs_list <- split(df_map$SYMBOL, df_map$PathwayName)
          
        } else {
          # >>> 2. 其他数据库 (GO, Wiki, Reactome): 使用本地 MSigDB <<<
          # 这里不再包含 KEGG 的旧版逻辑
          
          df <- NULL
          if (db == "GOBP") {
            message("  Loading GO BP gene sets (MSigDB)...")
            df <- msigdbr(species = "Homo sapiens", category = "C5") %>% dplyr::filter(gs_subcat == "GO:BP")
          } else if (db == "GOMF") {
            message("  Loading GO MF gene sets (MSigDB)...")
            df <- msigdbr(species = "Homo sapiens", category = "C5") %>% dplyr::filter(gs_subcat == "GO:MF")
          } else if (db == "GOCC") {
            message("  Loading GO CC gene sets (MSigDB)...")
            df <- msigdbr(species = "Homo sapiens", category = "C5") %>% dplyr::filter(gs_subcat == "GO:CC")
          } else if (db == "Wiki") {
            message("  Loading Wiki Pathways (MSigDB)...")
            df <- msigdbr(species = "Homo sapiens", category = "C2") %>% dplyr::filter(gs_subcat == "CP:WIKIPATHWAYS")
          } else if (db == "Reactome") {
            message("  Loading Reactome pathways (MSigDB)...")
            df <- msigdbr(species = "Homo sapiens", category = "C2") %>% dplyr::filter(gs_subcat == "CP:REACTOME")
          }
          
          if (!is.null(df) && nrow(df) > 0) {
            gs_list <- split(df$gene_symbol, df$gs_name)
          }
        }

        # >>> 3. 公共后处理: 去重与大小过滤 <<<
        if (length(gs_list) > 0) {
          gs_list <- lapply(gs_list, unique) # 确保基因不重复
          lens <- lengths(gs_list)
          gs_list <- gs_list[lens >= min_size & lens <= max_size]
          message(paste("  Loaded", length(gs_list), "gene sets for", db))
        } else {
          warning(paste("  Database", db, "returned no data or not recognized."))
        }

      }, error = function(e) {
        message(paste("  Error loading", db, ":", e$message))
        gs_list <- list() 
      })

      return(gs_list)
    },

    # --- 功能2: 运行GSVA分析 (适配 GSVA >= 1.52 新版 API) ---
    run_gsva = function(dbs = c("GOBP", "GOMF", "GOCC", "KEGG", "Wiki", "Reactome"),
                         min_size = 10, max_size = 500,
                         kcdf = "Gaussian", method = "gsva",
                         parallel_threads = 4,
                         save_dir = NULL) {
      message("--- Running GSVA Analysis (Per Database) ---")
      
      # 引入 BiocParallel 以支持并行
      if (!requireNamespace("BiocParallel", quietly = TRUE)) {
        stop("Package 'BiocParallel' is required for the new GSVA version.")
      }

      # 检查数据
      if (is.null(self$data_manager$imputed_se)) {
        stop("Please run perform_imputation() on data_manager first.")
      }

      # 获取表达矩阵
      expr_mat <- assay(self$data_manager$imputed_se)
      message(paste("  Expression matrix dimensions:", nrow(expr_mat), "x", ncol(expr_mat)))

      # 如果指定了保存目录，检查并创建
      if (!is.null(save_dir) && !dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
      }

      # --- 设置并行后端 (适配新版 GSVA 的 BPPARAM) ---
      if (parallel_threads > 1) {
        if (.Platform$OS.type == "windows") {
          bpparam <- BiocParallel::SnowParam(workers = parallel_threads, progressbar = TRUE)
        } else {
          bpparam <- BiocParallel::MulticoreParam(workers = parallel_threads, progressbar = TRUE)
        }
        message(paste("  Parallel processing enabled with", parallel_threads, "threads."))
      } else {
        bpparam <- BiocParallel::SerialParam()
      }

      # 遍历每个数据库，独立运行GSVA
      for (db in dbs) {
        message(paste("\n========== Processing:", db, "=========="))

        # 准备基因集
        if (is.null(self$pathway_genesets[[db]])) {
           gs_list <- self$prepare_geneset(db = db, min_size = min_size, max_size = max_size)
           self$pathway_genesets[[db]] <- gs_list
        } else {
           message("  Using cached gene sets.")
           gs_list <- self$pathway_genesets[[db]]
        }

        if (length(gs_list) == 0) {
          warning(paste("  No gene sets available for", db, "- Skipping"))
          next
        }

        # 运行GSVA (新版 API)
        message(paste("  Running GSVA for", db, "using method:", method, "..."))
        
        tryCatch({
            param_obj <- NULL
            
            # --- 1. 构建参数对象 (GSVA >= 1.52 核心变更) ---
            if (method == "gsva") {
              param_obj <- gsvaParam(
                exprData = expr_mat,
                geneSets = gs_list,
                minSize = min_size,
                maxSize = max_size,
                kcdf = kcdf,
                maxDiff = TRUE 
              )
            } else if (method == "ssgsea") {
              param_obj <- ssgseaParam(
                exprData = expr_mat,
                geneSets = gs_list,
                minSize = min_size,
                maxSize = max_size
              )
            } else if (method == "zscore") {
              param_obj <- zscoreParam(
                exprData = expr_mat,
                geneSets = gs_list,
                minSize = min_size,
                maxSize = max_size
              )
            } else if (method == "plage") {
              param_obj <- plageParam(
                exprData = expr_mat,
                geneSets = gs_list,
                minSize = min_size,
                maxSize = max_size
              )
            } else {
              stop(paste("Method", method, "not supported directly in this wrapper."))
            }

            # --- 2. 调用 gsva() ---
            gsva_result <- gsva(param_obj, verbose = FALSE, BPPARAM = bpparam)
            
            # 保存结果
            self$gsva_matrices[[db]] <- gsva_result
            message(paste("  GSVA matrix dimensions:", nrow(gsva_result), "x", ncol(gsva_result)))
    
            # 保存到Excel文件
            if (!is.null(save_dir)) {
              output_file <- file.path(save_dir, paste0("GSVA_", db, "_scores.xlsx"))
              write.xlsx(as.data.frame(gsva_result) %>% rownames_to_column("Pathway"),
                         file = output_file, overwrite = TRUE)
              message(paste("  Saved to:", output_file))
            }
        }, error = function(e) {
            message(paste("  GSVA execution failed for", db, ":", e$message))
            message("  Debugging: Ensure you have BiocParallel installed and GSVA version >= 1.52")
        })
      }

      message("\n--- All GSVA analyses completed ---")
      return(self$gsva_matrices)
    },

    # --- 功能3: 差异通路分析 (按数据库独立分析) ---
    run_pathway_diff_analysis = function(dbs = NULL, group_col, control_group, case_group,
                                          covariates = NULL,
                                          pval_cutoff = 0.05, fc_cutoff = 0.5,
                                          adjusted = TRUE,
                                          save_dir = NULL) {
      message(paste("--- Pathway Differential Analysis:", case_group, "vs", control_group, "---"))

      if (length(self$gsva_matrices) == 0) {
        stop("Please run run_gsva() first.")
      }

      if (is.null(dbs)) {
        dbs <- names(self$gsva_matrices)
      }
      dbs <- intersect(dbs, names(self$gsva_matrices))

      if (length(dbs) == 0) {
        stop("No valid databases specified or available.")
      }

      if (!is.null(save_dir) && !dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
      }

      all_results <- list()

      for (db in dbs) {
        message(paste("\n========== Differential Analysis:", db, "=========="))

        gsva_mat <- self$gsva_matrices[[db]]
        meta <- self$data_manager$meta_data

        if (!group_col %in% colnames(meta)) {
          warning(paste("Group column", group_col, "not found - skipping", db))
          next
        }

        # 验证协变量
        if (!is.null(covariates)) {
          missing_cov <- setdiff(covariates, colnames(meta))
          if (length(missing_cov) > 0) {
            warning(paste("Covariates not found:", paste(missing_cov, collapse = ", "), "- removing"))
            covariates <- setdiff(covariates, missing_cov)
          }
          if (length(covariates) == 0) covariates <- NULL
        }

        # 样本匹配逻辑 (Use Rownames)
        target_meta <- meta %>%
          filter(!!sym(group_col) %in% c(control_group, case_group))
        target_samples <- rownames(target_meta) 

        matrix_samples <- colnames(gsva_mat)
        valid_samples <- intersect(target_samples, matrix_samples)

        if (length(valid_samples) < 3) {
          warning(paste("Not enough valid samples found for", db, "- Skipping"))
          next
        }
        
        message(paste("  Samples matched:", length(valid_samples), "/", length(target_samples)))

        gsva_subset <- gsva_mat[, valid_samples, drop = FALSE]

        # 构建 Design Matrix
        design_df <- meta[valid_samples, ] %>%
          mutate(Group = factor(!!sym(group_col), levels = c(control_group, case_group)))
        
        if (!is.null(covariates)) {
          for (cov in covariates) {
            if (!is.numeric(design_df[[cov]])) {
              design_df[[cov]] <- as.factor(design_df[[cov]])
            }
          }
        }

        formula_str <- "~ 0 + Group"
        if (!is.null(covariates)) {
          formula_str <- paste0(formula_str, " + ", paste(covariates, collapse = " + "))
        }

        design <- model.matrix(as.formula(formula_str), data = design_df)
        
        colnames(design) <- gsub("^Group", "", colnames(design))
        colnames(design) <- make.names(colnames(design)) 

        fit <- lmFit(gsva_subset, design)

        safe_case <- make.names(case_group)
        safe_ctrl <- make.names(control_group)
        
        contrast_str <- paste0(safe_case, " - ", safe_ctrl)
        message(paste("  Contrast:", contrast_str))

        tryCatch({
            contrast_matrix <- makeContrasts(contrasts = contrast_str, levels = design)
            
            fit2 <- contrasts.fit(fit, contrast_matrix)
            fit2 <- eBayes(fit2)
    
            res_table <- topTable(fit2, number = Inf, sort.by = "none")
            res_table$adj.P.Val <- p.adjust(res_table$P.Value, method = "fdr")
    
            p_col <- if (adjusted) "adj.P.Val" else "P.Value"
    
            res_table <- res_table %>%
              as.data.frame() %>%
              rownames_to_column("Pathway") %>%
              mutate(
                Database = db,
                Used_PValue = !!sym(p_col),
                Diff_Status = case_when(
                  !!sym(p_col) < pval_cutoff & logFC > fc_cutoff ~ "Up-regulated",
                  !!sym(p_col) < pval_cutoff & logFC < -fc_cutoff ~ "Down-regulated",
                  TRUE ~ "Not Significant"
                )
              ) %>%
              arrange(!!sym(p_col))
    
            self$gsva_results[[db]] <- res_table
            all_results[[db]] <- res_table
    
            sig_pathways <- res_table %>% filter(!!sym(p_col) < pval_cutoff)
            self$diff_pathways[[db]] <- sig_pathways
    
            message(paste("  Significant pathways:", nrow(sig_pathways)))
    
            if (!is.null(save_dir)) {
              output_file <- file.path(save_dir, paste0("DiffPathways_", db, ".xlsx"))
              self$save_diff_db_to_excel(res_table, sig_pathways, output_file, db)
            }
        }, error = function(e) {
            message(paste("  Error during Limma contrast:", e$message))
        })
      }

      message("\n--- All differential analyses completed ---")
      return(all_results)
    },

    # 3b: 连续变量分析
    run_pathway_continuous_analysis = function(dbs = NULL, continuous_col,
                                                method = "spearman",
                                                covariates = NULL,
                                                pval_cutoff = 0.05,
                                                adjusted = TRUE,
                                                save_dir = NULL) {
      message(paste("--- Pathway Continuous Analysis on:", continuous_col, "---"))

      if (length(self$gsva_matrices) == 0) {
        stop("Please run run_gsva() first.")
      }

      if (is.null(dbs)) {
        dbs <- names(self$gsva_matrices)
      }
      dbs <- intersect(dbs, names(self$gsva_matrices))

      if (!is.null(save_dir) && !dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
      }

      all_results <- list()

      for (db in dbs) {
        message(paste("\n========== Continuous Analysis:", db, "=========="))

        gsva_mat <- self$gsva_matrices[[db]]
        meta <- self$data_manager$meta_data

        if (!continuous_col %in% colnames(meta)) {
          warning(paste("Column", continuous_col, "not found - skipping", db))
          next
        }

        if (!is.null(covariates)) {
          missing_cov <- setdiff(covariates, colnames(meta))
          if (length(missing_cov) > 0) {
            warning(paste("Covariates not found:", paste(missing_cov, collapse = ", "), "- removing"))
            covariates <- setdiff(covariates, missing_cov)
          }
          if (length(covariates) == 0) covariates <- NULL
        }

        samples <- colnames(gsva_mat)
        valid_samples <- intersect(samples, rownames(meta))
        
        if(length(valid_samples) == 0) {
           warning(paste("No matching samples found for", db))
           next
        }

        gsva_subset <- gsva_mat[, valid_samples, drop=FALSE]
        meta_subset <- meta[valid_samples, ] 

        cont_vals <- meta_subset[[continuous_col]]

        if (!is.numeric(cont_vals)) {
          warning(paste("Column", continuous_col, "is not numeric."))
          next
        }

        use_ppcor <- FALSE
        if (!is.null(covariates)) {
          if (requireNamespace("ppcor", quietly = TRUE)) {
            use_ppcor <- TRUE
          } else {
            warning("Package 'ppcor' not installed. Running simple correlation instead.")
          }
        }
        
        results <- list()
        for (pathway in rownames(gsva_mat)) {
          pathway_expr <- gsva_mat[pathway, ]

          if (is.null(covariates) || !use_ppcor) {
            cor_res <- suppressWarnings(cor.test(pathway_expr, cont_vals, method = method))
            results[[pathway]] <- data.frame(
              Pathway = pathway,
              Correlation = cor_res$estimate,
              P_Value = cor_res$p.value,
              Database = db
            )
          } else {
            partial_data <- data.frame(
              Pathway_Expr = pathway_expr,
              Continuous = cont_vals
            )
            for (cov in covariates) {
              partial_data[[cov]] <- meta_subset[[cov]]
            }
            
            partial_res <- ppcor::pcor.test(partial_data$Pathway_Expr, partial_data$Continuous,
                                      partial_data[, covariates, drop = FALSE],
                                      method = method)
            results[[pathway]] <- data.frame(
              Pathway = pathway,
              Correlation = partial_res$estimate,
              P_Value = partial_res$p.value,
              Database = db
            )
          } 
        }

        res_df <- bind_rows(results) %>%
          mutate(adj_P_Value = p.adjust(P_Value, method = "fdr")) %>%
          mutate(
            Used_PValue = if(adjusted) adj_P_Value else P_Value,
            Diff_Status = case_when(
              Used_PValue < pval_cutoff & Correlation > 0 ~ "Up-regulated",
              Used_PValue < pval_cutoff & Correlation < 0 ~ "Down-regulated",
              TRUE ~ "Not Significant"
          )
          ) %>% arrange(P_Value)

        self$gsva_results[[db]] <- res_df
        all_results[[db]] <- res_df

        p_col <- if (adjusted) "adj_P_Value" else "P_Value"
        sig_pathways <- res_df %>% filter(!!sym(p_col) < pval_cutoff)
        self$diff_pathways[[db]] <- sig_pathways

        message(paste("  Significant pathways:", nrow(sig_pathways)))

        if (!is.null(save_dir)) {
          output_file <- file.path(save_dir, paste0("ContinuousPathways_", db, ".xlsx"))
          write.xlsx(list(All_Results = res_df, Significant = sig_pathways),
                     file = output_file, overwrite = TRUE)
          message(paste("  Saved to:", output_file))
        }
      }

      return(all_results)
    },

    # --- 辅助功能: 保存单个DB的差异结果到Excel [美化版] ---
    save_diff_db_to_excel = function(res_table, sig_pathways, output_file, db_name) {
      wb <- createWorkbook()
      
      # --- [新增] 名称美化函数: 去除前缀 + Title Case ---
      clean_df <- function(d) {
        if("Pathway" %in% colnames(d)) {
          d$Pathway_ID <- d$Pathway # 保留原始ID
          # 1. Remove Prefix (up to first underscore)
          clean <- gsub("^[^_]+_", "", d$Pathway)
          # 2. Underscore to Space
          clean <- gsub("_", " ", clean)
          # 3. Title Case
          d$Pathway <- stringr::str_to_title(clean)
          
          # 把 Pathway 挪到第一列
          d <- d %>% dplyr::select(Pathway, Pathway_ID, everything())
        }
        return(d)
      }
      # ------------------------
      
      addWorksheet(wb, paste0(db_name, "_Summary"))
      writeData(wb, paste0(db_name, "_Summary"), clean_df(res_table))

      sig_up <- sig_pathways %>% filter(Diff_Status == "Up-regulated")
      sig_down <- sig_pathways %>% filter(Diff_Status == "Down-regulated")

      if (nrow(sig_up) > 0) {
        addWorksheet(wb, paste0(db_name, "_Up"))
        writeData(wb, paste0(db_name, "_Up"), clean_df(sig_up))
      }

      if (nrow(sig_down) > 0) {
        addWorksheet(wb, paste0(db_name, "_Down"))
        writeData(wb, paste0(db_name, "_Down"), clean_df(sig_down))
      }

      saveWorkbook(wb, output_file, overwrite = TRUE)
      message(paste("  Saved:", output_file))
    },

    # --- 功能4: 绘制通路火山图 (指定数据库) [美化版] ---
    plot_pathway_volcano = function(db = "GOBP", pval_cutoff = 0.05, fc_cutoff = 0,
                                     top_n = 15, adjusted = TRUE) {
      message(paste("--- Plotting Pathway Volcano Plot:", db, "---"))

      if (is.null(self$gsva_results[[db]])) {
        stop(paste("No differential results for", db))
      }

      df <- self$gsva_results[[db]]
      
      # --- [新增] 名称美化: 去除前缀 + Title Case ---
      clean_temp <- gsub("^[^_]+_", "", df$Pathway) # Remove Prefix
      clean_temp <- gsub("_", " ", clean_temp)      # Underscore to Space
      df$Clean_Pathway <- stringr::str_to_title(clean_temp) # Title Case
      # ---------------------

      # 准备绘图数据
      if ("adj.P.Val" %in% colnames(df)) {
        p_col <- if (adjusted) "adj.P.Val" else "P.Value"
        p_lab <- if (adjusted) "Adjusted P-value" else "P-value"
        df <- df %>%
          mutate(
            log_pval = -log10(!!sym(p_col)),
            expression = case_when(
              !!sym(p_col) < pval_cutoff & logFC > fc_cutoff ~ "Up-regulated",
              !!sym(p_col) < pval_cutoff & logFC < -fc_cutoff ~ "Down-regulated",
              TRUE ~ "Not Significant"
            )
          )
        x_lab <- expression(log[2] ~ "Fold Change")
        y_lab <- as.expression(bquote(-log[10] ~ .(p_lab)))
        x_col <- "logFC"
      } else {
        p_col <- if (adjusted) "adj_P_Value" else "P_Value"
        p_lab <- if (adjusted) "Adjusted P-value" else "P-value"
        df <- df %>%
          mutate(
            log_pval = -log10(!!sym(p_col)),
            expression = case_when(
              !!sym(p_col) < pval_cutoff & Correlation > fc_cutoff ~ "Positive",
              !!sym(p_col) < pval_cutoff & Correlation < -fc_cutoff ~ "Negative",
              TRUE ~ "Not Significant"
            )
          )
        x_lab <- "Spearman Correlation"
        y_lab <- as.expression(bquote(-log[10] ~ .(p_lab)))
        x_col <- "Correlation"
      }

      # 选择top N标签
      top_up <- df %>% filter(expression %in% c("Up-regulated", "Positive")) %>% arrange(desc(log_pval)) %>% slice_head(n = top_n)
      top_down <- df %>% filter(expression %in% c("Down-regulated", "Negative")) %>% arrange(desc(log_pval)) %>% slice_head(n = top_n)
      top_labels <- bind_rows(top_up, top_down) %>% distinct()

      # 绘图
      p <- ggplot(df, aes(x = !!sym(x_col), y = log_pval)) +
        geom_point(aes(color = expression), alpha = 0.6, size = 2) +
        scale_color_manual(values = c("Up-regulated" = "#B31B21", "Down-regulated" = "#1465AC",
                                       "Not Significant" = "grey80", "Positive" = "#B31B21", "Negative" = "#1465AC")) +
        geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "black", alpha = 0.5) +
        geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "black", alpha = 0.5) +
        
        # [修改] 使用 Clean_Pathway 作为标签
        geom_text_repel(data = top_labels, aes(label = Clean_Pathway), size = 3,
                        box.padding = 0.3, max.overlaps = Inf) +
                        
        labs(title = paste("Pathway Volcano Plot -", db),
             subtitle = paste(p_lab, "cutoff:", pval_cutoff),
             x = x_lab, y = y_lab) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              legend.position = "bottom",
              plot.title = element_text(face = "bold", hjust = 0.5))

      return(p)
    },

    # --- 功能5: 绘制差异通路热图 (指定数据库) [美化版] ---
    plot_pathway_heatmap = function(db = "GOBP", top_n = 30, pathways = NULL,
                                     group_col = "condition", 
                                     show_sample_annot = TRUE,
                                     cluster_rows = TRUE, cluster_cols = TRUE) {
      message(paste("--- Plotting Pathway Heatmap:", db, "---"))

      if (is.null(self$gsva_matrices[[db]])) {
        stop(paste("No GSVA matrix for", db))
      }

      # 选择通路
      if (!is.null(pathways)) {
        selected_pathways <- intersect(pathways, rownames(self$gsva_matrices[[db]]))
        if (length(selected_pathways) == 0) {
          stop("None of the specified pathways found.")
        }
      } else if (!is.null(self$diff_pathways[[db]])) {
        p_col_name <- if("adj.P.Val" %in% colnames(self$diff_pathways[[db]])) "adj.P.Val" else "P.Value"
        if(!p_col_name %in% colnames(self$diff_pathways[[db]])) p_col_name <- "P_Value"

        selected_pathways <- self$diff_pathways[[db]] %>%
              arrange(!!sym(p_col_name)) %>%
              slice_head(n = top_n) %>%
              pull(Pathway)
      } else {
        selected_pathways <- rownames(self$gsva_matrices[[db]])[1:top_n]
      }

      # 提取子矩阵
      mat <- self$gsva_matrices[[db]][selected_pathways, , drop = FALSE]

      # --- [新增] 名称美化步骤 ---
      # 1. Remove Prefix
      clean_names <- gsub("^[^_]+_", "", rownames(mat))
      # 2. Underscore to Space
      clean_names <- gsub("_", " ", clean_names)
      # 3. Title Case
      clean_names <- stringr::str_to_title(clean_names)
      
      # 替换矩阵行名
      rownames(mat) <- clean_names
      # ---------------------------

      mat_scaled <- t(scale(t(mat)))

      # 样本注释
      ha <- NULL
      if (show_sample_annot) {
        meta <- self$data_manager$meta_data
        
        if (group_col %in% colnames(meta)) {
          sample_names <- colnames(self$gsva_matrices[[db]]) 
          
          # 优先匹配 rownames
          if(all(sample_names %in% rownames(meta))) {
             annot_vals <- meta[sample_names, group_col]
          } else {
             idx <- match(sample_names, meta$label)
             annot_vals <- meta[[group_col]][idx]
          }

          if(any(is.na(annot_vals))) {
             annot_vals <- as.character(annot_vals)
             annot_vals[is.na(annot_vals)] <- "Unknown"
          }

          annot_df <- data.frame(Group = annot_vals)
          colnames(annot_df) <- group_col

          # 对 groups 排序以确保颜色分配一致
          groups <- sort(unique(as.character(annot_vals)))
          if(length(groups) <= 10) {
             colors <- ggsci::pal_npg()(length(groups))
          } else {
             colors <- ggsci::pal_d3("category20")(length(groups))
          }
          group_colors <- setNames(colors, groups)
          
          col_list <- list()
          col_list[[group_col]] <- group_colors

          ha <- HeatmapAnnotation(
            df = annot_df,
            col = col_list,
            show_annotation_name = TRUE
          )
        }
      }

      # 绘制热图
      ht <- Heatmap(
        mat_scaled,
        name = "Z-score",
        top_annotation = ha,
        cluster_rows = cluster_rows,
        cluster_columns = cluster_cols,
        show_row_names = length(selected_pathways) <= 50,
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize = ifelse(length(selected_pathways) > 30, 8, 10)),
        col = colorRamp2(c(-2, 0, 2), c("#3C5488B2", "white", "#E64B35B2")),
        heatmap_legend_param = list(title = "Z-score"),
        column_title = paste(db, "- Top", length(selected_pathways), "Pathways")
      )

      return(ht)
    },

    # --- 功能6: GSEA样式瀑布图 (Barplot, 指定数据库) [美化版] ---
    plot_pathway_gsea_bar = function(db = "GOBP", top_n = 20, sort_by = "logFC") {
      message(paste("--- Plotting GSEA-style Bar Plot:", db, "---"))

      if (is.null(self$gsva_results[[db]])) {
        stop(paste("No differential results for", db))
      }
      
      df <- self$gsva_results[[db]]
      
      if(!sort_by %in% colnames(df)) {
         if("Correlation" %in% colnames(df)) {
             sort_by <- "Correlation"
         } else {
             stop(paste("Sort column", sort_by, "not found in results."))
         }
      }
      
      # --- [新增] 名称美化: 去除前缀 + Title Case ---
      clean_temp <- gsub("^[^_]+_", "", df$Pathway) # Remove Prefix
      clean_temp <- gsub("_", " ", clean_temp)      # Underscore to Space
      df$Clean_Pathway <- stringr::str_to_title(clean_temp)
      # ---------------------

      df <- df %>%
        arrange(desc(abs(!!sym(sort_by)))) %>%
        slice_head(n = top_n) %>%
        mutate(Direction = ifelse(!!sym(sort_by) > 0, "Up/Pos", "Down/Neg"))

      # 绘图 [修改] x轴使用 Clean_Pathway
      p <- ggplot(df, aes(x = reorder(Clean_Pathway, !!sym(sort_by)), y = !!sym(sort_by), fill = Direction)) +
        geom_bar(stat = "identity", width = 0.7) +
        coord_flip() +
        scale_fill_manual(values = c("Up/Pos" = "#B31B21", "Down/Neg" = "#1465AC")) +
        labs(title = paste(db, "- Top", top_n, "Pathways"),
             x = "Pathway", y = sort_by) +
        theme_bw() +
        theme(
          panel.grid = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 9)
        )

      return(p)
    },

    # --- 功能7-11: 保持原样 ---
    # ... (此处省略未改动的辅助绘图函数 plot_pathway_group_comparison 等，如需使用请确保之前的类中包含)
    # 为保证完整性，以下补全剩余未改动函数

    plot_pathway_group_comparison = function(db, pathway_name, group_col,
                                              colors = NULL, add_stats = TRUE) {
      # ... (代码与之前相同，略)
      message(paste("Plotting group comparison for:", pathway_name, "(", db, ")"))

      if (is.null(self$gsva_matrices[[db]])) {
        stop(paste("No GSVA matrix for", db))
      }

      matched_pathway <- rownames(self$gsva_matrices[[db]])[grep(pathway_name, rownames(self$gsva_matrices[[db]]), ignore.case = TRUE)]
      if (length(matched_pathway) == 0) {
        stop(paste("Pathway", pathway_name, "not found in", db))
      }
      pathway_id <- matched_pathway[1]

      expr_vals <- self$gsva_matrices[[db]][pathway_id, ]
      meta <- self$data_manager$meta_data
      
      # Fix matching
      sample_names <- names(expr_vals)
      if(all(sample_names %in% rownames(meta))) {
          meta_subset <- meta[sample_names, ]
      } else {
          meta_subset <- meta[match(sample_names, meta$label), ]
      }

      plot_df <- data.frame(
        Sample = sample_names,
        Pathway_Score = expr_vals,
        Group = meta_subset[[group_col]]
      )

      p <- ggplot(plot_df, aes(x = Group, y = Pathway_Score, fill = Group)) +
        geom_violin(alpha = 0.5, trim = FALSE) +
        geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
        labs(title = paste(db, ":", pathway_name),
             subtitle = pathway_id,
             y = "GSVA Score", x = group_col) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              plot.title = element_text(face = "bold", hjust = 0.5))

      if (!is.null(colors)) {
        p <- p + scale_fill_manual(values = colors)
      } else {
        p <- p + scale_fill_npg()
      }

      if (add_stats && length(unique(plot_df$Group)) >= 2) {
        p <- p + stat_compare_method(method = ifelse(length(unique(plot_df$Group)) == 2, "t.test", "anova"),
                                      label = "p.format")
      }
      return(p)
    },

    plot_pathway_trend = function(db, pathway_name, continuous_col,
                                   add_ci = TRUE, color = "#1465AC") {
      # ... (代码与之前相同)
      message(paste("Plotting trend for:", pathway_name, "(", db, ")"))

      if (is.null(self$gsva_matrices[[db]])) {
        stop(paste("No GSVA matrix for", db))
      }

      matched_pathway <- rownames(self$gsva_matrices[[db]])[grep(pathway_name, rownames(self$gsva_matrices[[db]]), ignore.case = TRUE)]
      if (length(matched_pathway) == 0) {
        stop(paste("Pathway", pathway_name, "not found in", db))
      }
      pathway_id <- matched_pathway[1]

      expr_vals <- self$gsva_matrices[[db]][pathway_id, ]
      meta <- self$data_manager$meta_data
      
      # Fix matching
      sample_names <- names(expr_vals)
      if(all(sample_names %in% rownames(meta))) {
          meta_subset <- meta[sample_names, ]
      } else {
          meta_subset <- meta[match(sample_names, meta$label), ]
      }

      plot_df <- data.frame(
        Sample = sample_names,
        Pathway_Score = expr_vals,
        Cont_Value = meta_subset[[continuous_col]]
      ) %>% filter(!is.na(Cont_Value))

      p <- ggplot(plot_df, aes(x = Cont_Value, y = Pathway_Score)) +
        geom_point(alpha = 0.6, size = 2, color = color) +
        geom_smooth(method = "loess", color = "#B31B21", fill = "#B31B21",
                    alpha = 0.2, se = add_ci) +
        labs(title = paste(db, ":", pathway_name),
             x = continuous_col, y = "GSVA Score") +
        theme_bw() +
        theme(panel.grid = element_blank(),
              plot.title = element_text(face = "bold", hjust = 0.5))

      cor_res <- cor.test(plot_df$Pathway_Score, plot_df$Cont_Value, method = "spearman")
      p <- p + annotate("text", x = Inf, y = Inf,
                        label = paste0("rho = ", round(cor_res$estimate, 3), "\nP = ", formatC(cor_res$p.value, format = "e", digits = 2)),
                        hjust = 1.1, vjust = 1.5, size = 4)
      return(p)
    },

    plot_pathway_correlation = function(db, pathways = NULL, top_n = 20,
                                         method = "spearman", cluster = TRUE) {
      message(paste("--- Plotting Pathway Correlation Heatmap:", db, "---"))

      if (is.null(self$gsva_matrices[[db]])) {
        stop(paste("No GSVA matrix for", db))
      }

      if (!is.null(pathways)) {
        selected_pathways <- intersect(pathways, rownames(self$gsva_matrices[[db]]))
      } else if (!is.null(self$diff_pathways[[db]])) {
        p_col <- if("adj.P.Val" %in% colnames(self$diff_pathways[[db]])) "adj.P.Val" else "adj_P_Value"
        selected_pathways <- self$diff_pathways[[db]] %>%
          arrange(!!sym(p_col)) %>%
          slice_head(n = top_n) %>%
          pull(Pathway)
      } else {
        selected_pathways <- rownames(self$gsva_matrices[[db]])[1:top_n]
      }
      
      if(length(selected_pathways) < 2) stop("Need at least 2 pathways for correlation.")

      mat <- self$gsva_matrices[[db]][selected_pathways, , drop = FALSE]
      cor_mat <- cor(t(mat), method = method)

      ht <- Heatmap(
        cor_mat,
        name = paste0(toupper(method), " Corr"),
        col = colorRamp2(c(-1, 0, 1), c("#3C5488B2", "white", "#E64B35B2")),
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        cluster_rows = cluster,
        cluster_columns = cluster,
        heatmap_legend_param = list(title = "Correlation"),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(!is.na(cor_mat[i, j])) {
             grid.text(sprintf("%.2f", cor_mat[i, j]), x, y, gp = gpar(fontsize = 6))
          }
        }
      )
      return(ht)
    },

    save_all_gsva_matrices = function(save_dir) {
      if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
      }
      for (db in names(self$gsva_matrices)) {
        if (!is.null(self$gsva_matrices[[db]])) {
          output_file <- file.path(save_dir, paste0("GSVA_", db, "_scores.txt"))
          write.table(as.data.frame(self$gsva_matrices[[db]]) %>% rownames_to_column("Pathway"),
                      file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
          message(paste("Saved:", output_file))
        }
      }
    },

    save_all_diff_results = function(output_file) {
      wb <- createWorkbook()
      # Helper (copied for self-containment)
      clean_df <- function(d) {
        if("Pathway" %in% colnames(d)) {
          d$Pathway_ID <- d$Pathway
          clean <- gsub("^[^_]+_", "", d$Pathway) # Remove Prefix
          clean <- gsub("_", " ", clean)
          d$Pathway <- stringr::str_to_title(clean)
          d <- d %>% dplyr::select(Pathway, Pathway_ID, everything())
        }
        return(d)
      }

      for (db in names(self$gsva_results)) {
        if (!is.null(self$gsva_results[[db]])) {
          res <- self$gsva_results[[db]]
          safe_name <- gsub("[^a-zA-Z0-9]", "_", db)

          addWorksheet(wb, paste0(safe_name, "_Summary"))
          writeData(wb, paste0(safe_name, "_Summary"), clean_df(res))

          sig_up <- res %>% filter(Diff_Status == "Up-regulated")
          sig_down <- res %>% filter(Diff_Status == "Down-regulated")

          if (nrow(sig_up) > 0) {
            addWorksheet(wb, paste0(safe_name, "_Up"))
            writeData(wb, paste0(safe_name, "_Up"), clean_df(sig_up))
          }
          if (nrow(sig_down) > 0) {
            addWorksheet(wb, paste0(safe_name, "_Down"))
            writeData(wb, paste0(safe_name, "_Down"), clean_df(sig_down))
          }
        }
      }
      saveWorkbook(wb, output_file, overwrite = TRUE)
      message(paste("All results saved to:", output_file))
    },

    get_gsva_matrix = function(db) {
      if (is.null(self$gsva_matrices[[db]])) {
        stop(paste("No GSVA matrix for", db))
      }
      return(self$gsva_matrices[[db]])
    },

    get_diff_results = function(db) {
      if (is.null(self$gsva_results[[db]])) {
        stop(paste("No differential results for", db))
      }
      return(self$gsva_results[[db]])
    }
  ),

  private = list(
    organism = "hsa"
  )
)