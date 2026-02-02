library(R6)
library(tidyverse)
library(DEP)
library(SummarizedExperiment)
library(missForest)
library(doParallel)
library(Mfuzz)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(splines)
library(broom)
library(openxlsx)
library(ggrepel)
library(limma)
library(here)
library(patchwork)
library(glmnet)
library(caret)
library(ComplexHeatmap)
library(circlize)

# --- Class 1: 数据管理与预处理 ---
# --- Class 1: 数据管理与预处理 (Updated with Missing Filter) ---
ProteomicsDataManager <- R6Class("ProteomicsDataManager",
  public = list(
    raw_mat = NULL,
    meta_data = NULL,
    meta_data_backup = NULL,
    annot_data = NULL,
    se_obj = NULL,
    imputed_se = NULL,
    imputed_se_backup = NULL,
    valid_protein_numbers = NULL,
    tolerate_missing_percent = NULL,  # [新增] 缺失值容忍度 (0~1)，默认 0.5 表示允许 50% 缺失，超过则过滤
    imputation_method = NULL,  # [新增] 保存使用的 imputation 方法
    missing_mask = NULL,  # [新增] 保存imputation前的缺失位置mask (TRUE = 原始缺失值)
    
    # --- 初始化 ---
    initialize = function(mat_file, meta_file, annot_file, tolerate_missing_percent) {
      self$raw_mat <- mat_file
      self$meta_data <- meta_file
      self$annot_data <- annot_file
      self$tolerate_missing_percent <- tolerate_missing_percent
    },
    
    # --- 核心修正：数据预处理与构建 SE 对象 ---
    # [修改] 支持基于metadata中任意列进行筛选
    process_data = function(filter_col = NULL, filter_value = NULL) {

message("--- Step 1: Data Filtering & Matching ---")

      # 1. 识别矩阵结构
      mat_col_names <- colnames(self$raw_mat)
      id_col_name   <- mat_col_names[1]
      all_sample_cols <- mat_col_names[-1]

      # 2. 筛选 Metadata (支持任意列)
      if (!is.null(filter_col) && !is.null(filter_value)) {
        if (!filter_col %in% colnames(self$meta_data)) {
          stop(paste0("Error: '", filter_col, "' column not found in metadata."))
        }
        target_meta <- self$meta_data %>% dplyr::filter(.data[[filter_col]] == filter_value)
        message(paste0("  Filtering by ", filter_col, " = ", filter_value))
      } else {
        target_meta <- self$meta_data
        message("  No filtering applied, using all samples.")
      }
      
      # 3. 取交集
      valid_samples <- intersect(all_sample_cols, target_meta$label)
      if (length(valid_samples) == 0) stop("Error: No common samples found.")
      
      message(paste("  Selected", length(valid_samples), "samples."))
      
      # 4. 对齐数据
      self$raw_mat <- self$raw_mat %>% dplyr::select(all_of(c(id_col_name, valid_samples)))
      self$meta_data <- target_meta %>% dplyr::filter(label %in% valid_samples)
      
      # 5. ID 转换
      message("--- Step 2: Mapping Protein IDs to Gene Names ---")
      annot_join_col <- colnames(self$annot_data)[1]
      
      # 自动识别 annot_data 中的 Gene Name 列 (假设含有 'Gene' 字样，否则默认第二列)
      gene_name_col <- grep("Gene", colnames(self$annot_data), value = TRUE, ignore.case = TRUE)[1]
      if(is.na(gene_name_col)) gene_name_col <- colnames(self$annot_data)[2]
      
      # 重命名以便统一处理
      annot_clean <- self$annot_data %>% 
        dplyr::select(all_of(c(annot_join_col, gene_name_col))) %>%
        setNames(c(annot_join_col, "Gene.Name.Mapped"))
        
      merged_data <- self$raw_mat %>%
        left_join(annot_clean, by = setNames(annot_join_col, id_col_name)) %>%
        filter(!is.na(Gene.Name.Mapped) & Gene.Name.Mapped != "")
      
      # 处理重复基因名
      unique_genes <- make.unique(as.character(merged_data$Gene.Name.Mapped))
      
      # 6. 构建初始表达矩阵 (用于构建 SE)
      # 这里我们需要构建一个符合 DEP make_se 要求的 data.frame
      # 必须包含 "name" (Gene Name) 和 "ID" (Protein ID) 列
      
      df_for_se <- merged_data %>%
        dplyr::select(all_of(valid_samples)) %>%
        mutate(name = unique_genes, ID = merged_data[[id_col_name]]) %>%
        dplyr::select(name, ID, everything())
      
      # 7. 创建 SummarizedExperiment 对象
      # columns 应该是样本列的索引
      sample_cols_indices <- which(colnames(df_for_se) %in% valid_samples)
      
      self$se_obj <- make_se(df_for_se, columns = sample_cols_indices, expdesign = self$meta_data) %>%
        normalize_vsn()
      
      # 更新 meta_data 为 SE 中的 colData (确保顺序一致)
      self$meta_data <- colData(self$se_obj) %>% data.frame()
      
      # --- [修正] Step 2.5: 直接在 SE 对象层面过滤缺失值 ---
      message(paste0("--- Step 2.5: Filtering proteins with > ", self$tolerate_missing_percent*100, "% missing values ---"))
      
      se_mat <- assay(self$se_obj)
      missing_rate <- rowMeans(is.na(se_mat))
      keep_rows <- missing_rate <= self$tolerate_missing_percent
      
      old_count <- nrow(self$se_obj)
      self$se_obj <- self$se_obj[keep_rows, ]
      
      self$valid_protein_numbers <- nrow(self$se_obj)
      message(paste("  Removed", old_count - self$valid_protein_numbers, "proteins. Remaining:", self$valid_protein_numbers))
      
      message("--- SummarizedExperiment Object Created Successfully ---")
    },
    
    # --- 功能: 绘制缺失值模式 ---
    plot_missing_pattern = function() {
      if (is.null(self$se_obj)) stop("Please run process_data() first.")
      message("--- Plotting Missing Value Pattern using DEP::plot_missval ---")
      p <- plot_missval(self$se_obj)
    },

perform_imputation = function(method = c("missForest","bpca", "knn", "QRILC", "MLE", "MinDet", "MinProb",
      "man", "min", "zero", "mixed", "nbavg"), cores = 4) {
      
      if (is.null(self$se_obj)) stop("Please run process_data() first.")
      method <- match.arg(method)

      message(paste("--- Step 3: Running Imputation with Method:", method, "---"))

      data_matrix <- assay(self$se_obj)

      # [新增] 保存imputation前的缺失位置mask
      self$missing_mask <- is.na(data_matrix)
      message(paste0("  Detected ", sum(self$missing_mask), " missing values (",
                    round(sum(self$missing_mask) / prod(dim(data_matrix)) * 100, 2), "%)"))

      if (method == "missForest") {
        # 依赖包检查
        if (!requireNamespace("missForest", quietly = TRUE)) stop("Package 'missForest' is required.")
        if (!requireNamespace("doParallel", quietly = TRUE)) stop("Package 'doParallel' is required.")
        if (!requireNamespace("parallel", quietly = TRUE)) stop("Package 'parallel' is required.")

        # missForest 需要行是观测(Sample)，列是变量(Protein)，所以先转置
        data_matrix_t <- t(data_matrix)
        
        # --- [修正] 针对不同系统的并行处理 ---
        message(paste0("  Setting up parallel backend for ", .Platform$OS.type, " with ", cores, " cores..."))
        
        if (.Platform$OS.type == "windows") {
          # Windows: 必须使用 Socket Cluster (PSOCK)
          cl <- parallel::makePSOCKcluster(cores)
          doParallel::registerDoParallel(cl)
          
          # 确保函数结束时关闭集群，释放内存
          on.exit(parallel::stopCluster(cl), add = TRUE)
          
        } else {
          # Mac/Linux: 可以直接使用 Fork 机制
          doParallel::registerDoParallel(cores)
        }
        
        message("  Running missForest (this may take time)...")
        set.seed(123)
        # parallelize = 'forests' 通常比 'variables' 更快且更稳健
        mf_result <- missForest::missForest(data_matrix_t, verbose = TRUE, parallelize = "forests")
        
        # 转置回 SE 格式 (行=Protein, 列=Sample)
        imputed_matrix <- t(mf_result$ximp)

      } else {
        # DEP 包的其他填补方法
        # fun 参数可选: "MinProb", "MinDet", "QRILC", "Man", "Min" 等
        imputed_matrix <- assay(DEP::impute(self$se_obj, fun = method))
      }

      self$imputed_se <- self$se_obj
      self$imputed_se_backup <- self$se_obj

      assay(self$imputed_se) <- imputed_matrix
      assay(self$imputed_se_backup) <- imputed_matrix
      self$meta_data_backup  <- self$meta_data
      self$imputation_method <- method  # [新增] 保存使用的方法

      message("--- Imputation Complete ---")
      message("--- Backup created. You can restore full data using reset_data() ---")

      invisible(list(method = method))
    },

# --- [修正版] 功能 2.5: 评估填补效果 (修复绘图列名报错 + 修复删除outlier后维度不匹配) ---
    assess_imputation = function(plot_type = c("structure", "density"),
                                 check_integrity = TRUE,
                                 scale_data = TRUE,
                                 color_col = NULL) {

      if (is.null(self$se_obj)) stop("Please run process_data() first.")
      if (is.null(self$imputed_se)) stop("Please run perform_imputation() first.")

      # 1. 参数处理
      plot_type <- match.arg(plot_type, several.ok = TRUE)
      library(vegan)
      library(patchwork)

      # [新增] 修复维度不匹配问题: 如果删除了outlier，使用backup数据
      raw_mat <- assay(self$se_obj)
      imp_mat <- assay(self$imputed_se)

      # 检查维度是否匹配
      if (ncol(raw_mat) != ncol(imp_mat)) {
        warning(paste0("Sample numbers differ (raw: ", ncol(raw_mat), ", imputed: ", ncol(imp_mat), "). Using backup data for comparison."))

        if (!is.null(self$imputed_se_backup)) {
          # 使用backup数据，但只保留当前存在的样本
          common_samples <- intersect(colnames(self$imputed_se), colnames(self$imputed_se_backup))
          if (length(common_samples) == 0) {
            stop("No common samples found between current and backup data. Cannot assess imputation.")
          }

          # 从backup中提取se_obj对应的原始数据（imputation前的）
          # 注意：se_obj是imputation前的，imputed_se_backup是imputation后的完整数据
          # 我们需要使用imputed_se_backup来反推原始数据
          message(paste0("  Using ", length(common_samples), " common samples for assessment."))
          imp_mat <- assay(self$imputed_se)  # 当前的imputed数据（可能删除了outlier）

          # 从backup中获取相同样本的数据进行比较
          if (!is.null(self$imputed_se_backup) && ncol(self$imputed_se_backup) > ncol(imp_mat)) {
            raw_mat_backup <- assay(self$imputed_se_backup)[, common_samples, drop = FALSE]
            # 注意：这里我们实际上是在比较删除outlier前后的数据，而不是imputation前后
            # 为了正确评估，我们应该跳过这个分析，或者只做density plot
            message("  Note: Structure comparison is not available when outliers have been removed.")
            message("  Showing distribution comparison only.")
            plot_type <- "density"  # 只显示density plot
          }
        } else {
          stop("Dimension mismatch and no backup available. Please re-run imputation.")
        }
      }

      # 2. 准备基础数据
      n_imputed <- sum(is.na(raw_mat))
      prop_imputed <- round(n_imputed / prod(dim(raw_mat)) * 100, 2)
      message(paste0("--- Imputation Assessment (", prop_imputed, "% values imputed) ---"))

      stats_list <- list()
      plot_list <- list()
      
      # ==========================================================
      # 模块 A: 结构稳健性分析 (Structure / MDS-Procrustes)
      # ==========================================================
      if ("structure" %in% plot_type) {
        message("  1. Calculating Sample Distances (Robustness Check)...")
        
        # A.1 数据准备与标准化
        raw_for_dist <- raw_mat
        row_means <- rowMeans(raw_for_dist, na.rm = TRUE)
        idx <- which(is.na(raw_for_dist), arr.ind = TRUE)
        raw_for_dist[idx] <- row_means[idx[,1]]
        
        imp_for_dist <- imp_mat
        
        if (scale_data) {
          raw_for_dist <- t(scale(t(raw_for_dist)))
          imp_for_dist <- t(scale(t(imp_for_dist)))
          raw_for_dist[is.nan(raw_for_dist)] <- 0
          imp_for_dist[is.nan(imp_for_dist)] <- 0
        }
        
        # A.2 MDS 计算
        dist_raw <- dist(t(raw_for_dist))       
        dist_imp <- dist(t(imp_for_dist))            
        mds_raw <- cmdscale(dist_raw, k = 2)    
        mds_imp <- cmdscale(dist_imp, k = 2)    
        
        # A.3 Procrustes 分析
        pro_res <- vegan::protest(mds_raw, mds_imp, scores = "sites", permutations = 999)
        m2_val <- round(pro_res$ss, 4)
        
        message(paste0("  -> Procrustes M2: ", m2_val, " (Lower is better)"))
        stats_list$robustness <- list(m2 = m2_val, pval = pro_res$signif, correlation = pro_res$t0)
        
        # A.4 [关键修复] 强制重命名列，确保 left_join 能正确产生后缀
        df_raw <- as.data.frame(pro_res$X)
        colnames(df_raw)[1:2] <- c("Dim1", "Dim2") # 强制命名
        df_raw <- df_raw %>% mutate(Type = "Raw") %>% rownames_to_column("Sample")
        
        df_imp <- as.data.frame(pro_res$Yrot)
        colnames(df_imp)[1:2] <- c("Dim1", "Dim2") # 强制命名
        df_imp <- df_imp %>% mutate(Type = "Imputed") %>% rownames_to_column("Sample")
        
        # 现在 Dim1 和 Dim2 会发生冲突，left_join 会自动添加后缀
        plot_df <- left_join(df_raw, df_imp, by = "Sample", suffix = c(".raw", ".imp"))
        
        # 处理颜色列
        meta <- as.data.frame(colData(self$imputed_se)) %>% rownames_to_column("Sample")
        plot_df <- left_join(plot_df, meta, by = "Sample")
        
        target_color <- color_col
        if (is.null(target_color)) {
          if ("condition" %in% colnames(meta)) target_color <- "condition"
          else if ("group" %in% colnames(meta)) target_color <- "group"
          else target_color <- colnames(meta)[2]
        }
        
        if (!target_color %in% colnames(meta)) {
          warning(paste("Color column", target_color, "not found. Using default."))
          target_color <- colnames(meta)[2]
        }
        
        # 更新绘图代码使用新的列名 (Dim1.raw, Dim1.imp)
        p_struc <- ggplot(plot_df, aes(x = Dim1.raw, y = Dim2.raw)) +
          geom_segment(aes(xend = Dim1.imp, yend = Dim2.imp), 
                       arrow = arrow(length = unit(0.2, "cm")), 
                       color = "grey50", alpha = 0.5) +
          geom_point(aes(color = .data[[target_color]]), size = 3, alpha = 0.8) +
          geom_point(aes(x = Dim1.imp, y = Dim2.imp, color = .data[[target_color]]), shape = 1, size = 3) +
          labs(
            title = "Structure Drift (MDS)",
            subtitle = paste0("M2 = ", m2_val, " | Corr = ", round(pro_res$t0, 3)),
            x = "Dimension 1", y = "Dimension 2",
            color = target_color
          ) +
          theme_bw() +
          theme(plot.title = element_text(face = "bold", hjust = 0.5))
        
        plot_list$structure <- p_struc
      }
      
      # ==========================================================
      # 模块 B: 分布一致性 (Density Plot)
      # ==========================================================
      if ("density" %in% plot_type) {
        plot_df_dens <- data.frame(
          Intensity = as.vector(imp_mat),
          Type = ifelse(as.vector(is.na(raw_mat)), "Imputed", "Observed")
        )
        
        p_dist <- ggplot(plot_df_dens, aes(x = Intensity, fill = Type, color = Type)) +
          geom_density(alpha = 0.4) +
          scale_fill_manual(values = c("Observed" = "#1465AC", "Imputed" = "#B31B21")) +
          scale_color_manual(values = c("Observed" = "#1465AC", "Imputed" = "#B31B21")) +
          labs(
            title = "Distribution Check",
            subtitle = paste0("Total Imputed: ", prop_imputed, "%"),
            x = "Intensity (Log2)", y = "Density"
          ) +
          theme_bw() + theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "top")
        
        plot_list$density <- p_dist
      }

      # ==========================================================
      # 模块 C: 数据完整性检查
      # ==========================================================
      if (check_integrity) {
        mask <- !is.na(raw_mat)
        diff_vals <- abs(raw_mat[mask] - imp_mat[mask])
        if(max(diff_vals) < 1e-9) {
          message("  [Pass] Data Integrity: Original values preserved.")
        } else {
          warning(paste0("  [Alert] Data Integrity: Original values changed! Max diff: ", max(diff_vals)))
        }
      }

      # ==========================================================
      # 最终展示
      # ==========================================================
      final_plot <- NULL
      if (length(plot_list) == 2) {
        final_plot <- plot_list$structure + plot_list$density
      } else if (length(plot_list) == 1) {
        final_plot <- plot_list[[1]]
      }
      
      if (!is.null(final_plot)) print(final_plot)
      
      return(invisible(list(plots = final_plot, stats = stats_list)))
    },
    
    # --- 功能 3: 绘制 PCA 图 ---
    plot_pca = function(color_col = NULL) {
      if (is.null(self$imputed_se)) stop("Please run perform_imputation() first.")
      
      pca_res <- prcomp(t(assay(self$imputed_se)), scale. = TRUE)
      
      pca_df <- as.data.frame(pca_res$x) %>% 
        rownames_to_column("sample") %>%
        left_join(as.data.frame(colData(self$imputed_se)) %>% rownames_to_column("sample"), by="sample")
      
      color_col <- color_col
      
      p <- ggplot(pca_df, aes(x=PC1, y=PC2, color = .data[[color_col]], label=sample)) +
        geom_point(size=4, alpha=0.8) +
        geom_text_repel(max.overlaps = 10) +
        labs(title = "PCA Analysis (Imputed Data)", 
             subtitle = "Check for outliers here") +
        theme_bw() +
        theme(plot.title = element_text(face="bold", hjust=0.5))
      
      return(p)
    },

    # --- [新增] 功能 4: 自动化离群值检测 (PCA 或 相关性) ---
    detect_outliers = function(method = c("pca", "correlation"), 
                               sd_threshold = 3, 
                               remove = FALSE, 
                               do_plot = TRUE) {
      
      if (is.null(self$imputed_se)) stop("Please run perform_imputation() first.")
      method <- match.arg(method)
      
      mat <- assay(self$imputed_se)
      samples <- colnames(mat)
      outlier_samples <- c()
      stats_df <- data.frame(Sample = samples)
      
      message(paste0("--- Detecting Outliers using [", method, "] method (Threshold: ", sd_threshold, " SD) ---"))
      
      if (method == "pca") {
        # --- 策略 A: PCA 距离法 ---
        # 适用于检测整体表达分布偏离群体的样本
        pca_res <- prcomp(t(mat), scale. = TRUE)
        coords <- as.data.frame(pca_res$x[, 1:2]) # 取前两轴
        
        # 计算 PC1 和 PC2 的 Z-score
        coords$z_pc1 <- (coords$PC1 - mean(coords$PC1)) / sd(coords$PC1)
        coords$z_pc2 <- (coords$PC2 - mean(coords$PC2)) / sd(coords$PC2)
        
        # 定义离群: 任意一轴超过阈值
        is_outlier <- abs(coords$z_pc1) > sd_threshold | abs(coords$z_pc2) > sd_threshold
        outlier_samples <- rownames(coords)[is_outlier]
        
        stats_df$Score1 <- coords$PC1
        stats_df$Score2 <- coords$PC2
        stats_df$Is_Outlier <- is_outlier
        
        # 绘图
        p <- ggplot(stats_df, aes(x = Score1, y = Score2, color = Is_Outlier, label = Sample)) +
          geom_point(size = 3, alpha = 0.8) +
          geom_text_repel(data = subset(stats_df, Is_Outlier), max.overlaps = 20, size = 3) +
          stat_ellipse(level = 0.99, linetype = "dashed", color = "grey50") + # 99% 也就是大概 2.5-3 SD 范围
          scale_color_manual(values = c("FALSE"="#1465AC", "TRUE"="#B31B21")) +
          labs(title = "Outlier Detection (PCA)", 
               x = "PC1", y = "PC2", 
               subtitle = paste("Outliers:", length(outlier_samples))) +
          theme_bw()
        
      } else {
        # --- 策略 B: 平均相关性法 (Connectivity) ---
        # 适用于检测与大部队“格格不入”的样本 (如实验失败)
        
        # 计算相关矩阵
        cor_mat <- cor(mat, method = "pearson")
        
        # 计算每个样本与其他样本的平均相关系数 (Connectivity)
        # diag=NA 防止自己和自己相关性为1拉高均值
        diag(cor_mat) <- NA
        mean_cors <- colMeans(cor_mat, na.rm = TRUE)
        
        # 计算 Z-score (越低越差)
        z_scores <- (mean_cors - mean(mean_cors)) / sd(mean_cors)
        
        # 定义离群: 相关性显著过低 (单尾检验，只看负方向)
        is_outlier <- z_scores < -sd_threshold
        outlier_samples <- names(mean_cors)[is_outlier]
        
        stats_df$Mean_Cor <- mean_cors
        stats_df$Z_Score <- z_scores
        stats_df$Is_Outlier <- is_outlier
        
        # 绘图
        stats_df <- stats_df %>% arrange(Mean_Cor)
        stats_df$Sample <- factor(stats_df$Sample, levels = stats_df$Sample)
        
        p <- ggplot(stats_df, aes(x = Sample, y = Mean_Cor, color = Is_Outlier)) +
          geom_point(size = 3) +
          geom_segment(aes(x = Sample, xend = Sample, y = min(Mean_Cor)-0.02, yend = Mean_Cor)) +
          geom_hline(yintercept = mean(mean_cors) - sd_threshold * sd(mean_cors), 
                     linetype = "dashed", color = "red") +
          scale_color_manual(values = c("FALSE"="#1465AC", "TRUE"="#B31B21")) +
          labs(title = "Outlier Detection (Sample Connectivity)", 
               y = "Average Correlation with Peers", x = "Sample",
               subtitle = paste("Outliers (Z < -", sd_threshold, "): ", length(outlier_samples), sep="")) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
      }
      
      if (do_plot) print(p)
      
      # 执行删除
      if (length(outlier_samples) > 0) {
        message(paste0("  Identified ", length(outlier_samples), " outliers: ", paste(outlier_samples, collapse = ", ")))
        if (remove) {
          self$remove_outliers(outlier_samples) # 调用已有的删除函数
          message("  Outliers removed from imputed_se.")
        } else {
          message("  (Set remove = TRUE to verify and delete them)")
        }
      } else {
        message("  No outliers detected based on current threshold.")
      }
      
      return(list(outliers = outlier_samples, plot = p, stats = stats_df))
    },
    
    remove_outliers = function(outlier_samples_list) {
      if (is.null(self$imputed_se)) stop("Data not initialized.")

      if(is.null(self$imputed_se_backup)){
        self$imputed_se_backup <- self$imputed_se
      }

      current_samples <- colnames(self$imputed_se)
      valid_keep <- setdiff(current_samples, outlier_samples_list)

      if (length(valid_keep) < length(current_samples)) {
        # 1. 更新 imputed_se 对象
        self$imputed_se <- self$imputed_se[, valid_keep]

        # [新增] 2. 同时更新 se_obj 以保持一致性
        if (!is.null(self$se_obj)) {
          self$se_obj <- self$se_obj[, valid_keep]
        }

        # 3. 更新 meta_data (使用 base R 写法更安全，利用 rownames 匹配)
        self$meta_data <- self$meta_data[valid_keep, , drop = FALSE]

        message(paste("Removed Outliers:", paste(outlier_samples_list, collapse=", ")))
        message(paste("Remaining samples:", length(valid_keep)))
      } else {
        message("No matching outliers found to remove.")
      }
    },

    # --- [新增] 功能: 基于元数据筛选样本 (用于下游分析) ---
    subset_samples = function(column_name, values_to_keep) {
      if (is.null(self$imputed_se)) stop("Data not initialized or imputed.")
      
      message(paste("--- Subsetting Data based on:", column_name, "---"))
      
      current_meta <- colData(self$imputed_se)
      
      if (!column_name %in% colnames(current_meta)) {
        stop(paste("Error: Column '", column_name, "' not found.", sep=""))
      }
      
      keep_mask <- current_meta[[column_name]] %in% values_to_keep


      if (sum(keep_mask) == 0) {
        stop("Error: No samples found matching criteria.")
      }

      if(is.null(self$imputed_se_backup)){
        self$imputed_se_backup = self$imputed_se
      }

      if(is.null(self$meta_data_backup)){
        self$meta_data_backup = self$meta_data
      }

      
      # 直接修改 imputed_se and meta_data
      self$imputed_se <- self$imputed_se[, keep_mask]
      self$meta_data <- self$meta_data[keep_mask, ]
      
      message(paste("  Subset successful."))
      message(paste("  Kept values:", paste(values_to_keep, collapse=", ")))
      message(paste("  Remaining samples:", ncol(self$imputed_se)))
      message("  (Tip: Use reset_data() to restore the full dataset)")
    },

    reset_data = function() {
      if (is.null(self$imputed_se_backup)) {
        stop("No backup found. Please run perform_imputation() first.")
      }
      
      message("--- Resetting data to the state immediately after imputation ---")
      
      # 将当前工作对象恢复为备份对象
      self$imputed_se <- self$imputed_se_backup
      self$meta_data <- self$meta_data_backup
      
      # 同时尝试恢复 se_obj (虽然分析主要用 imputed_se，但保持一致性较好)
      # 注意：如果您的 se_obj 在 impute 前后没有结构变化，可以不恢复它，
      # 但为了安全，通常让 SE 也回到全量样本状态（尽管里面的值是未填补的）
      # 这里为了简单，我们主要保证 downstream 分析用的 imputed_se 是全量的
      
      message(paste("  Data restored. Current samples:", ncol(self$imputed_se)))
    },

    # --- [新增] 功能: 自动绘制蛋白表达趋势 (自动识别 连续vs分类 变量) ---
    plot_protein_expression = function(proteins, variable, color_var = NULL, add_labels = FALSE) {
      if (is.null(self$imputed_se)) stop("Data not ready. Please run perform_imputation() first.")
      
      # 1. 检查变量是否存在
      meta <- as.data.frame(colData(self$imputed_se))
      if (!variable %in% colnames(meta)) {
        stop(paste("Variable '", variable, "' not found in metadata.", sep=""))
      }
      
      # 2. 检查蛋白是否存在
      valid_prots <- intersect(proteins, rownames(self$imputed_se))
      if (length(valid_prots) == 0) stop("None of the specified proteins found in the dataset.")
      if (length(valid_prots) < length(proteins)) {
        warning(paste("Some proteins were not found:", paste(setdiff(proteins, valid_prots), collapse=", ")))
      }
      
      # 3. 准备绘图数据 (Long Format)
      # 使用 drop=FALSE 防止单蛋白时矩阵变成向量导致报错
      expr_mat <- assay(self$imputed_se)[valid_prots, , drop = FALSE]
      
      plot_df <- t(expr_mat) %>%
        as.data.frame() %>%
        rownames_to_column("sample_id") %>%
        left_join(meta %>% rownames_to_column("sample_id"), by = "sample_id") %>%
        pivot_longer(cols = all_of(valid_prots), names_to = "Protein", values_to = "Expression")
      
      # 4. 判断变量类型并绘图
      # 默认颜色变量与 x轴变量一致，除非用户指定
      if (is.null(color_var)) color_var <- variable
      
      x_vals <- plot_df[[variable]]
      
      # 初始化 ggplot
      p <- ggplot(plot_df, aes(x = .data[[variable]], y = Expression))
      
      # --- 逻辑分支 A: 连续变量 (Numeric) ---
      if (is.numeric(x_vals) && length(unique(x_vals)) > 2) {
        message(paste("Detected numeric variable '", variable, "'. Using Scatter plot + Loess.", sep=""))
        
        p <- p +
          geom_point(aes(color = .data[[color_var]]), size = 2.5, alpha = 0.7) +
          geom_smooth(method = "loess", color = "#B31B21", fill = "#B31B21", alpha = 0.2, linewidth = 1) +
          scale_color_viridis_c(option = "D", begin = 0.2, end = 0.8) + # 默认连续配色
          theme_bw()
          
      } else {
        # --- 逻辑分支 B: 分类变量 (Factor/Character) ---
        message(paste("Detected categorical variable '", variable, "'. Using Violin + Boxplot.", sep=""))
        
        # 强制转换为因子以确保离散绘图
        plot_df[[variable]] <- as.factor(plot_df[[variable]])
        
        p <- p +
          geom_violin(aes(fill = .data[[color_var]]), alpha = 0.5, trim = FALSE, scale = "width") +
          geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.9) +
          geom_jitter(width = 0.1, size = 1.5, alpha = 0.6, color = "grey30") +
          theme_bw()
      }
      
      # 5. 通用修饰
      # 如果有多个蛋白，使用分面显示
      if (length(valid_prots) > 1) {
        p <- p + facet_wrap(~Protein, scales = "free_y")
      } else {
        p <- p + labs(title = valid_prots[1])
      }
      
      # 添加样本标签 (可选)
      if (add_labels) {
         p <- p + ggrepel::geom_text_repel(aes(label = sample_id), size = 3, max.overlaps = 10)
      }
      
      p <- p + 
        labs(y = "Normalized Intensity (log2)", x = variable) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              strip.background = element_rect(fill = "white", color = "grey"),
              strip.text = element_text(face = "bold"),
              legend.position = "right")
      
      return(p)
    }
  )
)

# --- Class 2: 差异分析 (工厂模式：支持 DEP 和 Regression) ---
DiffExpAnalyst <- R6Class("DiffExpAnalyst",
  public = list(
    se_obj = NULL,
    diff_results = NULL, # 存储最终差异表
    sig_results = NULL,
    
    initialize = function(se_object) {
      self$se_obj <- se_object
    },
    
    # 模式 A: 分组比较 (利用 DEP 流程)
    # --- 核心修改：完全模拟 DEP test_diff 的行为，但支持协变量 ---
    run_dep_analysis = function(condition_col = "condition", control_group, case_group, covariates = NULL) {
      
      message(paste0("--- Running DEP-style Analysis with Covariates: ", case_group, " vs ", control_group, " ---"))
      
      # 1. 准备数据
      data_se <- self$se_obj
      meta_data <- colData(data_se)
      meta_data[[condition_col]] = as.factor(as.character(meta_data[[condition_col]]))
      
      # 2. 构建设计矩阵 (含协变量)
      # ---------------------------------------------------------
      # 确保分组变量是因子
      if(!is.factor(meta_data[[condition_col]])) {
         meta_data[[condition_col]] <- factor(meta_data[[condition_col]])
      }
      
      # 基础公式: ~ 0 + Group 
      formula_str <- paste0("~ 0 + ", condition_col)
      
      # 添加协变量
      if (!is.null(covariates)) {
        # 检查协变量是否存在
        missing_cov <- setdiff(covariates, colnames(meta_data))
        if (length(missing_cov) > 0) stop(paste("Covariates not found:", paste(missing_cov, collapse=", ")))
        
        formula_str <- paste0(formula_str, " + ", paste(covariates, collapse = " + "))
        message(paste("  Adjusting for:", paste(covariates, collapse=", ")))
      }
      
      # 生成设计矩阵
      design <- model.matrix(as.formula(formula_str), data = meta_data)
      
      # 3. 清洗列名 (适配 Limma Contrast)
      # ---------------------------------------------------------
      # model.matrix 会生成 "conditionD36" 这样的列名，我们需要改回 "D36"
      group_levels <- levels(meta_data[[condition_col]])
      current_cols <- colnames(design)
      
      for(lvl in group_levels) {
        # 正则匹配：^condition_col + level$ (例如 ^GroupD36$)
        pattern <- paste0("^", condition_col, lvl, "$")
        idx <- grep(pattern, current_cols)
        if(length(idx) > 0) current_cols[idx] <- lvl
      }
      colnames(design) <- current_cols
      
      # 4. 运行 Limma 分析
      # ---------------------------------------------------------
      contrast_formula <- paste0(case_group, " - ", control_group)
      contrast_name    <- paste0(case_group, "_vs_", control_group) # DEP 要求的命名格式


      
      message(paste("  Contrast:", contrast_formula))
      
      fit <- lmFit(assay(data_se), design)
      cont.matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
      fit2 <- contrasts.fit(fit, cont.matrix)
      fit2 <- eBayes(fit2)
      
      # 获取统计结果
      res_table <- topTable(fit2, number = Inf, sort.by = "none") # 必须 sort.by="none" 以保持顺序与 se 一致
      res_table$adj.P.Val = p.adjust(res_table$P.Value, method = "fdr")

      # 5. 关键步骤：构造符合 DEP 要求的 rowData 结构
      # ---------------------------------------------------------
      # DEP test_diff 会向 rowData 添加三列：
      # {contrast}_diff   (即 logFC)
      # {contrast}_p.val  (即 P.Value)
      # {contrast}_p.adj  (即 adj.P.Val)
      
      # 提取当前 rowData
      row_data <- rowData(data_se)
      
      # 准备新列 (列名必须严格符合 DEP 规范)
      col_diff <- paste0(contrast_name, "_diff")
      col_pval <- paste0(contrast_name, "_p.val")
      col_padj <- paste0(contrast_name, "_p.adj")
      
      # 将 Limma 结果注入 row_data
      # 确保行顺序一致 (topTable sort.by='none' 且 se未变动，通常是一致的，但为了安全用 match)
      mm <- match(rownames(row_data), rownames(res_table))
      
      row_data[[col_diff]] <- res_table$logFC[mm]
      row_data[[col_pval]] <- res_table$P.Value[mm]
      row_data[[col_padj]] <- res_table$adj.P.Val[mm]
      
      # 6. 更新并返回 SE 对象
      # ---------------------------------------------------------
      rowData(data_se) <- row_data
      
      # 存储一份 table 版本在内部 (可选，方便查看)
      self$diff_results <- res_table %>% rownames_to_column("Protein")
      
      message("  Analysis complete. Results added to SummarizedExperiment rowData.")
      message(paste("  Columns added:", col_diff, ",", col_pval, ",", col_padj))
      
      return(data_se) # 返回 SE 对象，以便后续直接接 add_rejections()
    },
    
    # 模式 B: 连续变量回归 (Linear or Spline)
    run_continuous_analysis = function(time_col = "time_num", method = "spline", df = 3) {
      message(paste("Running", method, "regression on", time_col))
      
      # 准备长数据
      expr_long <- assay(self$se_obj) %>%
        as.data.frame() %>%
        rownames_to_column("Protein") %>%
        pivot_longer(-Protein, names_to = "sample", values_to = "value") %>%
        left_join(as.data.frame(colData(self$se_obj)) %>% rownames_to_column("sample"), by = "sample")
      
      # 确保时间列是数值
      current_vals = expr_long[[time_col]]
      if (!is.numeric(current_vals)) {
        message(paste0("Notice: Column '", time_col, "' is not numeric (Type: ", class(current_vals)[1], "). Attempting to convert..."))
        
        # 1. 先转为字符
        vals_char <- as.character(current_vals)
        
        # 2. 使用 readr::parse_number 智能提取数字
        vals_num <- readr::parse_number(vals_char)
        
        # 3. 检查转换结果
        if (all(is.na(vals_num))) {
          stop(paste0("Error: Failed to convert column '", time_col, "' to numeric. ",
                      "It implies that the labels do not contain extractable numbers (e.g., 'Control', 'Treat')."))
        }
        
        # 4. 覆盖原列
        expr_long[[time_col]] <- vals_num
        message(paste0("  Conversion successful. Example: '", vals_char[1], "' -> ", vals_num[1]))
      }
      
      # 并行计算 (建议，因为对每个蛋白做回归很慢)
      results <- expr_long %>%
        group_by(Protein) %>%
        do({
          dat <- .
          if (method == "spline") {
            # Spline Regression (Helper 2 中的逻辑)
            fit_full <- lm(value ~ ns(get(time_col), df = df), data = dat)
            fit_null <- lm(value ~ 1, data = dat)
          } else {
            # Linear Regression
            fit_full <- lm(value ~ get(time_col), data = dat)
            fit_null <- lm(value ~ 1, data = dat)
          }
          
          lrt_res <- anova(fit_null, fit_full)
          model_stats <- glance(fit_full)
          
          # 计算相关性方向
          cor_res <- cor.test(dat[[time_col]], dat$value, method = "spearman", exact = FALSE)
          
          data.frame(
            p_value = lrt_res$`Pr(>F)`[2],
            adj_r_squared = model_stats$adj.r.squared,
            spearman_rho = cor_res$estimate,
            direction = ifelse(cor_res$estimate > 0, "Up", "Down")
          )
        }) %>%
        ungroup() %>%
        mutate(adj_pval = p.adjust(p_value, method = "BH")) %>%
        arrange(adj_pval)
      
      self$diff_results <- results
      return(results)
    },

    # 模式 C: 组间 ANOVA 分析 (多组比较)
    run_anova_analysis = function(condition_col = "condition", covariates = NULL, selected_groups = NULL) {
      message(paste("Running ANOVA analysis on", condition_col))

      # 1. 准备数据
      data_se <- self$se_obj
      meta_data <- colData(data_se)

      # 确保分组变量是因子
      if (!is.factor(meta_data[[condition_col]])) {
        meta_data[[condition_col]] <- as.factor(as.character(meta_data[[condition_col]]))
      }
      meta_data[[condition_col]] <- droplevels(meta_data[[condition_col]])

      # 2. 筛选特定组 (如果指定)
      all_group_levels <- levels(meta_data[[condition_col]])
      if (!is.null(selected_groups)) {
        invalid_groups <- setdiff(selected_groups, all_group_levels)
        if (length(invalid_groups) > 0) {
          stop(paste("Selected groups not found:", paste(invalid_groups, collapse = ", ")))
        }
        message(paste("  Filtering to selected groups:", paste(selected_groups, collapse = ", ")))
        
        keep_mask <- meta_data[[condition_col]] %in% selected_groups
        data_se <- data_se[, keep_mask]
        meta_data <- colData(data_se) # 更新 meta_data
        meta_data[[condition_col]] <- droplevels(meta_data[[condition_col]])
      }
      
      # 检查组数
      group_levels <- levels(meta_data[[condition_col]])
      n_groups <- length(group_levels)
      if (n_groups < 2) stop("ANOVA requires at least 2 groups.")
      message(paste("  Analyzing", n_groups, "groups using F-test."))

      # 3. 构建设计矩阵 (有截距模型 ~ Condition)
      # [修正点 1] 移除 '0 +'，使用默认的截距模型
      # 这样 coef=1 是截距(Reference)，coef=2..N 是各组相对于 Reference 的差异
      formula_str <- paste0("~ ", condition_col)

      if (!is.null(covariates)) {
        missing_cov <- setdiff(covariates, colnames(meta_data))
        if (length(missing_cov) > 0) stop(paste("Covariates not found:", paste(missing_cov, collapse=", ")))
        formula_str <- paste0(formula_str, " + ", paste(covariates, collapse = " + "))
        message(paste("  Adjusting for covariates:", paste(covariates, collapse=", ")))
      }

      design <- model.matrix(as.formula(formula_str), data = meta_data)

      # 4. 运行 Limma
      fit <- lmFit(assay(data_se), design)
      
      fit2 <- eBayes(fit)

      # 5. 提取 ANOVA 结果 (F-test)
      # [修正点 3] 只要检验除了 Intercept (第1列) 以外的与分组相关的列
      # 通常 condition 相关的系数是第 2 到 第 n_groups 列
      # 注意：如果有协变量，协变量的列会在 condition 之后，我们通常只检验 condition 的差异

      # 智能识别需要检验的系数 (列名中包含 condition_col 的列，但排除截距)
      coefs_to_test <- grep(condition_col, colnames(design))

      # 如果 grep 没找到 (比如列名被转义了)，回退到测试 2:n_groups
      if(length(coefs_to_test) == 0) {
         coefs_to_test <- 2:n_groups
      }

      # [关键修复] 再次检查 coefs_to_test 是否为空 (可能 n_groups < 2)
      if(length(coefs_to_test) == 0) {
        stop(paste("Cannot run ANOVA: Unable to identify coefficients for condition column '", condition_col, "'. Check if the condition column has at least 2 unique groups.", sep=""))
      }

      message(paste("  Testing coefficients columns:", paste(coefs_to_test, collapse=", ")))
      
      # 对选定的系数进行排序
      anova_res <- topTable(fit2, coef = coefs_to_test, number = Inf, sort.by = "p")

      # 整理结果
      res_table <- as.data.frame(anova_res) %>%
        rownames_to_column("Protein") %>%
        mutate(
          anova_F = F,
          anova_pval = P.Value # 这里的 P.Value 就是 F-test 的 P 值
        ) %>%
        dplyr::select(Protein, anova_F, P.Value, adj.P.Val, anova_pval)

      # 6. 存储结果
      self$diff_results <- res_table

      # 简单的结果诊断
      sig_count <- sum(res_table$adj.P.Val < 0.05, na.rm = TRUE)
      message("  ANOVA analysis complete.")
      message(paste("  Significant proteins (adj.P.Val < 0.05):", sig_count, "/", nrow(res_table)))
      
      # 如果依然全显著，打印警告
      if(sig_count > 0.8 * nrow(res_table)) {
        warning("  [Warning] >80% proteins are significant. Check if normalization was performed (e.g., log transformation).")
      }

      return(res_table)
    },
    # --- 2. 独立绘图函数: 绘制 ANOVA 差异热图 ---

    plot_anova_heatmap = function(se_obj, sig_proteins, 
                                  group_col, time_col = NULL,
                                  scale_rows = TRUE, 
                                  title = "ANOVA Significant Proteins") {
      
      message(paste("Drawing heatmap for", length(sig_proteins), "proteins..."))
      
      if(length(sig_proteins) == 0) {
        warning("No significant proteins to plot.")
        return(NULL)
      }
      
      if(length(sig_proteins) > 2000) {
        warning("Too many proteins (>2000). Selecting top 2000 based on variance for visibility.")
        # 如果蛋白太多，热图会看不清，这里增加一个自动下采样的保护机制
        mat_full <- assay(se_obj)[sig_proteins, ]
        vars <- apply(mat_full, 1, var)
        sig_proteins <- names(sort(vars, decreasing = TRUE))[1:2000]
      }
      
      # 1. 提取数据
      mat <- assay(se_obj)[sig_proteins, ]
      
      # 2. Z-score 标准化 (Row scaling)
      if(scale_rows) {
        mat <- t(scale(t(mat)))
      }
      
      # 3. 构建列注释 (样本信息)
      meta <- as.data.frame(colData(se_obj))
      
      # 准备注释数据框
      anno_df <- meta %>% dplyr::select(all_of(group_col))
      if(!is.null(time_col) && time_col %in% colnames(meta)) {
        anno_df[[time_col]] <- meta[[time_col]]
      }
      
      # 自动生成颜色
      # 为 Group 生成离散颜色
      groups <- unique(anno_df[[group_col]])
      group_colors <- setNames(ggsci::pal_npg()(length(groups)), groups)
      
      anno_colors <- list()
      anno_colors[[group_col]] <- group_colors
      
      # 创建 ComplexHeatmap 注释对象
      col_anno <- HeatmapAnnotation(
        df = anno_df,
        col = anno_colors,
        show_annotation_name = TRUE
      )
      
      # 4. 绘制热图
      ht <- Heatmap(
        mat,
        name = "Z-score",
        top_annotation = col_anno,

        # 行列聚类设置
        cluster_rows = TRUE,           # 聚类蛋白
        cluster_columns = TRUE,        # 聚类样本

        # [修复] 添加聚类距离方法和树状图设置
        clustering_distance_columns = "euclidean",  # 列聚类距离计算方法
        clustering_method_columns = "ward.D2",      # 列聚类方法

        # 显示列树状图
        show_column_dend = TRUE,
        column_dend_height = unit(10, "mm"),

        # 视觉调整
        show_row_names = length(sig_proteins) < 50, # 蛋白少时显示名字，多时不显示
        show_column_names = TRUE,

        # 颜色 (蓝-白-红)
        col = colorRamp2(c(-2, 0, 2), c("#3C5488B2", "white", "#E64B35B2")),

        column_title = title,

        # 自动切分行 (K-means): 自动把蛋白分成 2-4 类模式，便于观察
        row_km = ifelse(length(sig_proteins) > 100, 4, 1)
      )
      
      return(ht)
    },

    # 获取显著蛋白列表
    get_sig_proteins = function(pval_cutoff = 0.05,
                                    r2_cutoff = 0.5,
                                    corr_cutoff = 0.5,
                                    logfc_cutoff = 1,
                                    use_adjusted_pval = TRUE) {

          if (is.null(self$diff_results)) stop("No analysis results found. Run analysis first.")

          # 模式 1: ANOVA 分析模式 (检查 anova_pval 列)
          if ("anova_pval" %in% colnames(self$diff_results)) {
            message("Detected ANOVA analysis results.")

            if (use_adjusted_pval) {
              message(paste("Filtering ANOVA: adj.P.Val <", pval_cutoff))
              self$sig_results = self$diff_results %>%
                       filter(adj.P.Val < pval_cutoff)
            } else {
              message(paste("Filtering ANOVA: Raw P.Value <", pval_cutoff))
              self$sig_results = self$diff_results %>%
                       filter(P.Value < pval_cutoff)
            }

          # 模式 2: 连续变量回归模式 (检查 adj_r_squared 列)
          } else if ("adj_r_squared" %in% colnames(self$diff_results)) {
            message("Detected Continuous Regression results.")

            if (use_adjusted_pval) {
              message(paste("Filtering Regression: adj_pval <", pval_cutoff, "& |rho| >", corr_cutoff, "& R2 >", r2_cutoff))
              self$sig_results = self$diff_results %>%
                       filter(adj_pval < pval_cutoff & abs(spearman_rho) > corr_cutoff & adj_r_squared > r2_cutoff)
            } else {
              message(paste("Filtering Regression: Raw p_value <", pval_cutoff, "& |rho| >", corr_cutoff, "& R2 >", r2_cutoff))
              self$sig_results = self$diff_results %>%
                       filter(p_value < pval_cutoff & abs(spearman_rho) > corr_cutoff & adj_r_squared > r2_cutoff)
            }

          # 模式 3: 两组间比较模式 (DEP/Limma)
          } else {
            message("Detected Group Comparison results.")

            if (use_adjusted_pval) {
              message(paste("Filtering Group: adj.P.Val <", pval_cutoff, "& abs(logFC) >", logfc_cutoff))
              self$sig_results = self$diff_results %>%
                       filter(adj.P.Val < pval_cutoff & abs(logFC) > logfc_cutoff)
            } else {
              message(paste("Filtering Group: Raw p.val <", pval_cutoff, "& abs(logFC) >", logfc_cutoff))
              self$sig_results = self$diff_results %>%
                       filter(P.Value < pval_cutoff & abs(logFC) > logfc_cutoff)
            }
          }

          return(self$sig_results %>% pull(Protein))
    },

# --- [Intelligent Volcano Plot ---
    # 自动根据 diff_results 的内容判断是 Group 比较还是 Continuous 回归
    plot_volcano = function(logfc_cutoff = 1,       # 仅用于 Group 模式
                            corr_cutoff = 0.5,      # 仅用于 Regression 模式
                            r2_cutoff = 0.5,        # 仅用于 Regression 模式
                            pval_cutoff = 0.05,
                            top_n_labels = 10,
                            use_adjusted_pval = TRUE) {
      if(is.null(self$diff_results)) stop("No analysis results found. Please run run_dep_analysis() or run_continuous_analysis() first.")

      df <- self$diff_results

      # --- 0. 检查是否为 ANOVA 结果 (ANOVA 没有 volcano plot) ---
      if ("anova_F" %in% colnames(df)) {
        stop("Volcano plot is not available for ANOVA analysis. ANOVA compares multiple groups simultaneously and does not produce pairwise log fold changes. Please use the Heatmap or other visualization options.")
      }

      # --- 1. 判断分析类型并准备绘图数据 ---

      if ("adj_r_squared" %in% colnames(df)) {
        # >>> 模式 A: 连续变量回归 (Regression) <<<
        message("Detected Continuous Regression results. Plotting: X=Spearman Rho, Y=-log10(Pval)")

        # 确定 P 值列
        p_col <- if(use_adjusted_pval) "adj_pval" else "p_value"

        # 准备数据
        plot_data <- df %>%
          dplyr::select(Protein,
                        x_val = spearman_rho,      # X轴：相关系数 (表示方向)
                        y_val = !!sym(p_col),      # Y轴：P值
                        r2 = adj_r_squared) %>%    # 辅助：R2用于显著性判断
          filter(!is.na(x_val) & !is.na(y_val)) %>%
          mutate(
            log_pval = -log10(y_val),
            # 逻辑：P值显著 且 |Rho| 达标 且 R2 达标，根据 Rho 正负定上下调
            expression = case_when(
              y_val < pval_cutoff & abs(x_val) > corr_cutoff & r2 > r2_cutoff & x_val > 0 ~ "Up-regulated",
              y_val < pval_cutoff & abs(x_val) > corr_cutoff & r2 > r2_cutoff & x_val < 0 ~ "Down-regulated",
              TRUE ~ "Not Significant"
            )
          )

        # 设置坐标轴标签和阈值线位置
        x_lab <- "Spearman Correlation (rho)"
        y_lab <- if(use_adjusted_pval) expression(-log[10] ~ "Adjusted P-value") else expression(-log[10] ~ "P-value")
        v_lines <- c(-corr_cutoff, corr_cutoff) # 使用 corr_cutoff 作为垂直线阈值
        title_suffix <- "(Regression)"

      } else {
        # >>> 模式 B: 分组差异分析 (DEP/Limma) <<<
        message("Detected Group Comparison results. Plotting: X=LogFC, Y=-log10(Pval)")
        
        # 确定 P 值列 (run_dep_analysis 结果列名为 P.Value 和 adj.P.Val)
        p_col <- if(use_adjusted_pval) "adj.P.Val" else "P.Value"
        
        plot_data <- df %>%
          dplyr::select(Protein, 
                        x_val = logFC, 
                        y_val = !!sym(p_col)) %>%
          filter(!is.na(x_val) & !is.na(y_val)) %>%
          mutate(
            log_pval = -log10(y_val),
            expression = case_when(
              y_val < pval_cutoff & x_val >= logfc_cutoff ~ "Up-regulated",
              y_val < pval_cutoff & x_val <= -logfc_cutoff ~ "Down-regulated",
              TRUE ~ "Not Significant"
            )
          )
        
        x_lab <- expression(log[2] ~ "Fold Change")
        y_lab <- if(use_adjusted_pval) expression(-log[10] ~ "Adjusted P-value") else expression(-log[10] ~ "P-value")
        v_lines <- c(-logfc_cutoff, logfc_cutoff)
        title_suffix <- "(Group Comparison)"
      }
      
      # --- 2. 提取 Top Genes (用于标记) ---
      top_up <- plot_data %>%
        filter(expression == "Up-regulated") %>%
        arrange(desc(x_val)) %>% # 无论是 LogFC 还是 Rho，正向越大越显著
        slice_head(n = top_n_labels)
      
      top_down <- plot_data %>%
        filter(expression == "Down-regulated") %>%
        arrange(x_val) %>%       # 负向越小越显著
        slice_head(n = top_n_labels)
      
      top_genes <- bind_rows(top_up, top_down)
      
      # --- 3. 绘制火山图 ---
      p <- ggplot(plot_data, aes(x = x_val, y = log_pval)) +
        # 散点
        geom_point(aes(color = expression), alpha = 0.6, size = 1.5) +
        
        # 颜色定义
        scale_color_manual(values = c("Up-regulated" = "#B31B21", 
                                      "Down-regulated" = "#1465AC", 
                                      "Not Significant" = "grey80")) +
        
        # 阈值线
        geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "black", alpha = 0.5) +
        geom_vline(xintercept = v_lines, linetype = "dashed", color = "black", alpha = 0.5) +
        
        # 标签
        geom_text_repel(data = top_genes, 
                        aes(label = Protein),
                        size = 3,
                        box.padding = 0.5,
                        max.overlaps = Inf, 
                        segment.color = "grey50") +
        
        # 标题与轴
        labs(title = paste("Volcano Plot", title_suffix),
             subtitle = paste("Up:", nrow(filter(plot_data, expression == "Up-regulated")), 
                              "| Down:", nrow(filter(plot_data, expression == "Down-regulated"))),
             x = x_lab,
             y = y_lab) +
        
        # 主题
        theme_bw() +
        theme(panel.grid = element_blank(),
              legend.position = "bottom",
              legend.title = element_blank(),
              plot.title = element_text(face = "bold", hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5))
      
      return(p)
    }

  )
)

# --- Class 3: Mfuzz 聚类封装 ---
MfuzzClusterer <- R6Class("MfuzzClusterer",
  public = list(
    se_obj = NULL,
    mfuzz_set = NULL,
    cl_obj = NULL,
    enrich_worker = NULL,
    
    initialize = function(se_object, org_db = "org.Hs.eg.db") {
      self$se_obj <- se_object
      # 组合模式：MfuzzClusterer "拥有" 一个 EnrichmentAnalyst
      self$enrich_worker <- EnrichmentAnalyst$new(org_db)
    },
    
    run_mfuzz = function(sig_prots, time_col = "time_num", centers = 4) {
      library(Mfuzz)
      
      # 1. 准备均值数据 (按时间点取平均)
      meta <- as.data.frame(colData(self$se_obj))
      expr <- assay(self$se_obj)[sig_prots, ]
      
      df_mean <- as.data.frame(expr) %>%
        rownames_to_column("Protein") %>%
        pivot_longer(-Protein, names_to = "sample", values_to = "val") %>%
        left_join(meta %>% rownames_to_column("sample"), by="sample") %>%
        group_by(Protein, !!sym(time_col)) %>%
        summarise(mean_val = mean(val, na.rm=TRUE), .groups = "drop") %>%
        pivot_wider(names_from = !!sym(time_col), values_from = mean_val) %>%
        column_to_rownames("Protein") %>%
        as.matrix()
      
      # 2. 构建 ExpressionSet 并标准化
      self$mfuzz_set <- ExpressionSet(assayData = df_mean)
      self$mfuzz_set <- standardise(self$mfuzz_set)
      
      # 3. 聚类
      m <- mestimate(self$mfuzz_set)
      set.seed(123)
      self$cl_obj <- mfuzz(self$mfuzz_set, c = centers, m = m)
      message("Mfuzz clustering complete.")
    },


    # --- 核心功能: 自动完成富集并绘图 ---
    plot_clusters = function(mfrow = c(2,2)) {
    mfuzz.plot2(self$mfuzz_set, cl = self$cl_obj, min.mem = 0.7, mfrow = mfrow, 
                time.labels = colnames(exprs(self$mfuzz_set)), colo = "fancy", x11=FALSE)
    },
    
    get_cluster_df = function(min.acore_threshold = 0.7) {
      acore(self$mfuzz_set, self$cl_obj, min.acore = min.acore_threshold) %>% 
        bind_rows(.id="Cluster") %>% 
        rename(Gene = NAME)
    },

    # --- 绘图逻辑: 趋势图 + ORA 气泡图 ---
    draw_composite_plot = function(cluster_id, gene_list,  time_col = time_col, ora_res, db = c("GO_BP", "KEGG")) {
      # 1. 趋势图 (Loess)      
      meta <- as.data.frame(colData(self$se_obj))
      expr <- assay(self$se_obj)[gene_list, ]
      
      plot_df <- as.data.frame(expr) %>%
        rownames_to_column("Protein") %>%
        pivot_longer(-Protein, names_to = "sample", values_to = "val") %>%
        left_join(meta %>% rownames_to_column("sample"), by="sample") %>%
        group_by(Protein, !!sym(time_col)) %>%
        summarise(mean_val = mean(val, na.rm=TRUE), .groups = "drop") 


      p_trend <- ggplot(plot_df, aes(x = !!sym(time_col), y = mean_val)) +
        geom_line(aes(group = Protein), color = "grey70", alpha = 0.2, linewidth = 0.3) +
        geom_smooth(aes(group = 1), method = "loess", color = "#B31B21", fill = "#B31B21", alpha = 0.2, linewidth = 1.5) +
        labs(title = paste("Cluster", cluster_id, "Trend"), x = "Time", y = "Z-score") +
        theme_bw() + theme(panel.grid = element_blank())
      
      # 2. db1 气泡图
      p_db1 <- if(!is.null(ora_res[[db[1]]]) && nrow(ora_res[[db[1]]]) > 0) {
        clusterProfiler::dotplot(ora_res[[db[1]]], showCategory = 10, color = "pvalue") + ggtitle(db[1]) + theme(axis.text.y = element_text(size = 8))
      } else { ggplot() + theme_void() + annotate("text", x=0.5, y=0.5, label="No Sig pathways") }
      
      # 3. db2 气泡图
      p_db2 <- if(!is.null(ora_res[[db[2]]]) && nrow(ora_res[[db[2]]]) > 0) {
        clusterProfiler::dotplot(ora_res[[db[2]]], showCategory = 10, color = "pvalue") + ggtitle(db[2]) + theme(axis.text.y = element_text(size = 8))
      } else { ggplot() + theme_void() + annotate("text", x=0.5, y=0.5, label="No Sig pathways") }
      
      # 4. 拼图
      layout <- "
      AABBB
      AACCC
      "
      p = p_trend + p_db1 + p_db2 + patchwork::plot_layout(design = layout) + 
        plot_annotation(title = paste("Cluster", cluster_id, "Analysis"))

      return(p)
    },

    analyze_and_plot_clusters = function(min.acore_threshold = 0.7, view_db = c("GO_BP", "KEGG"), pvalueCutoff = 0.05, time_col = "time_num", target_dir = NULL) {
      
      if(is.null(self$cl_obj)) stop("Run run_mfuzz() first.")
      message("=== Mfuzz: Internal Enrichment & Plotting ===")
      
      # 1. 获取聚类信息
      cluster_df <- self$get_cluster_df(min.acore_threshold = min.acore_threshold)
      
      clusters <- unique(cluster_df$Cluster)
      excel_list <- list()
      plot_list <- list()
      
      # 2. 循环处理每个 Cluster
        if (!dir.exists(target_dir)) {
          dir.create(target_dir)
        }

      for(cl in clusters) {
        message(paste("  Processing Cluster", cl, "..."))
        target_genes <- cluster_df %>% filter(Cluster == cl) %>% pull(Gene)
        
        # [关键] 调用内部 worker 进行分析
        universe_genes = rownames(self$se_obj)
        ora_res <- self$enrich_worker$run_comprehensive_ora(target_genes, universe_genes, pval_cutoff = pvalueCutoff)
        
        # 保存表格

        self$enrich_worker$enrich_to_excel(ora_res, direction = "Mfuzz", output_prefix = paste0("Mfuzz_cluster", "_", cl), target_dir = target_dir)
        
        # 绘制组合图
        p_combined <- self$draw_composite_plot(cl, target_genes, time_col, ora_res, db = view_db)
        plot_list[[as.character(cl)]] <- p_combined
      }

        # 保存图片
        pdf(here::here(target_dir, "Mfuzz_Cluster_Analysis_Summary.pdf"), width = 14, height = 10)
        for(p in plot_list) { print(p) }
        dev.off()

        return(plot_list)
      }

  )
)


# --- Class 4: 富集分析封装 ---
EnrichmentAnalyst <- R6Class("EnrichmentAnalyst",
  public = list(
    organism_db = NULL,
    ora_up = NULL,
    ora_down = NULL,
    gsea_res = NULL,


    initialize = function(org_db = "org.Hs.eg.db") {
      self$organism_db <- org_db
    },
    
    # 静态辅助方法：ID 转换
    convert_to_entrez = function(gene_symbols) {
      fix_map <- c(
        "GATD3B" = "GATD3",
        "SMAP" = "KIFAP3",
        "MT-CYB" = "CYTB",
        "MT-CO1" = "COX1",
        "MT-CO2" = "COX2",
        "MT-CO3" = "COX3",
        "MT-ND1" = "ND1",
        "MT-ND2" = "ND2",
        "MT-ND3" = "ND3",
        "MT-ND4" = "ND4",
        "MT-ND5" = "ND5",
        "MT-ATP8" = "ATP8",
        "PRPF4B" = "PRP4K",
        "BAP18" = "BACC1",
        "KCT2" = "C5orf15",
        "CUSTOS" = "C12orf43",
        "SARG" = "C1orf116",
        "IMUP" = "C19orf33",
        "C12orf4" = "FERRY3",
        "BARGIN" = "SH3BP1",
        "NMES1" = "COXFA4L3"

        )
      
      genes_fixed <- gene_symbols
      for(old_name in names(fix_map)) {
        genes_fixed[genes_fixed == old_name] <- fix_map[old_name]
      }

      bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = self$organism_db)$ENTREZID
    },
    
    # 综合 ORA 分析 (Top 10 上下调)
    run_comprehensive_ora = function(gene_list, universe_genes, pval_cutoff, padjust_method = "none", qval_cutoff = 1) {
      # gene_list: 字符向量 (Gene Symbols)
      
      entrez_ids <- self$convert_to_entrez(gene_list)
      universe_entrez <- self$convert_to_entrez(universe_genes)
      
      res_list <- list()
      
      # 1. GO (BP, MF, CC)
      for (ont in c("BP", "MF", "CC")) {
        res_list[[paste0("GO_", ont)]] <- enrichGO(
          gene = entrez_ids, universe = universe_entrez, OrgDb = self$organism_db,
          ont = ont, readable = TRUE, pvalueCutoff = pval_cutoff, pAdjustMethod = padjust_method, qvalueCutoff = qval_cutoff)
      }
      
      # 2. KEGG
      res_list[["KEGG"]] <- enrichKEGG(gene = entrez_ids, organism = "hsa", universe = universe_entrez, 
        pvalueCutoff = pval_cutoff, pAdjustMethod = padjust_method, qvalueCutoff = qval_cutoff) 

      # SetReadable later
      if (!is.null(res_list[["KEGG"]])) {
      res_list[["KEGG"]] <- setReadable(res_list[["KEGG"]], OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      }


      # 3. Reactome
      res_list[["Reactome"]] <- enrichPathway(gene = entrez_ids, organism = "human", universe = universe_entrez, readable = TRUE, 
        pvalueCutoff = pval_cutoff, pAdjustMethod = padjust_method, qvalueCutoff = qval_cutoff)

      
      # 4. Wiki (Need logical check if installed, skip if not)
       res_list[["Wiki"]] <- enrichWP(gene = entrez_ids, universe = universe_entrez, organism = "Homo sapiens", 
        pvalueCutoff = pval_cutoff, pAdjustMethod = padjust_method, qvalueCutoff = qval_cutoff) 
    
      if (!is.null(res_list[["Wiki"]])) {
        res_list[["Wiki"]] <- setReadable(res_list[["Wiki"]], OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      }

      return(res_list)
    },
    
    # GSEA 分析
    run_comprehensive_gsea = function(ranked_gene_list, pval_cutoff, padjust_method = "none") {
      # ranked_gene_list: named vector (names=Symbol, values=stat/logFC), sorted desc
      
      ids <- bitr(names(ranked_gene_list), fromType="SYMBOL", toType="ENTREZID", OrgDb=self$organism_db)
      # Merge LogFC
      df_merge <- data.frame(SYMBOL = names(ranked_gene_list), val = ranked_gene_list) %>%
        inner_join(ids, by="SYMBOL") %>%
        arrange(desc(val))
      
      gene_vec <- df_merge$val
      names(gene_vec) <- df_merge$ENTREZID
      
      gsea_res <- list()
      gsea_res[["GO_BP"]] <- gseGO(gene_vec, OrgDb=self$organism_db, ont="BP", pvalueCutoff = pval_cutoff, pAdjustMethod = padjust_method, minGSSize = 10, maxGSSize = 500, seed=123)
      gsea_res[["GO_MF"]] <- gseGO(gene_vec, OrgDb=self$organism_db, ont="MF", pvalueCutoff = pval_cutoff, pAdjustMethod = padjust_method, minGSSize = 10, maxGSSize = 500, seed=123)
      gsea_res[["GO_CC"]] <- gseGO(gene_vec, OrgDb=self$organism_db, ont="CC", pvalueCutoff = pval_cutoff, pAdjustMethod = padjust_method, minGSSize = 10, maxGSSize = 500, seed=123)
      gsea_res[["Wiki"]] <-  gseWP(gene_vec, organism = "Homo sapiens", pvalueCutoff = pval_cutoff, pAdjustMethod = padjust_method, minGSSize = 10, seed = 123)
      gsea_res[["KEGG"]] <- gseKEGG(gene_vec, organism="hsa", pvalueCutoff = pval_cutoff, pAdjustMethod = padjust_method, minGSSize = 10, seed=123)
      gsea_res[["Reactome"]] <- gsePathway(gene_vec, organism="human", pvalueCutoff = pval_cutoff, pAdjustMethod = padjust_method, minGSSize = 10, seed=123)

      return(gsea_res)
    },

    # --- 业务逻辑: 处理 DiffExpAnalyst 对象 ---
    analyze_diff_obj = function(diff_obj, pval_cutoff = 1) {
      if(is.null(diff_obj$sig_results) || is.null(diff_obj$diff_results)) stop("Run DE analysis first.")
      
      message("=== EnrichmentAnalyst: Processing DiffExp Object ===")
      
      # 1. 提取数据
      universe <- diff_obj$diff_results$Protein
      sig_df <- diff_obj$sig_results
      
      # 2. 区分上下调
      if ("logFC" %in% colnames(sig_df)) {
        up_genes <- sig_df %>% filter(logFC > 0) %>% pull(Protein)
        down_genes <- sig_df %>% filter(logFC < 0) %>% pull(Protein)
        rank_vec <- setNames(diff_obj$diff_results$logFC, diff_obj$diff_results$Protein)
      } else { # 假设是相关性
        up_genes <- sig_df %>% filter(spearman_rho > 0) %>% pull(Protein)
        down_genes <- sig_df %>% filter(spearman_rho < 0) %>% pull(Protein)
        rank_vec <- setNames(diff_obj$diff_results$spearman_rho, diff_obj$diff_results$Protein)
      }
      
      # 3. 运行 ORA
      self$ora_up <- self$run_comprehensive_ora(up_genes, universe, pval_cutoff = pval_cutoff)
      self$ora_down <- self$run_comprehensive_ora(down_genes, universe, pval_cutoff = pval_cutoff)
      
      # 4. 运行 GSEA
      self$gsea_res <- self$run_comprehensive_gsea(sort(rank_vec, decreasing = TRUE), pval_cutoff = pval_cutoff)      
      message("DiffExp Enrichment Done.")
    },

# 修改后的 enrich_to_excel 函数
    enrich_to_excel = function(enrich_obj = NULL, direction = "UP", output_prefix = "DiffExp", target_dir = NULL){
      
      # --- 逻辑 1: 自动根据 direction 选择内部对象 ---
      if(is.null(enrich_obj)){
        if(toupper(direction) == "UP"){
          enrich_obj <- self$ora_up
          message("Selected self$ora_up for export.")
        } else if (toupper(direction) == "DOWN"){
          enrich_obj <- self$ora_down
          message("Selected self$ora_down for export.")
        } else {
          stop("Direction must be 'UP' or 'DOWN' when enrich_obj is not provided.")
        }
      }
      
      # 检查对象是否为空（可能尚未运行分析）
      if(is.null(enrich_obj)) {
        warning(paste("No enrichment results found for direction:", direction, "- Skipping Excel export."))
        return(NULL)
      }

      # --- 逻辑 2: 准备 Excel 数据 ---
      excel_list <- list()
      databases = c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Wiki", "Reactome")

      for(db in databases){
        # 列表名增加 direction 前缀，如 UP_GO_BP
        list_name = paste0(direction, "_", db)
        
        # 检查该数据库是否有结果
        if(!is.null(enrich_obj[[db]]) && nrow(enrich_obj[[db]]) > 0) {
          excel_list[[list_name]] <- as.data.frame(enrich_obj[[db]])
        }
      }
      
      if(length(excel_list) == 0) {
        warning("Enrichment object is empty (no significant pathways). No Excel file created.")
        return(NULL)
      }
      
      # --- 逻辑 3: 处理保存路径与 target_dir ---
      # 构建文件名 (建议加上 direction 防止覆盖)
      file_name <- paste0(output_prefix, "_", direction, ".xlsx")
      
      final_path <- file_name
      if(!is.null(target_dir)) {
        # 如果目录不存在，自动创建
        if(!dir.exists(target_dir)) {
          dir.create(target_dir, recursive = TRUE)
          message(paste("Created directory:", target_dir))
        }
        final_path <- here::here(target_dir, file_name)
      }
      
      # 保存
      write.xlsx(excel_list, file = final_path, overwrite = TRUE)
      message(paste("Enrichment results saved to:", final_path))
    },

    gsea_to_excel = function(enrich_obj = NULL, output_prefix = "DiffExp", target_dir = NULL){
      # --- 逻辑 1: 自动选择内部 GSEA 对象 ---
      if(is.null(enrich_obj)){
        enrich_obj <- self$gsea_res
        message("Selected self$gsea_res for export.")
      }
      
      # 检查是否依然为空 (既没传参，内部也没算过)
      if(is.null(enrich_obj) ) {
        stop("GSEA results are empty. Please run analyze_diff_obj() first or provide an enrichment object.")
      }

      # --- 逻辑 2: 准备 Excel 数据 ---
      databases = c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Wiki", "Reactome")
      excel_gsea_list <- list()

      for(db in databases){
        list_name = paste0("GSEA", "_", db)
        # 增加 nrow > 0 的判断，防止写入空 sheet
        if(!is.null(enrich_obj[[db]]) && nrow(enrich_obj[[db]]) > 0) {
          excel_gsea_list[[list_name]] <- as.data.frame(enrich_obj[[db]])
        }
      }
      
      if(length(excel_gsea_list) == 0) {
        warning("GSEA object is empty (no significant pathways found). No Excel file created.")
        return(NULL)
      }

      # --- 逻辑 3: 处理保存路径与 target_dir ---
      file_name <- paste0(output_prefix, "_GSEA.xlsx")
      final_path <- file_name
      
      if(!is.null(target_dir)) {
        # 如果目录不存在，自动创建
        if(!dir.exists(target_dir)) {
          dir.create(target_dir, recursive = TRUE)
          message(paste("Created directory:", target_dir))
        }
        final_path <- file.path(target_dir, file_name)
      }

      # 保存文件
      write.xlsx(excel_gsea_list, file = final_path, overwrite = TRUE)
      message(paste("GSEA results saved to:", final_path))
    },


    draw_enrich_plot = function(enrich_obj = NULL, db_name, top_n = 10) {

      if (is.null(enrich_obj[[db_name]]) || nrow(enrich_obj[[db_name]]) == 0) {
        # Return placeholder if empty
        return(ggdraw() + 
                 draw_label(paste("No significant", db_name, "\npathways found"), 
                            fontface = 'italic', size = 10))
      } else {
        # Create dotplot
        # We create a title like "GO (BP)"
        p <- clusterProfiler::dotplot(enrich_obj[[db_name]], showCategory = top_n, color = "pvalue") + 
          ggtitle(db_name) +
          theme(axis.text.y = element_text(size = 8),
                plot.title = element_text(size = 11, face = "bold"))
        return(p)
      }
    },

    draw_gsea_plot = function(enrich_obj, db_name, top_n = 15) {

      p <- clusterProfiler::dotplot(enrich_obj[[db_name]], showCategory = 15, split = ".sign", label_format = 50, color = "pvalue") + 
      facet_grid(.~.sign) +
      ggtitle(paste("GSEA:", db_name)) +
      theme(plot.title = element_text(size=10))

      return(p)
    },

    # --- [新增] 功能: 绘制特定 Pathway 的蛋白表达趋势图 ---
    plot_pathway_trend = function(pathway_id,    # Pathway ID (如 "hsa04110") 或 Description (部分匹配)
                                  se_obj,        # 必须传入包含表达矩阵的 SE 对象
                                  meta_col,      # metadata 中的列名 (如 "age", "stage")
                                  direction = "UP", # "UP" 或 "DOWN"
                                  db_name = "KEGG"  # "GO_BP", "KEGG", "Reactome" 等
                                  ) {
      
      # 1. 获取富集结果对象
      target_res_list <- if (toupper(direction) == "UP") self$ora_up else self$ora_down
      
      if (is.null(target_res_list) || is.null(target_res_list[[db_name]])) {
        stop(paste0("No enrichment results found for direction: ", direction, " and DB: ", db_name))
      }
      
      res_df <- as.data.frame(target_res_list[[db_name]])
      
      # 2. 查找 Pathway
      # 尝试精确匹配 ID，如果不行尝试匹配 Description
      target_row <- res_df %>% filter(ID == pathway_id)
      
      if (nrow(target_row) == 0) {
        # 尝试 Description 模糊匹配
        target_row <- res_df %>% filter(grepl(pathway_id, Description, ignore.case = TRUE))
      }
      
      if (nrow(target_row) == 0) {
        stop(paste0("Pathway '", pathway_id, "' not found in ", direction, " ", db_name, " results."))
      } else if (nrow(target_row) > 1) {
        warning(paste0("Multiple pathways matched '", pathway_id, "'. Using the first one: ", target_row$Description[1]))
        target_row <- target_row[1, ]
      }
      
      pathway_name <- target_row$Description
      gene_str <- target_row$geneID # 假设是 Symbol (readable=TRUE)
      
      message(paste0("Plotting trend for: ", pathway_name, " (", target_row$ID, ")"))
      
      # 3. 解析基因并提取表达数据
      target_genes <- unlist(strsplit(gene_str, "/")) # clusterProfiler 默认用 "/" 分隔
      
      # 检查 SE 对象中的基因
      expr_mat <- assay(se_obj)
      valid_genes <- intersect(target_genes, rownames(expr_mat))
      
      if (length(valid_genes) == 0) stop("None of the pathway genes found in the expression matrix.")
      if (length(valid_genes) < 3) warning("Less than 3 genes found for this pathway. Trends might be unstable.")
      
      # 4. Z-score 标准化 (按行/蛋白 Scale)
      # t(apply(..., 1, scale)) 使得每个蛋白在样本间均值为0，方差为1
      subset_mat <- expr_mat[valid_genes, , drop=FALSE]
      scaled_mat <- t(apply(subset_mat, 1, scale))
      colnames(scaled_mat) <- colnames(expr_mat) # scale 会丢失列名，补回
      
      # 5. 整合 Metadata
      meta <- as.data.frame(colData(se_obj))
      
      if (!meta_col %in% colnames(meta)) stop(paste0("Column '", meta_col, "' not found in metadata."))
      
      plot_df <- as.data.frame(scaled_mat) %>%
        rownames_to_column("Gene") %>%
        pivot_longer(cols = -Gene, names_to = "SampleID", values_to = "Z_Score")
      
      # 添加 metadata 信息
      # 确保 SampleID 能匹配到 metadata 的行名或 ID 列
      if(!all(plot_df$SampleID %in% rownames(meta))) {
         # 尝试匹配 metadata 中的 ID 列（如果有的话），否则假设 rowname 就是 sample id
         # 这里简单处理：假设 metadata 的 rowname 就是 sample id
         warning("Sample IDs in matrix do not fully match metadata rownames. Attempting matching...")
      }
      
      plot_df$MetaVal <- meta[plot_df$SampleID, meta_col]
      
      # 6. 绘图 (自动判断 X 轴类型)
      p <- ggplot(plot_df, aes(x = MetaVal, y = Z_Score))
      
      if (is.numeric(plot_df$MetaVal)) {
        # --- 连续变量 (如 Age): 散点 + Loess 曲线 ---
        p <- p +
          # 背景线条：每个基因的独立趋势
          geom_line(aes(group = Gene), alpha = 0.15, color = "grey60", linewidth = 0.3) +
          # 整体趋势线
          geom_smooth(method = "loess", color = "#B31B21", fill = "#B31B21", alpha = 0.2, linewidth = 1.5) +
          # 散点
          geom_point(alpha = 0.3, size = 1, color = "#1465AC") +
          labs(x = meta_col)
          
      } else {
        # --- 分类变量 (如 Stage): 箱线图 ---
        p <- p +
          geom_violin(aes(fill = MetaVal), alpha = 0.2, color = NA) +
          geom_boxplot(aes(color = MetaVal), width = 0.2, fill = "white", outlier.shape = NA) +
          geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) +
          scale_color_brewer(palette = "Set1") +
          scale_fill_brewer(palette = "Set1") +
          labs(x = meta_col)
      }
      
      # 7. 通用修饰
      p <- p +
        theme_minimal(base_size = 14) +
        labs(
          title = stringr::str_trunc(paste("Pathway Trend:", pathway_name), 60),
          subtitle = str_wrap(paste("Based on", length(valid_genes), "significant genes: ", paste(valid_genes, collapse = "; ")), width = 80),
          y = "Normalized Expression (Z-score)"
        ) +
        theme(
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 10, color = "grey30"),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black")
        )
      
      return(p)
    }

  )
)
