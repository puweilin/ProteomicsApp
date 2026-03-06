# =============================================================================
# ProteomicsAnalysis: Core R6 classes for proteomics data analysis
# =============================================================================
# Note: All packages are loaded via app.R setup_environment().
# Script files use namespace-qualified calls (pkg::fn) for clarity.
# =============================================================================

# --- Class 1: Data management and preprocessing ---
# --- Class 1: Data management and preprocessing (Updated with Missing Filter) ---
ProteomicsDataManager <- R6::R6Class("ProteomicsDataManager",
  public = list(
    raw_mat = NULL,
    meta_data = NULL,
    meta_data_backup = NULL,
    annot_data = NULL,
    se_obj = NULL,
    imputed_se = NULL,
    imputed_se_backup = NULL,
    valid_protein_numbers = NULL,
    tolerate_missing_percent = NULL,  # Missing value tolerance (0~1), default 0.5 = allow 50% missing, filter above
    imputation_method = NULL,  # Store the imputation method used
    missing_mask = NULL,  # Store missing value mask before imputation (TRUE = originally missing)
    
    # --- Initialize ---
    initialize = function(mat_file, meta_file, annot_file, tolerate_missing_percent) {
      self$raw_mat <- mat_file
      self$meta_data <- meta_file
      self$annot_data <- annot_file
      self$tolerate_missing_percent <- tolerate_missing_percent
    },
    
    # --- Core: Data preprocessing and SE object construction ---
    # [Modified] Support filtering by any metadata column
    process_data = function(filter_col = NULL, filter_value = NULL) {

message("--- Step 1: Data Filtering & Matching ---")

      # 1. Identify matrix structure
      mat_col_names <- colnames(self$raw_mat)
      id_col_name   <- mat_col_names[1]
      all_sample_cols <- mat_col_names[-1]

      # 2. Filter metadata (supports any column)
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
      
      # 3. Take intersection
      valid_samples <- intersect(all_sample_cols, target_meta$label)
      if (length(valid_samples) == 0) stop("Error: No common samples found.")
      
      message(paste("  Selected", length(valid_samples), "samples."))
      
      # 4. Align data
      self$raw_mat <- self$raw_mat %>% dplyr::select(all_of(c(id_col_name, valid_samples)))
      self$meta_data <- target_meta %>% dplyr::filter(label %in% valid_samples)
      
      # 5. ID conversion
      message("--- Step 2: Mapping Protein IDs to Gene Names ---")
      annot_join_col <- colnames(self$annot_data)[1]
      
      # Auto-detect Gene Name column in annot_data (look for "Gene" keyword, fallback to 2nd column)
      gene_name_col <- grep("Gene", colnames(self$annot_data), value = TRUE, ignore.case = TRUE)[1]
      if(is.na(gene_name_col)) gene_name_col <- colnames(self$annot_data)[2]
      
      # Rename for consistent processing
      annot_clean <- self$annot_data %>% 
        dplyr::select(all_of(c(annot_join_col, gene_name_col))) %>%
        setNames(c(annot_join_col, "Gene.Name.Mapped"))
        
      merged_data <- self$raw_mat %>%
        left_join(annot_clean, by = setNames(annot_join_col, id_col_name)) %>%
        filter(!is.na(Gene.Name.Mapped) & Gene.Name.Mapped != "")
      
      # Handle duplicate gene names
      unique_genes <- make.unique(as.character(merged_data$Gene.Name.Mapped))
      
      # 6. Build initial expression matrix (for SE construction)
      # Build a data.frame conforming to DEP make_se requirements
      # Must contain "name" (Gene Name) and "ID" (Protein ID) columns
      
      df_for_se <- merged_data %>%
        dplyr::select(all_of(valid_samples)) %>%
        mutate(name = unique_genes, ID = merged_data[[id_col_name]]) %>%
        dplyr::select(name, ID, everything())
      
      # 7. Create SummarizedExperiment object
      # columns should be sample column indices
      sample_cols_indices <- which(colnames(df_for_se) %in% valid_samples)
      
      self$se_obj <- make_se(df_for_se, columns = sample_cols_indices, expdesign = self$meta_data) %>%
        normalize_vsn()
      
      # Update meta_data to SE colData (ensure consistent order)
      self$meta_data <- colData(self$se_obj) %>% data.frame()
      
      # --- [Fixed] Step 2.5: Filter missing values at SE object level ---
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
    
    # --- Feature: Plot missing value pattern ---
    plot_missing_pattern = function() {
      if (is.null(self$se_obj)) stop("Please run process_data() first.")
      message("--- Plotting Missing Value Pattern using DEP::plot_missval ---")
      p <- plot_missval(self$se_obj)
    },

perform_imputation = function(method = c("missForest","bpca", "knn", "QRILC", "MLE", "MinDet", "MinProb",
      "man", "min", "zero", "mixed", "nbavg"), cores = 4,
      ntree = 100, maxiter = 10) {
      
      if (is.null(self$se_obj)) stop("Please run process_data() first.")
      method <- match.arg(method)

      message(paste("--- Step 3: Running Imputation with Method:", method, "---"))

      data_matrix <- assay(self$se_obj)

      # [Added] Save missing value mask before imputation
      self$missing_mask <- is.na(data_matrix)
      message(paste0("  Detected ", sum(self$missing_mask), " missing values (",
                    round(sum(self$missing_mask) / prod(dim(data_matrix)) * 100, 2), "%)"))

      if (method == "missForest") {
        # Dependency check
        if (!requireNamespace("missForest", quietly = TRUE)) stop("Package 'missForest' is required.")
        if (!requireNamespace("doParallel", quietly = TRUE)) stop("Package 'doParallel' is required.")
        if (!requireNamespace("parallel", quietly = TRUE)) stop("Package 'parallel' is required.")

        # missForest requires rows=observations(Sample), cols=variables(Protein), so transpose first
        data_matrix_t <- t(data_matrix)
        
        # --- [Fixed] Platform-specific parallel processing ---
        message(paste0("  Setting up parallel backend for ", .Platform$OS.type, " with ", cores, " cores..."))
        
        if (.Platform$OS.type == "windows") {
          # Windows: must use Socket Cluster (PSOCK)
          cl <- parallel::makePSOCKcluster(cores)
          doParallel::registerDoParallel(cl)
          
          # Ensure cluster is closed and memory released when function exits
          on.exit(parallel::stopCluster(cl), add = TRUE)
          
        } else {
          # Mac/Linux: can use Fork mechanism directly
          doParallel::registerDoParallel(cores)
        }
        
        message("  Running missForest (this may take time)...")
        set.seed(123)
        # parallelize = 'forests' is typically faster and more robust than 'variables'
        mf_result <- missForest::missForest(data_matrix_t,
                                             ntree = ntree,
                                             maxiter = maxiter,
                                             verbose = TRUE,
                                             parallelize = "forests")
        
        # Transpose back to SE format (rows=Protein, cols=Sample)
        imputed_matrix <- t(mf_result$ximp)

      } else {
        # DEP package other imputation methods
        # fun parameter options: "MinProb", "MinDet", "QRILC", "Man", "Min" etc.
        imputed_matrix <- assay(DEP::impute(self$se_obj, fun = method))
      }

      self$imputed_se <- self$se_obj
      self$imputed_se_backup <- self$se_obj

      assay(self$imputed_se) <- imputed_matrix
      assay(self$imputed_se_backup) <- imputed_matrix
      self$meta_data_backup  <- self$meta_data
      self$imputation_method <- method  # [Added] Store the method used

      message("--- Imputation Complete ---")
      message("--- Backup created. You can restore full data using reset_data() ---")

      invisible(list(method = method))
    },

# --- [Fixed] Feature 2.5: Assess imputation quality (fixed plot column name error + dimension mismatch after outlier removal) ---
    assess_imputation = function(plot_type = c("structure", "density"),
                                 check_integrity = TRUE,
                                 scale_data = TRUE,
                                 color_col = NULL) {

      if (is.null(self$se_obj)) stop("Please run process_data() first.")
      if (is.null(self$imputed_se)) stop("Please run perform_imputation() first.")

      # 1. Parameter handling
      plot_type <- match.arg(plot_type, several.ok = TRUE)

      # Get data matrices
      raw_mat_full <- assay(self$se_obj)
      imp_mat_full <- assay(self$imputed_se)

      # Get common proteins and samples
      common_proteins <- intersect(rownames(raw_mat_full), rownames(imp_mat_full))
      common_samples <- intersect(colnames(raw_mat_full), colnames(imp_mat_full))

      if (length(common_proteins) == 0 || length(common_samples) == 0) {
        stop("No common proteins or samples found between raw and imputed data.")
      }

      # Use common proteins and samples
      raw_mat <- raw_mat_full[common_proteins, common_samples, drop = FALSE]
      imp_mat <- imp_mat_full[common_proteins, common_samples, drop = FALSE]

      message(paste0("  Using ", length(common_proteins), " proteins and ", length(common_samples), " samples for assessment."))

      # 2. Prepare base data
      n_imputed <- sum(is.na(raw_mat))
      prop_imputed <- round(n_imputed / prod(dim(raw_mat)) * 100, 2)
      message(paste0("--- Imputation Assessment (", prop_imputed, "% values imputed) ---"))

      stats_list <- list()
      plot_list <- list()
      
      # ==========================================================
      # Module A: Structural robustness analysis (Structure / MDS-Procrustes)
      # ==========================================================
      if ("structure" %in% plot_type) {
        message("  1. Calculating Sample Distances (Robustness Check)...")
        
        # A.1 Data preparation and standardization
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
        
        # A.2 MDS computation
        dist_raw <- dist(t(raw_for_dist))       
        dist_imp <- dist(t(imp_for_dist))            
        mds_raw <- cmdscale(dist_raw, k = 2)    
        mds_imp <- cmdscale(dist_imp, k = 2)    
        
        # A.3 Procrustes analysis
        pro_res <- vegan::protest(mds_raw, mds_imp, scores = "sites", permutations = 999)
        m2_val <- round(pro_res$ss, 4)
        
        message(paste0("  -> Procrustes M2: ", m2_val, " (Lower is better)"))
        stats_list$robustness <- list(m2 = m2_val, pval = pro_res$signif, correlation = pro_res$t0)
        
        # A.4 [Key Fix] Force-rename columns to ensure left_join produces correct suffixes
        df_raw <- as.data.frame(pro_res$X)
        colnames(df_raw)[1:2] <- c("Dim1", "Dim2") # Force naming
        df_raw <- df_raw %>% mutate(Type = "Raw") %>% rownames_to_column("Sample")
        
        df_imp <- as.data.frame(pro_res$Yrot)
        colnames(df_imp)[1:2] <- c("Dim1", "Dim2") # Force naming
        df_imp <- df_imp %>% mutate(Type = "Imputed") %>% rownames_to_column("Sample")
        
        # Now Dim1 and Dim2 will conflict, left_join auto-adds suffixes
        plot_df <- left_join(df_raw, df_imp, by = "Sample", suffix = c(".raw", ".imp"))
        
        # Handle color column
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
        
        # Update plot code to use new column names (Dim1.raw, Dim1.imp)
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
      # Module B: Distribution consistency (Density Plot)
      # ==========================================================
      if ("density" %in% plot_type) {
        # Check dimension match; if mismatch, show only imputed data distribution
        if (ncol(raw_mat) != ncol(imp_mat)) {
          # Dimension mismatch, show only current data distribution
          plot_df_dens <- data.frame(
            Intensity = as.vector(imp_mat),
            Type = "Current Data"
          )
          p_dist <- ggplot(plot_df_dens, aes(x = Intensity)) +
            geom_density(fill = "#1465AC", color = "#1465AC", alpha = 0.4) +
            labs(
              title = "Distribution (After Outlier Removal)",
              subtitle = paste0("Samples: ", ncol(imp_mat)),
              x = "Intensity (Log2)", y = "Density"
            ) +
            theme_bw() + theme(plot.title = element_text(face = "bold", hjust = 0.5))
        } else {
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
        }

        plot_list$density <- p_dist
      }

      # ==========================================================
      # Module C: Data integrity check (only when dimensions match)
      # ==========================================================
      if (check_integrity && ncol(raw_mat) == ncol(imp_mat)) {
        mask <- !is.na(raw_mat)
        diff_vals <- abs(raw_mat[mask] - imp_mat[mask])
        if(max(diff_vals) < 1e-9) {
          message("  [Pass] Data Integrity: Original values preserved.")
        } else {
          warning(paste0("  [Alert] Data Integrity: Original values changed! Max diff: ", max(diff_vals)))
        }
      }

      # ==========================================================
      # Final display
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
    
    # --- Feature 3: Plot PCA ---
    plot_pca = function(color_col = NULL) {
      if (is.null(self$imputed_se)) stop("Please run perform_imputation() first.")

      pca_res <- prcomp(t(assay(self$imputed_se)), scale. = TRUE)

      pca_df <- as.data.frame(pca_res$x) %>%
        rownames_to_column("sample") %>%
        left_join(as.data.frame(colData(self$imputed_se)) %>% rownames_to_column("sample"), by="sample")

      color_col <- color_col

      # Get group count to select appropriate color palette
      n_groups <- length(unique(pca_df[[color_col]]))

      p <- ggplot(pca_df, aes(x=PC1, y=PC2, color = .data[[color_col]], label=sample)) +
        geom_point(size=4, alpha=0.8) +
        geom_text_repel(max.overlaps = 10) +
        labs(title = "PCA Analysis (Imputed Data)",
             subtitle = "Check for outliers here") +
        theme_bw() +
        theme(plot.title = element_text(face="bold", hjust=0.5))

      # Use ggsci color palette
      if (n_groups <= 10) {
        p <- p + ggsci::scale_color_npg()
      } else {
        p <- p + ggsci::scale_color_d3(palette = "category20")
      }

      return(p)
    },

    # --- [Added] Feature 4: Automated outlier detection (PCA or correlation) ---
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
        # --- Strategy A: PCA distance method ---
        # Suitable for detecting samples with overall expression deviating from the group
        pca_res <- prcomp(t(mat), scale. = TRUE)
        coords <- as.data.frame(pca_res$x[, 1:2]) # Take first two axes
        
        # Calculate Z-scores for PC1 and PC2
        coords$z_pc1 <- (coords$PC1 - mean(coords$PC1)) / sd(coords$PC1)
        coords$z_pc2 <- (coords$PC2 - mean(coords$PC2)) / sd(coords$PC2)
        
        # Define outlier: any axis exceeds threshold
        is_outlier <- abs(coords$z_pc1) > sd_threshold | abs(coords$z_pc2) > sd_threshold
        outlier_samples <- rownames(coords)[is_outlier]
        
        stats_df$Score1 <- coords$PC1
        stats_df$Score2 <- coords$PC2
        stats_df$Is_Outlier <- is_outlier
        
        # Plot
        p <- ggplot(stats_df, aes(x = Score1, y = Score2, color = Is_Outlier, label = Sample)) +
          geom_point(size = 3, alpha = 0.8) +
          geom_text_repel(data = subset(stats_df, Is_Outlier), max.overlaps = 20, size = 3) +
          stat_ellipse(level = 0.99, linetype = "dashed", color = "grey50") + # 99% (approx 2.5-3 SD range)
          scale_color_manual(values = c("FALSE"="#1465AC", "TRUE"="#B31B21")) +
          labs(title = "Outlier Detection (PCA)", 
               x = "PC1", y = "PC2", 
               subtitle = paste("Outliers:", length(outlier_samples))) +
          theme_bw()
        
      } else {
        # --- Strategy B: Average correlation method (Connectivity) ---
        # Suitable for detecting samples inconsistent with the group (e.g., failed experiments)
        
        # Calculate correlation matrix
        cor_mat <- cor(mat, method = "pearson")
        
        # Calculate each sample's average correlation with others (Connectivity)
        # diag=NA prevents self-correlation of 1 inflating the mean
        diag(cor_mat) <- NA
        mean_cors <- colMeans(cor_mat, na.rm = TRUE)
        
        # Calculate Z-score (lower is worse)
        z_scores <- (mean_cors - mean(mean_cors)) / sd(mean_cors)
        
        # Define outlier: significantly low correlation (one-tailed test, negative direction only)
        is_outlier <- z_scores < -sd_threshold
        outlier_samples <- names(mean_cors)[is_outlier]
        
        stats_df$Mean_Cor <- mean_cors
        stats_df$Z_Score <- z_scores
        stats_df$Is_Outlier <- is_outlier
        
        # Plot
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
      
      # Execute removal
      if (length(outlier_samples) > 0) {
        message(paste0("  Identified ", length(outlier_samples), " outliers: ", paste(outlier_samples, collapse = ", ")))
        if (remove) {
          self$remove_outliers(outlier_samples) # Call existing removal function
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
        # 1. Update imputed_se object
        self$imputed_se <- self$imputed_se[, valid_keep]

        # [Added] 2. Also update se_obj for consistency
        if (!is.null(self$se_obj)) {
          self$se_obj <- self$se_obj[, valid_keep]
        }

        # 3. Update meta_data (using base R for safety, matching by rownames)
        self$meta_data <- self$meta_data[valid_keep, , drop = FALSE]

        message(paste("Removed Outliers:", paste(outlier_samples_list, collapse=", ")))
        message(paste("Remaining samples:", length(valid_keep)))
      } else {
        message("No matching outliers found to remove.")
      }
    },

    # --- [Added] Feature: Filter samples by metadata (for downstream analysis) ---
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

      
      # Directly modify imputed_se and meta_data
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
      
      # Restore current working object from backup
      self$imputed_se <- self$imputed_se_backup
      self$meta_data <- self$meta_data_backup
      
      # Also try to restore se_obj (analysis mainly uses imputed_se, but consistency is better)
      # Note: If se_obj has no structural changes before/after imputation, restoring is optional,
      # but for safety, we usually restore SE to full sample state (even with pre-imputation values)
      # For simplicity, we mainly ensure the imputed_se used for downstream analysis is complete
      
      message(paste("  Data restored. Current samples:", ncol(self$imputed_se)))
    },

    # --- [Added] Feature: Auto-plot protein expression trends (auto-detect continuous vs categorical) ---
    plot_protein_expression = function(proteins, variable, color_var = NULL, add_labels = FALSE) {
      if (is.null(self$imputed_se)) stop("Data not ready. Please run perform_imputation() first.")
      
      # 1. Check if variable exists
      meta <- as.data.frame(colData(self$imputed_se))
      if (!variable %in% colnames(meta)) {
        stop(paste("Variable '", variable, "' not found in metadata.", sep=""))
      }
      
      # 2. Check if proteins exist
      valid_prots <- intersect(proteins, rownames(self$imputed_se))
      if (length(valid_prots) == 0) stop("None of the specified proteins found in the dataset.")
      if (length(valid_prots) < length(proteins)) {
        warning(paste("Some proteins were not found:", paste(setdiff(proteins, valid_prots), collapse=", ")))
      }
      
      # 3. Prepare plot data (Long Format)
      # Use drop=FALSE to prevent matrix-to-vector conversion with single protein
      expr_mat <- assay(self$imputed_se)[valid_prots, , drop = FALSE]
      
      plot_df <- t(expr_mat) %>%
        as.data.frame() %>%
        rownames_to_column("sample_id") %>%
        left_join(meta %>% rownames_to_column("sample_id"), by = "sample_id") %>%
        pivot_longer(cols = all_of(valid_prots), names_to = "Protein", values_to = "Expression")
      
      # 4. Determine variable type and plot
      # Default color variable matches x-axis variable, unless user specifies
      if (is.null(color_var)) color_var <- variable
      
      x_vals <- plot_df[[variable]]
      
      # Initialize ggplot
      p <- ggplot(plot_df, aes(x = .data[[variable]], y = Expression))

      # Get group count
      n_groups <- length(unique(plot_df[[color_var]]))

      # --- Branch A: Continuous variable (Numeric) ---
      if (is.numeric(x_vals) && length(unique(x_vals)) > 2) {
        message(paste("Detected numeric variable '", variable, "'. Using Scatter plot + Loess.", sep=""))

        p <- p +
          geom_point(aes(color = .data[[color_var]]), size = 2.5, alpha = 0.7) +
          geom_smooth(method = "loess", color = "#B31B21", fill = "#B31B21", alpha = 0.2, linewidth = 1) +
          scale_color_viridis_c(option = "D", begin = 0.2, end = 0.8) + # Default continuous palette
          theme_bw()

      } else {
        # --- Branch B: Categorical variable (Factor/Character) ---
        message(paste("Detected categorical variable '", variable, "'. Using Violin + Boxplot.", sep=""))

        # Force convert to factor for discrete plotting
        plot_df[[variable]] <- as.factor(plot_df[[variable]])

        p <- p +
          geom_violin(aes(fill = .data[[color_var]]), alpha = 0.5, trim = FALSE, scale = "width") +
          geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.9) +
          geom_jitter(width = 0.1, size = 1.5, alpha = 0.6, color = "grey30") +
          theme_bw()

        # Use ggsci color palette
        if (n_groups <= 10) {
          p <- p + ggsci::scale_fill_npg()
        } else {
          p <- p + ggsci::scale_fill_d3(palette = "category20")
        }
      }
      
      # 5. Common formatting
      # If multiple proteins, use faceted display
      if (length(valid_prots) > 1) {
        p <- p + facet_wrap(~Protein, scales = "free_y")
      } else {
        p <- p + labs(title = valid_prots[1])
      }
      
      # Add sample labels (optional)
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

# --- Class 2: Differential analysis (Factory pattern: supports DEP and Regression) ---
DiffExpAnalyst <- R6::R6Class("DiffExpAnalyst",
  public = list(
    se_obj = NULL,
    diff_results = NULL, # Store final differential results table
    sig_results = NULL,
    
    initialize = function(se_object) {
      self$se_obj <- se_object
    },
    
    # Mode A: Group comparison (using DEP workflow)
    # --- Core modification: fully emulate DEP test_diff behavior, with covariate support ---
    run_dep_analysis = function(condition_col = "condition", control_group, case_group, covariates = NULL, paired_col = NULL) {
      
      message(paste0("--- Running DEP-style Analysis with Covariates: ", case_group, " vs ", control_group, " ---"))
      
      # 1. Prepare data
      data_se <- self$se_obj
      meta_data <- colData(data_se)
      meta_data[[condition_col]] = as.factor(as.character(meta_data[[condition_col]]))
      
      # 2. Build design matrix (with covariates)
      # ---------------------------------------------------------
      # Ensure group variable is a factor
      if(!is.factor(meta_data[[condition_col]])) {
         meta_data[[condition_col]] <- factor(meta_data[[condition_col]])
      }
      
      # Base formula: ~ 0 + Group 
      formula_str <- paste0("~ 0 + ", condition_col)
      
      # Add covariates
      if (!is.null(covariates)) {
        # Check if covariates exist
        missing_cov <- setdiff(covariates, colnames(meta_data))
        if (length(missing_cov) > 0) stop(paste("Covariates not found:", paste(missing_cov, collapse=", ")))
        
        formula_str <- paste0(formula_str, " + ", paste(covariates, collapse = " + "))
        message(paste("  Adjusting for:", paste(covariates, collapse=", ")))
      }
      
      # Generate design matrix
      design <- model.matrix(as.formula(formula_str), data = meta_data)
      
      # 3. Clean column names (for Limma Contrast compatibility)
      # ---------------------------------------------------------
      # model.matrix generates names like "conditionD36", we need to clean back to "D36"
      group_levels <- levels(meta_data[[condition_col]])
      current_cols <- colnames(design)
      
      for(lvl in group_levels) {
        # Regex match: ^condition_col + level$ (e.g. ^GroupD36$)
        pattern <- paste0("^", condition_col, lvl, "$")
        idx <- grep(pattern, current_cols)
        if(length(idx) > 0) current_cols[idx] <- lvl
      }
      colnames(design) <- current_cols
      
      # 4. Run Limma analysis
      # ---------------------------------------------------------
      contrast_formula <- paste0(case_group, " - ", control_group)
      contrast_name    <- paste0(case_group, "_vs_", control_group) # DEP required naming format


      
      message(paste("  Contrast:", contrast_formula))

      if (!is.null(paired_col)) {
        if (!paired_col %in% colnames(meta_data)) stop(paste("Paired column not found:", paired_col))
        block_var <- as.factor(meta_data[[paired_col]])
        message(paste("  Paired/blocked on:", paired_col))
        corfit <- duplicateCorrelation(assay(data_se), design, block = block_var)
        message(paste("  Consensus correlation:", round(corfit$consensus, 4)))
        fit <- lmFit(assay(data_se), design, block = block_var, correlation = corfit$consensus)
      } else {
        fit <- lmFit(assay(data_se), design)
      }
      cont.matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
      fit2 <- contrasts.fit(fit, cont.matrix)
      fit2 <- eBayes(fit2)
      
      # Get statistical results
      res_table <- topTable(fit2, number = Inf, sort.by = "none") # Must use sort.by="none" to maintain order consistent with se
      res_table$adj.P.Val = p.adjust(res_table$P.Value, method = "fdr")

      # 5. Key step: construct rowData structure conforming to DEP requirements
      # ---------------------------------------------------------
      # DEP test_diff adds three columns to rowData:
      # {contrast}_diff   (i.e. logFC)
      # {contrast}_p.val  (i.e. P.Value)
      # {contrast}_p.adj  (i.e. adj.P.Val)
      
      # Extract current rowData
      row_data <- rowData(data_se)
      
      # Prepare new columns (column names must strictly follow DEP specification)
      col_diff <- paste0(contrast_name, "_diff")
      col_pval <- paste0(contrast_name, "_p.val")
      col_padj <- paste0(contrast_name, "_p.adj")
      
      # Inject Limma results into row_data
      # Ensure row order consistency (topTable sort.by='none' and se unchanged, usually consistent, but use match for safety)
      mm <- match(rownames(row_data), rownames(res_table))
      
      row_data[[col_diff]] <- res_table$logFC[mm]
      row_data[[col_pval]] <- res_table$P.Value[mm]
      row_data[[col_padj]] <- res_table$adj.P.Val[mm]
      
      # 6. Update and return SE object
      # ---------------------------------------------------------
      rowData(data_se) <- row_data
      
      # Store a table version internally (optional, for easy inspection)
      self$diff_results <- res_table %>% rownames_to_column("Protein")
      
      message("  Analysis complete. Results added to SummarizedExperiment rowData.")
      message(paste("  Columns added:", col_diff, ",", col_pval, ",", col_padj))
      
      return(data_se) # Return SE object for direct use with add_rejections()
    },
    
    # Mode B: Continuous variable regression (Linear or Spline)
    run_continuous_analysis = function(time_col = "time_num", method = "spline", df = 3) {
      message(paste("Running", method, "regression on", time_col))
      
      # Prepare long-format data
      expr_long <- assay(self$se_obj) %>%
        as.data.frame() %>%
        rownames_to_column("Protein") %>%
        pivot_longer(-Protein, names_to = "sample", values_to = "value") %>%
        left_join(as.data.frame(colData(self$se_obj)) %>% rownames_to_column("sample"), by = "sample")
      
      # Ensure time column is numeric
      current_vals = expr_long[[time_col]]
      if (!is.numeric(current_vals)) {
        message(paste0("Notice: Column '", time_col, "' is not numeric (Type: ", class(current_vals)[1], "). Attempting to convert..."))
        
        # 1. Convert to character first
        vals_char <- as.character(current_vals)
        
        # 2. Use readr::parse_number to intelligently extract numbers
        vals_num <- readr::parse_number(vals_char)
        
        # 3. Check conversion result
        if (all(is.na(vals_num))) {
          stop(paste0("Error: Failed to convert column '", time_col, "' to numeric. ",
                      "It implies that the labels do not contain extractable numbers (e.g., 'Control', 'Treat')."))
        }
        
        # 4. Overwrite original column
        expr_long[[time_col]] <- vals_num
        message(paste0("  Conversion successful. Example: '", vals_char[1], "' -> ", vals_num[1]))
      }
      
      # Parallel computation (recommended, as per-protein regression is slow)
      results <- expr_long %>%
        group_by(Protein) %>%
        do({
          dat <- .
          if (method == "spline") {
            # Spline Regression (logic from Helper 2)
            fit_full <- lm(value ~ ns(get(time_col), df = df), data = dat)
            fit_null <- lm(value ~ 1, data = dat)
          } else {
            # Linear Regression
            fit_full <- lm(value ~ get(time_col), data = dat)
            fit_null <- lm(value ~ 1, data = dat)
          }
          
          lrt_res <- anova(fit_null, fit_full)
          model_stats <- glance(fit_full)
          
          # Calculate correlation direction
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

    # Mode C: Between-group ANOVA analysis (multi-group comparison)
    run_anova_analysis = function(condition_col = "condition", covariates = NULL, selected_groups = NULL) {
      message(paste("Running ANOVA analysis on", condition_col))

      # 1. Prepare data
      data_se <- self$se_obj
      meta_data <- colData(data_se)

      # Ensure group variable is a factor
      if (!is.factor(meta_data[[condition_col]])) {
        meta_data[[condition_col]] <- as.factor(as.character(meta_data[[condition_col]]))
      }
      meta_data[[condition_col]] <- droplevels(meta_data[[condition_col]])

      # 2. Filter specific groups (if specified)
      all_group_levels <- levels(meta_data[[condition_col]])
      if (!is.null(selected_groups)) {
        invalid_groups <- setdiff(selected_groups, all_group_levels)
        if (length(invalid_groups) > 0) {
          stop(paste("Selected groups not found:", paste(invalid_groups, collapse = ", ")))
        }
        message(paste("  Filtering to selected groups:", paste(selected_groups, collapse = ", ")))
        
        keep_mask <- meta_data[[condition_col]] %in% selected_groups
        data_se <- data_se[, keep_mask]
        meta_data <- colData(data_se) # Update meta_data
        meta_data[[condition_col]] <- droplevels(meta_data[[condition_col]])
      }
      
      # Check group count
      group_levels <- levels(meta_data[[condition_col]])
      n_groups <- length(group_levels)
      if (n_groups < 2) stop("ANOVA requires at least 2 groups.")
      message(paste("  Analyzing", n_groups, "groups using F-test."))

      # 3. Build design matrix (intercept model ~ Condition)
      # [Fix 1] Remove '0 +', use default intercept model
      # This way coef=1 is intercept(Reference), coef=2..N are group differences relative to Reference
      formula_str <- paste0("~ ", condition_col)

      if (!is.null(covariates)) {
        missing_cov <- setdiff(covariates, colnames(meta_data))
        if (length(missing_cov) > 0) stop(paste("Covariates not found:", paste(missing_cov, collapse=", ")))
        formula_str <- paste0(formula_str, " + ", paste(covariates, collapse = " + "))
        message(paste("  Adjusting for covariates:", paste(covariates, collapse=", ")))
      }

      design <- model.matrix(as.formula(formula_str), data = meta_data)

      # 4. Run Limma
      fit <- lmFit(assay(data_se), design)
      
      fit2 <- eBayes(fit)

      # 5. Extract ANOVA results (F-test)
      # [Fix 3] Only test group-related columns excluding Intercept (column 1)
      # Typically condition-related coefficients are columns 2 through n_groups
      # Note: If covariates exist, their columns follow condition; we usually only test condition differences

      # Smart identification of coefficients to test (columns containing condition_col, excluding intercept)
      coefs_to_test <- grep(condition_col, colnames(design))

      # If grep finds nothing (e.g., column names were escaped), fallback to testing 2:n_groups
      if(length(coefs_to_test) == 0) {
         coefs_to_test <- 2:n_groups
      }

      # [Key Fix] Check again if coefs_to_test is empty (possibly n_groups < 2)
      if(length(coefs_to_test) == 0) {
        stop(paste("Cannot run ANOVA: Unable to identify coefficients for condition column '", condition_col, "'. Check if the condition column has at least 2 unique groups.", sep=""))
      }

      message(paste("  Testing coefficients columns:", paste(coefs_to_test, collapse=", ")))
      
      # Sort selected coefficients
      anova_res <- topTable(fit2, coef = coefs_to_test, number = Inf, sort.by = "p")

      # Organize results
      res_table <- as.data.frame(anova_res) %>%
        rownames_to_column("Protein") %>%
        mutate(
          anova_F = F,
          anova_pval = P.Value # This P.Value is the F-test p-value
        ) %>%
        dplyr::select(Protein, anova_F, P.Value, adj.P.Val, anova_pval)

      # 6. Store results
      self$diff_results <- res_table

      # Simple result diagnostics
      sig_count <- sum(res_table$adj.P.Val < 0.05, na.rm = TRUE)
      message("  ANOVA analysis complete.")
      message(paste("  Significant proteins (adj.P.Val < 0.05):", sig_count, "/", nrow(res_table)))
      
      # If still all significant, print warning
      if(sig_count > 0.8 * nrow(res_table)) {
        warning("  [Warning] >80% proteins are significant. Check if normalization was performed (e.g., log transformation).")
      }

      return(res_table)
    },
    # --- 2. Standalone plot function: ANOVA differential heatmap ---

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
        # If too many proteins, heatmap becomes unreadable; auto-downsample as protection
        mat_full <- assay(se_obj)[sig_proteins, ]
        vars <- apply(mat_full, 1, var)
        sig_proteins <- names(sort(vars, decreasing = TRUE))[1:2000]
      }
      
      # 1. Extract data
      mat <- assay(se_obj)[sig_proteins, ]
      
      # 2. Z-score standardization (Row scaling)
      if(scale_rows) {
        mat <- t(scale(t(mat)))
      }
      
      # 3. Build column annotation (sample info)
      meta <- as.data.frame(colData(se_obj))
      
      # Prepare annotation data frame
      anno_df <- meta %>% dplyr::select(all_of(group_col))
      if(!is.null(time_col) && time_col %in% colnames(meta)) {
        anno_df[[time_col]] <- meta[[time_col]]
      }
      
      # Auto-generate colors
      # Generate discrete colors for Group
      groups <- unique(anno_df[[group_col]])
      group_colors <- setNames(ggsci::pal_npg()(length(groups)), groups)
      
      anno_colors <- list()
      anno_colors[[group_col]] <- group_colors
      
      # Create ComplexHeatmap annotation object
      col_anno <- HeatmapAnnotation(
        df = anno_df,
        col = anno_colors,
        show_annotation_name = TRUE
      )
      
      # 4. Draw heatmap
      ht <- Heatmap(
        mat,
        name = "Z-score",
        top_annotation = col_anno,

        # Row and column clustering settings
        cluster_rows = TRUE,           # Cluster proteins
        cluster_columns = TRUE,        # Cluster samples

        # [Fixed] Add clustering distance method and dendrogram settings
        clustering_distance_columns = "euclidean",  # Column clustering distance method
        clustering_method_columns = "ward.D2",      # Column clustering method

        # Show column dendrogram
        show_column_dend = TRUE,
        column_dend_height = unit(10, "mm"),

        # Visual adjustments
        show_row_names = length(sig_proteins) < 50, # Show names when few proteins, hide when many
        show_column_names = TRUE,

        # Colors (blue-white-red)
        col = colorRamp2(c(-2, 0, 2), c("#3C5488B2", "white", "#E64B35B2")),

        column_title = title,

        # Auto-split rows (K-means): automatically group proteins into 2-4 pattern clusters
        row_km = ifelse(length(sig_proteins) > 100, 4, 1)
      )
      
      return(ht)
    },

    # Get significant protein list
    get_sig_proteins = function(pval_cutoff = 0.05,
                                    r2_cutoff = 0.5,
                                    corr_cutoff = 0.5,
                                    logfc_cutoff = 1,
                                    use_adjusted_pval = TRUE) {

          if (is.null(self$diff_results)) stop("No analysis results found. Run analysis first.")

          # Mode 1: ANOVA analysis (check anova_pval column)
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

          # Mode 2: Continuous regression (check adj_r_squared column)
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

          # Mode 3: Two-group comparison (DEP/Limma)
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
    # Auto-detect Group comparison vs Continuous regression based on diff_results content
    plot_volcano = function(logfc_cutoff = 1,       # Group mode only
                            corr_cutoff = 0.5,      # Regression mode only
                            r2_cutoff = 0.5,        # Regression mode only
                            pval_cutoff = 0.05,
                            top_n_labels = 10,
                            use_adjusted_pval = TRUE) {
      if(is.null(self$diff_results)) stop("No analysis results found. Please run run_dep_analysis() or run_continuous_analysis() first.")

      df <- self$diff_results

      # --- 0. Check if ANOVA result (ANOVA has no volcano plot) ---
      if ("anova_F" %in% colnames(df)) {
        stop("Volcano plot is not available for ANOVA analysis. ANOVA compares multiple groups simultaneously and does not produce pairwise log fold changes. Please use the Heatmap or other visualization options.")
      }

      # --- 1. Determine analysis type and prepare plot data ---

      if ("adj_r_squared" %in% colnames(df)) {
        # >>> Mode A: Continuous variable regression (Regression) <<<
        message("Detected Continuous Regression results. Plotting: X=Spearman Rho, Y=-log10(Pval)")

        # Determine P-value column
        p_col <- if(use_adjusted_pval) "adj_pval" else "p_value"

        # Prepare data
        plot_data <- df %>%
          dplyr::select(Protein,
                        x_val = spearman_rho,      # X-axis: correlation coefficient (direction)
                        y_val = !!sym(p_col),      # Y-axis: P-value
                        r2 = adj_r_squared) %>%    # Auxiliary: R2 for significance filtering
          filter(!is.na(x_val) & !is.na(y_val)) %>%
          mutate(
            log_pval = -log10(y_val),
            # Logic: P-value significant AND |Rho| meets cutoff AND R2 meets cutoff, direction by Rho sign
            expression = case_when(
              y_val < pval_cutoff & abs(x_val) > corr_cutoff & r2 > r2_cutoff & x_val > 0 ~ "Up-regulated",
              y_val < pval_cutoff & abs(x_val) > corr_cutoff & r2 > r2_cutoff & x_val < 0 ~ "Down-regulated",
              TRUE ~ "Not Significant"
            )
          )

        # Set axis labels and threshold line positions
        x_lab <- "Spearman Correlation (rho)"
        y_lab <- if(use_adjusted_pval) expression(-log[10] ~ "Adjusted P-value") else expression(-log[10] ~ "P-value")
        v_lines <- c(-corr_cutoff, corr_cutoff) # Use corr_cutoff as vertical line threshold
        title_suffix <- "(Regression)"

      } else {
        # >>> Mode B: Group differential analysis (DEP/Limma) <<<
        message("Detected Group Comparison results. Plotting: X=LogFC, Y=-log10(Pval)")
        
        # Determine P-value column (run_dep_analysis columns: P.Value and adj.P.Val)
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
      
      # --- 2. Extract top genes (for labeling) ---
      top_up <- plot_data %>%
        filter(expression == "Up-regulated") %>%
        arrange(desc(x_val)) %>% # Whether LogFC or Rho, larger positive = more significant
        slice_head(n = top_n_labels)
      
      top_down <- plot_data %>%
        filter(expression == "Down-regulated") %>%
        arrange(x_val) %>%       # More negative = more significant
        slice_head(n = top_n_labels)
      
      top_genes <- bind_rows(top_up, top_down)
      
      # --- 3. Draw volcano plot ---
      p <- ggplot(plot_data, aes(x = x_val, y = log_pval)) +
        # Scatter points
        geom_point(aes(color = expression), alpha = 0.6, size = 1.5) +
        
        # Color definitions
        scale_color_manual(values = c("Up-regulated" = "#B31B21", 
                                      "Down-regulated" = "#1465AC", 
                                      "Not Significant" = "grey80")) +
        
        # Threshold lines
        geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "black", alpha = 0.5) +
        geom_vline(xintercept = v_lines, linetype = "dashed", color = "black", alpha = 0.5) +
        
        # Labels
        geom_text_repel(data = top_genes, 
                        aes(label = Protein),
                        size = 3,
                        box.padding = 0.5,
                        max.overlaps = Inf, 
                        segment.color = "grey50") +
        
        # Title and axes
        labs(title = paste("Volcano Plot", title_suffix),
             subtitle = paste("Up:", nrow(filter(plot_data, expression == "Up-regulated")), 
                              "| Down:", nrow(filter(plot_data, expression == "Down-regulated"))),
             x = x_lab,
             y = y_lab) +
        
        # Theme
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

# --- Class 3: Mfuzz clustering wrapper ---
MfuzzClusterer <- R6::R6Class("MfuzzClusterer",
  public = list(
    se_obj = NULL,
    mfuzz_set = NULL,
    cl_obj = NULL,
    enrich_worker = NULL,
    
    initialize = function(se_object, cache_manager = NULL) {
      self$se_obj <- se_object
      # Composition pattern: MfuzzClusterer "has-a" EnrichmentAnalyst
      self$enrich_worker <- EnrichmentAnalyst$new(cache_manager)
    },
    
    run_mfuzz = function(sig_prots, time_col = "time_num", centers = 4) {
      
      # 1. Prepare mean data (average by time point)
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
      
      # 2. Build ExpressionSet and standardize
      self$mfuzz_set <- ExpressionSet(assayData = df_mean)
      self$mfuzz_set <- standardise(self$mfuzz_set)
      
      # 3. Clustering
      m <- mestimate(self$mfuzz_set)
      set.seed(123)
      self$cl_obj <- mfuzz(self$mfuzz_set, c = centers, m = m)
      message("Mfuzz clustering complete.")
    },


    # --- Core feature: Auto-complete enrichment and plotting ---
    plot_clusters = function(mfrow = c(2,2)) {
    mfuzz.plot2(self$mfuzz_set, cl = self$cl_obj, min.mem = 0.7, mfrow = mfrow, 
                time.labels = colnames(exprs(self$mfuzz_set)), colo = "fancy", x11=FALSE)
    },
    
    get_cluster_df = function(min.acore_threshold = 0.7) {
      acore(self$mfuzz_set, self$cl_obj, min.acore = min.acore_threshold) %>% 
        bind_rows(.id="Cluster") %>% 
        rename(Gene = NAME)
    },

    # --- Plot logic: Trend plot + ORA dotplot ---
    draw_composite_plot = function(cluster_id, gene_list,  time_col = time_col, ora_res, db = c("GO_BP", "KEGG")) {
      # 1. Trend plot (Loess)      
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
      
      # 2. db1 dotplot
      p_db1 <- if(!is.null(ora_res[[db[1]]]) && nrow(ora_res[[db[1]]]) > 0) {
        clusterProfiler::dotplot(ora_res[[db[1]]], showCategory = 10, color = "pvalue") + ggtitle(db[1]) + theme(axis.text.y = element_text(size = 8))
      } else { ggplot() + theme_void() + annotate("text", x=0.5, y=0.5, label="No Sig pathways") }
      
      # 3. db2 dotplot
      p_db2 <- if(!is.null(ora_res[[db[2]]]) && nrow(ora_res[[db[2]]]) > 0) {
        clusterProfiler::dotplot(ora_res[[db[2]]], showCategory = 10, color = "pvalue") + ggtitle(db[2]) + theme(axis.text.y = element_text(size = 8))
      } else { ggplot() + theme_void() + annotate("text", x=0.5, y=0.5, label="No Sig pathways") }
      
      # 4. Compose plots
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
      
      # 1. Get clustering info
      cluster_df <- self$get_cluster_df(min.acore_threshold = min.acore_threshold)
      
      clusters <- unique(cluster_df$Cluster)
      excel_list <- list()
      plot_list <- list()
      
      # 2. Loop through each cluster
        if (!dir.exists(target_dir)) {
          dir.create(target_dir)
        }

      for(cl in clusters) {
        message(paste("  Processing Cluster", cl, "..."))
        target_genes <- cluster_df %>% filter(Cluster == cl) %>% pull(Gene)
        
        # [Key] Call internal worker for analysis
        universe_genes = rownames(self$se_obj)
        ora_res <- self$enrich_worker$run_comprehensive_ora(target_genes, universe_genes, pval_cutoff = pvalueCutoff)
        
        # Save tables

        self$enrich_worker$enrich_to_excel(ora_res, direction = "Mfuzz", output_prefix = paste0("Mfuzz_cluster", "_", cl), target_dir = target_dir)
        
        # Draw composite plot
        p_combined <- self$draw_composite_plot(cl, target_genes, time_col, ora_res, db = view_db)
        plot_list[[as.character(cl)]] <- p_combined
      }

        # Save plots
        pdf(here::here(target_dir, "Mfuzz_Cluster_Analysis_Summary.pdf"), width = 14, height = 10)
        for(p in plot_list) { print(p) }
        dev.off()

        return(plot_list)
      }

  )
)


# =============================================================================
# Utility: Convert fgsea result to clusterProfiler gseaResult S4 object
# Enables enrichplot::dotplot() and enrichplot::gseaplot2() compatibility
# =============================================================================
fgsea_to_gseaResult <- function(fgsea_res, gene_list, gene_sets, term2name = NULL) {
  if (nrow(fgsea_res) == 0) return(NULL)

  # Build description mapping
  if (!is.null(term2name) && nrow(term2name) > 0) {
    desc_map <- stats::setNames(term2name$name, term2name$term)
  } else {
    desc_map <- stats::setNames(fgsea_res$pathway, fgsea_res$pathway)
  }

  # Build result data.frame matching gseaResult@result format
  result_df <- data.frame(
    ID              = fgsea_res$pathway,
    Description     = ifelse(fgsea_res$pathway %in% names(desc_map),
                             desc_map[fgsea_res$pathway],
                             fgsea_res$pathway),
    setSize         = fgsea_res$size,
    enrichmentScore = fgsea_res$ES,
    NES             = fgsea_res$NES,
    pvalue          = fgsea_res$pval,
    p.adjust        = stats::p.adjust(fgsea_res$pval, method = "BH"),
    qvalue          = stats::p.adjust(fgsea_res$pval, method = "BH"),
    rank            = vapply(seq_len(nrow(fgsea_res)), function(i) {
      pw <- fgsea_res$pathway[i]
      genes_in_set <- gene_sets[[pw]]
      ranked_positions <- match(genes_in_set, names(gene_list))
      ranked_positions <- ranked_positions[!is.na(ranked_positions)]
      if (length(ranked_positions) == 0) return(NA_real_)
      as.numeric(ranked_positions[which.max(abs(cumsum(
        ifelse(seq_along(gene_list) %in% ranked_positions, 1, 0) -
        length(ranked_positions) / length(gene_list)
      )[ranked_positions]))])
    }, numeric(1)),
    leading_edge    = vapply(fgsea_res$leadingEdge, function(le) {
      paste0("tags=", length(le), "%")
    }, character(1)),
    core_enrichment = vapply(fgsea_res$leadingEdge, function(le) {
      paste(le, collapse = "/")
    }, character(1)),
    stringsAsFactors = FALSE
  )

  # Add sign column for dotplot faceting
  result_df$.sign <- ifelse(result_df$NES > 0, "activated", "suppressed")
  rownames(result_df) <- result_df$ID

  # Construct gseaResult S4 object
  params_list <- list(
    pvalueCutoff  = 1,
    nPerm         = 0,
    pAdjustMethod = "BH",
    exponent      = 1,
    minGSSize     = 10,
    maxGSSize     = 500
  )

  new("gseaResult",
    result     = result_df,
    geneSets   = gene_sets,
    geneList   = gene_list,
    params     = params_list,
    readable   = FALSE
  )
}

# --- Class 4: Enrichment analysis wrapper (Cache-based: enricher / fgsea) ---
EnrichmentAnalyst <- R6::R6Class("EnrichmentAnalyst",
  public = list(
    cache_manager = NULL,
    ora_up = NULL,
    ora_down = NULL,
    gsea_res = NULL,

    initialize = function(cache_manager = NULL) {
      if (is.null(cache_manager)) {
        cache_manager <- GeneSetCacheManager$new()
      }
      self$cache_manager <- cache_manager
    },

    # Comprehensive ORA using cached TERM2GENE + enricher()
    run_comprehensive_ora = function(gene_list, universe_genes, pval_cutoff,
                                     padjust_method = "none", qval_cutoff = 1) {
      dbs <- c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome", "Wiki")
      res_list <- list()

      for (db in dbs) {
        tryCatch({
          cached <- self$cache_manager$get_term2gene(db)
          res_list[[db]] <- clusterProfiler::enricher(
            gene         = gene_list,
            universe     = universe_genes,
            TERM2GENE    = cached$TERM2GENE,
            TERM2NAME    = cached$TERM2NAME,
            pvalueCutoff = pval_cutoff,
            pAdjustMethod = padjust_method,
            qvalueCutoff = qval_cutoff,
            minGSSize    = 10,
            maxGSSize    = 500
          )
        }, error = function(e) {
          message("  ORA failed for ", db, ": ", e$message)
          res_list[[db]] <<- NULL
        })
      }

      return(res_list)
    },

    # Comprehensive GSEA using fgsea + cached TERM2GENE (parallel across databases)
    run_comprehensive_gsea = function(ranked_gene_list, pval_cutoff,
                                      padjust_method = "BH",
                                      n_threads = 4) {
      dbs <- c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome", "Wiki")

      # 1. Pre-load all caches and convert to fgsea-compatible format
      all_caches <- list()
      all_gene_sets <- list()
      expressed_genes <- names(ranked_gene_list)

      for (db in dbs) {
        tryCatch({
          cached <- self$cache_manager$get_term2gene(db)
          t2g <- cached$TERM2GENE
          # Filter to expressed genes only
          t2g <- t2g[t2g$gene %in% expressed_genes, ]
          gs <- split(t2g$gene, t2g$term)
          gs <- gs[lengths(gs) >= 10 & lengths(gs) <= 500]
          all_caches[[db]] <- cached
          all_gene_sets[[db]] <- gs
        }, error = function(e) {
          message("  Failed to load cache for ", db, ": ", e$message)
        })
      }

      valid_dbs <- names(all_gene_sets)
      if (length(valid_dbs) == 0) {
        warning("No gene set databases loaded successfully.")
        return(list())
      }

      # 2. Run fgsea in parallel across databases
      run_single_db <- function(db) {
        gs <- all_gene_sets[[db]]
        if (is.null(gs) || length(gs) == 0) return(NULL)

        tryCatch({
          fgsea_res <- fgsea::fgseaMultilevel(
            pathways  = gs,
            stats     = ranked_gene_list,
            minSize   = 10,
            maxSize   = 500,
            nPermSimple = 10000
          )

          term2name <- all_caches[[db]]$TERM2NAME
          fgsea_to_gseaResult(fgsea_res, ranked_gene_list, gs, term2name)
        }, error = function(e) {
          message("  GSEA (fgsea) failed for ", db, ": ", e$message)
          NULL
        })
      }

      if (n_threads > 1 && length(valid_dbs) > 1) {
        message("  Running fgsea in parallel across ", length(valid_dbs),
                " databases with ", n_threads, " threads...")

        if (.Platform$OS.type == "windows") {
          cl <- parallel::makeCluster(min(n_threads, length(valid_dbs)))
          on.exit(parallel::stopCluster(cl), add = TRUE)
          parallel::clusterExport(cl, c("all_gene_sets", "all_caches",
                                         "ranked_gene_list", "fgsea_to_gseaResult"),
                                 envir = environment())
          parallel::clusterEvalQ(cl, {
            requireNamespace("fgsea", quietly = TRUE)
          })
          results <- parallel::parLapply(cl, valid_dbs, run_single_db)
        } else {
          results <- parallel::mclapply(valid_dbs, run_single_db,
                                         mc.cores = min(n_threads, length(valid_dbs)))
        }
        names(results) <- valid_dbs
      } else {
        message("  Running fgsea sequentially across ", length(valid_dbs), " databases...")
        results <- lapply(valid_dbs, run_single_db)
        names(results) <- valid_dbs
      }

      # Remove NULLs
      results <- results[!vapply(results, is.null, logical(1))]
      return(results)
    },

    # --- Business logic: Process DiffExpAnalyst object ---
    analyze_diff_obj = function(diff_obj, pval_cutoff = 1) {
      if(is.null(diff_obj$sig_results) || is.null(diff_obj$diff_results)) stop("Run DE analysis first.")

      message("=== EnrichmentAnalyst: Processing DiffExp Object ===")

      # 1. Extract data
      universe <- diff_obj$diff_results$Protein
      sig_df <- diff_obj$sig_results

      # 2. Separate up and down regulated
      if ("logFC" %in% colnames(sig_df)) {
        up_genes <- sig_df %>% filter(logFC > 0) %>% pull(Protein)
        down_genes <- sig_df %>% filter(logFC < 0) %>% pull(Protein)
        rank_vec <- setNames(diff_obj$diff_results$logFC, diff_obj$diff_results$Protein)
      } else { # Assume correlation-based
        up_genes <- sig_df %>% filter(spearman_rho > 0) %>% pull(Protein)
        down_genes <- sig_df %>% filter(spearman_rho < 0) %>% pull(Protein)
        rank_vec <- setNames(diff_obj$diff_results$spearman_rho, diff_obj$diff_results$Protein)
      }

      # 3. Run ORA (uses gene symbols directly via cached TERM2GENE)
      self$ora_up <- self$run_comprehensive_ora(up_genes, universe, pval_cutoff = pval_cutoff)
      self$ora_down <- self$run_comprehensive_ora(down_genes, universe, pval_cutoff = pval_cutoff)

      # 4. Run GSEA (uses gene symbols directly)
      rank_vec <- sort(rank_vec, decreasing = TRUE)
      self$gsea_res <- self$run_comprehensive_gsea(rank_vec, pval_cutoff = pval_cutoff)
      message("DiffExp Enrichment Done.")
    },

# Modified enrich_to_excel function
    enrich_to_excel = function(enrich_obj = NULL, direction = "UP", output_prefix = "DiffExp", target_dir = NULL){
      
      # --- Logic 1: Auto-select internal object by direction ---
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
      
      # Check if object is null (analysis may not have been run yet)
      if(is.null(enrich_obj)) {
        warning(paste("No enrichment results found for direction:", direction, "- Skipping Excel export."))
        return(NULL)
      }

      # --- Logic 2: Prepare Excel data ---
      excel_list <- list()
      databases = c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Wiki", "Reactome")

      for(db in databases){
        # Add direction prefix to sheet names, e.g. UP_GO_BP
        list_name = paste0(direction, "_", db)
        
        # Check if this database has results
        if(!is.null(enrich_obj[[db]]) && nrow(enrich_obj[[db]]) > 0) {
          excel_list[[list_name]] <- as.data.frame(enrich_obj[[db]])
        }
      }
      
      if(length(excel_list) == 0) {
        warning("Enrichment object is empty (no significant pathways). No Excel file created.")
        return(NULL)
      }
      
      # --- Logic 3: Handle save path and target_dir ---
      # Build filename (add direction to prevent overwrite)
      file_name <- paste0(output_prefix, "_", direction, ".xlsx")
      
      final_path <- file_name
      if(!is.null(target_dir)) {
        # Auto-create directory if it doesn't exist
        if(!dir.exists(target_dir)) {
          dir.create(target_dir, recursive = TRUE)
          message(paste("Created directory:", target_dir))
        }
        final_path <- here::here(target_dir, file_name)
      }
      
      # Save
      write.xlsx(excel_list, file = final_path, overwrite = TRUE)
      message(paste("Enrichment results saved to:", final_path))
    },

    gsea_to_excel = function(enrich_obj = NULL, output_prefix = "DiffExp", target_dir = NULL){
      # --- Logic 1: Auto-select internal GSEA object ---
      if(is.null(enrich_obj)){
        enrich_obj <- self$gsea_res
        message("Selected self$gsea_res for export.")
      }
      
      # Check if still null (neither passed as parameter nor computed internally)
      if(is.null(enrich_obj) ) {
        stop("GSEA results are empty. Please run analyze_diff_obj() first or provide an enrichment object.")
      }

      # --- Logic 2: Prepare Excel data ---
      databases = c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Wiki", "Reactome")
      excel_gsea_list <- list()

      for(db in databases){
        list_name = paste0("GSEA", "_", db)
        # Check nrow > 0 to prevent writing empty sheets
        if(!is.null(enrich_obj[[db]]) && nrow(enrich_obj[[db]]) > 0) {
          excel_gsea_list[[list_name]] <- as.data.frame(enrich_obj[[db]])
        }
      }
      
      if(length(excel_gsea_list) == 0) {
        warning("GSEA object is empty (no significant pathways found). No Excel file created.")
        return(NULL)
      }

      # --- Logic 3: Handle save path and target_dir ---
      file_name <- paste0(output_prefix, "_GSEA.xlsx")
      final_path <- file_name
      
      if(!is.null(target_dir)) {
        # Auto-create directory if it doesn't exist
        if(!dir.exists(target_dir)) {
          dir.create(target_dir, recursive = TRUE)
          message(paste("Created directory:", target_dir))
        }
        final_path <- file.path(target_dir, file_name)
      }

      # Save file
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

    # --- [Added] Feature: Plot protein expression trends for a specific pathway ---
    plot_pathway_trend = function(pathway_id,    # Pathway ID (e.g. "hsa04110") or Description (partial match)
                                  se_obj,        # Must provide SE object containing expression matrix
                                  meta_col,      # Column name in metadata (e.g. "age", "stage")
                                  direction = "UP", # "UP" or "DOWN"
                                  db_name = "KEGG"  # "GO_BP", "KEGG", "Reactome" etc.
                                  ) {
      
      # 1. Get enrichment result object
      target_res_list <- if (toupper(direction) == "UP") self$ora_up else self$ora_down
      
      if (is.null(target_res_list) || is.null(target_res_list[[db_name]])) {
        stop(paste0("No enrichment results found for direction: ", direction, " and DB: ", db_name))
      }
      
      res_df <- as.data.frame(target_res_list[[db_name]])
      
      # 2. Find pathway
      # Try exact ID match, fallback to Description matching
      target_row <- res_df %>% filter(ID == pathway_id)
      
      if (nrow(target_row) == 0) {
        # Try Description partial matching
        target_row <- res_df %>% filter(grepl(pathway_id, Description, ignore.case = TRUE))
      }
      
      if (nrow(target_row) == 0) {
        stop(paste0("Pathway '", pathway_id, "' not found in ", direction, " ", db_name, " results."))
      } else if (nrow(target_row) > 1) {
        warning(paste0("Multiple pathways matched '", pathway_id, "'. Using the first one: ", target_row$Description[1]))
        target_row <- target_row[1, ]
      }
      
      pathway_name <- target_row$Description
      gene_str <- target_row$geneID # Assumed to be Symbol (readable=TRUE)
      
      message(paste0("Plotting trend for: ", pathway_name, " (", target_row$ID, ")"))
      
      # 3. Parse genes and extract expression data
      target_genes <- unlist(strsplit(gene_str, "/")) # clusterProfiler default separator is "/"
      
      # Check genes in SE object
      expr_mat <- assay(se_obj)
      valid_genes <- intersect(target_genes, rownames(expr_mat))
      
      if (length(valid_genes) == 0) stop("None of the pathway genes found in the expression matrix.")
      if (length(valid_genes) < 3) warning("Less than 3 genes found for this pathway. Trends might be unstable.")
      
      # 4. Z-score standardization (row/protein scaling)
      # t(apply(..., 1, scale)) scales each protein to mean=0, sd=1 across samples
      subset_mat <- expr_mat[valid_genes, , drop=FALSE]
      scaled_mat <- t(apply(subset_mat, 1, scale))
      colnames(scaled_mat) <- colnames(expr_mat) # scale drops column names, restore them
      
      # 5. Integrate metadata
      meta <- as.data.frame(colData(se_obj))
      
      if (!meta_col %in% colnames(meta)) stop(paste0("Column '", meta_col, "' not found in metadata."))
      
      plot_df <- as.data.frame(scaled_mat) %>%
        rownames_to_column("Gene") %>%
        pivot_longer(cols = -Gene, names_to = "SampleID", values_to = "Z_Score")
      
      # Add metadata info
      # Ensure SampleID matches metadata rownames or ID column
      if(!all(plot_df$SampleID %in% rownames(meta))) {
         # Try matching metadata ID column (if available), otherwise assume rowname is sample ID
         # Simple approach: assume metadata rownames are sample IDs
         warning("Sample IDs in matrix do not fully match metadata rownames. Attempting matching...")
      }
      
      plot_df$MetaVal <- meta[plot_df$SampleID, meta_col]
      
      # 6. Plot (auto-detect X-axis type)
      p <- ggplot(plot_df, aes(x = MetaVal, y = Z_Score))
      
      if (is.numeric(plot_df$MetaVal)) {
        # --- Continuous variable (e.g. Age): scatter + Loess curve ---
        p <- p +
          # Background lines: individual trend per gene
          geom_line(aes(group = Gene), alpha = 0.15, color = "grey60", linewidth = 0.3) +
          # Overall trend line
          geom_smooth(method = "loess", color = "#B31B21", fill = "#B31B21", alpha = 0.2, linewidth = 1.5) +
          # Scatter points
          geom_point(alpha = 0.3, size = 1, color = "#1465AC") +
          labs(x = meta_col)
          
      } else {
        # --- Categorical variable (e.g. Stage): boxplot ---
        p <- p +
          geom_violin(aes(fill = MetaVal), alpha = 0.2, color = NA) +
          geom_boxplot(aes(color = MetaVal), width = 0.2, fill = "white", outlier.shape = NA) +
          geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) +
          scale_color_brewer(palette = "Set1") +
          scale_fill_brewer(palette = "Set1") +
          labs(x = meta_col)
      }
      
      # 7. Common formatting
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
