#!/usr/bin/env Rscript
# =============================================================================
# ProteomicsApp 安装脚本
# =============================================================================
#
# 用途：自动安装ProteomicsApp及其所有依赖包
# 使用：在R中运行 source("install_package.R")
#      或在终端运行 Rscript install_package.R
#
# =============================================================================

cat("\n=============================================================================\n")
cat("ProteomicsApp 安装程序\n")
cat("=============================================================================\n\n")

# 设置CRAN镜像
options(repos = c(CRAN = "https://cloud.r-project.org"))

# 步骤1: 安装devtools
cat("步骤 1/4: 检查并安装 devtools...\n")
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
  cat("✓ devtools 安装完成\n\n")
} else {
  cat("✓ devtools 已安装\n\n")
}

# 步骤2: 安装BiocManager
cat("步骤 2/4: 检查并安装 BiocManager...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  cat("✓ BiocManager 安装完成\n\n")
} else {
  cat("✓ BiocManager 已安装\n\n")
}

# 步骤3: 安装CRAN依赖包
cat("步骤 3/4: 安装 CRAN 依赖包...\n")
cran_packages <- c(
  "shiny", "shinydashboard", "shinyWidgets", "shinyjs", "shinybusy", "shinythemes",
  "tidyverse", "dplyr", "tidyr", "readxl", "openxlsx", "DT", "here",
  "ggplot2", "ggrepel", "patchwork", "R6", "broom",
  "missForest", "doParallel", "glmnet", "caret", "splines",
  "pheatmap", "viridis", "digest", "jsonlite", "ggplotify", "grid"
)

missing_cran <- cran_packages[!sapply(cran_packages, requireNamespace, quietly = TRUE)]

if (length(missing_cran) > 0) {
  cat(sprintf("需要安装 %d 个 CRAN 包...\n", length(missing_cran)))
  install.packages(missing_cran, dependencies = TRUE)
  cat("✓ CRAN 包安装完成\n\n")
} else {
  cat("✓ 所有 CRAN 包已安装\n\n")
}

# 步骤4: 安装Bioconductor依赖包
cat("步骤 4/4: 安装 Bioconductor 依赖包...\n")
bioc_packages <- c(
  "DEP", "SummarizedExperiment", "limma",
  "clusterProfiler", "enrichplot", "org.Hs.eg.db", "ReactomePA",
  "Mfuzz", "ComplexHeatmap", "GSVA", "GSEABase"
)

missing_bioc <- bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]

if (length(missing_bioc) > 0) {
  cat(sprintf("需要安装 %d 个 Bioconductor 包...\n", length(missing_bioc)))
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
  cat("✓ Bioconductor 包安装完成\n\n")
} else {
  cat("✓ 所有 Bioconductor 包已安装\n\n")
}

# 步骤5: 安装ProteomicsApp
cat("步骤 5/5: 构建并安装 ProteomicsApp...\n")

# 检测是否在ProteomicsApp目录中
if (file.exists("DESCRIPTION") && file.exists("R/run_app.R")) {
  # 在包目录中，直接安装
  devtools::install(dependencies = TRUE, upgrade = "never")
  cat("✓ ProteomicsApp 安装完成\n\n")
} else if (file.exists("ProteomicsApp/DESCRIPTION")) {
  # 在包的父目录中
  devtools::install("ProteomicsApp", dependencies = TRUE, upgrade = "never")
  cat("✓ ProteomicsApp 安装完成\n\n")
} else {
  cat("⚠ 警告: 未找到ProteomicsApp包目录\n")
  cat("请确保您在正确的目录中运行此脚本\n\n")
  stop("安装中止")
}

# 验证安装
cat("=============================================================================\n")
cat("验证安装...\n")
cat("=============================================================================\n\n")

if (requireNamespace("ProteomicsApp", quietly = TRUE)) {
  cat("✓ ProteomicsApp 安装成功！\n\n")
  cat("使用方法:\n")
  cat("  library(ProteomicsApp)\n")
  cat("  run_app()\n\n")
  cat("更多信息请查看: ?run_app\n\n")
  cat("=============================================================================\n")
  cat("安装完成！\n")
  cat("=============================================================================\n\n")
} else {
  cat("✗ 安装验证失败\n")
  cat("请检查错误信息并重试\n\n")
  stop("安装验证失败")
}
