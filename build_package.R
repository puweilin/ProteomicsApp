#!/usr/bin/env Rscript
# =============================================================================
# ProteomicsApp 构建和检查脚本
# =============================================================================

cat("\n=============================================================================\n")
cat("ProteomicsApp 构建和检查\n")
cat("=============================================================================\n\n")

# 设置工作目录
if (!file.exists("DESCRIPTION")) {
  if (file.exists("ProteomicsApp/DESCRIPTION")) {
    setwd("ProteomicsApp")
  } else {
    stop("请在ProteomicsApp包目录或其父目录中运行此脚本")
  }
}

# 加载devtools
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("请先安装devtools: install.packages('devtools')")
}

library(devtools)

# 1. 生成文档
cat("步骤 1/5: 生成文档...\n")
tryCatch({
  document()
  cat("✓ 文档生成完成\n\n")
}, error = function(e) {
  cat("✗ 文档生成失败:\n", e$message, "\n\n")
})

# 2. 检查包
cat("步骤 2/5: 运行 R CMD check...\n")
tryCatch({
  check_results <- check(
    document = FALSE,
    args = c("--no-manual", "--no-build-vignettes"),
    error_on = "warning"
  )
  cat("✓ 检查完成\n\n")
}, error = function(e) {
  cat("⚠ 检查发现问题:\n", e$message, "\n\n")
})

# 3. 构建源码包
cat("步骤 3/5: 构建源码包...\n")
tryCatch({
  pkg_file <- build(manual = FALSE, vignettes = FALSE)
  cat("✓ 源码包构建完成:", pkg_file, "\n\n")
}, error = function(e) {
  cat("✗ 构建失败:\n", e$message, "\n\n")
})

# 4. 构建二进制包（可选）
cat("步骤 4/5: 构建二进制包...\n")
tryCatch({
  bin_file <- build(binary = TRUE, manual = FALSE, vignettes = FALSE)
  cat("✓ 二进制包构建完成:", bin_file, "\n\n")
}, error = function(e) {
  cat("⚠ 二进制包构建失败（可忽略）:\n", e$message, "\n\n")
})

# 5. 测试安装
cat("步骤 5/5: 测试安装...\n")
tryCatch({
  install(upgrade = "never", quick = TRUE, build = FALSE)
  cat("✓ 安装成功\n\n")

  # 验证
  if (requireNamespace("ProteomicsApp", quietly = TRUE)) {
    cat("✓ 包验证成功\n\n")
  }
}, error = function(e) {
  cat("✗ 安装测试失败:\n", e$message, "\n\n")
})

# 总结
cat("=============================================================================\n")
cat("构建完成！\n")
cat("=============================================================================\n\n")
cat("生成的文件：\n")
if (exists("pkg_file")) cat("  - 源码包:", pkg_file, "\n")
if (exists("bin_file")) cat("  - 二进制包:", bin_file, "\n")
cat("\n")
cat("使用方法:\n")
cat("  library(ProteomicsApp)\n")
cat("  run_app()\n\n")
