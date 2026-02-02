#!/usr/bin/env Rscript
# =============================================================================
# ProteomicsApp 快速测试脚本
# =============================================================================

cat("\n=============================================================================\n")
cat("ProteomicsApp 快速测试\n")
cat("=============================================================================\n\n")

# 测试是否已安装
cat("测试 1/3: 检查包是否已安装...\n")
if (!requireNamespace("ProteomicsApp", quietly = TRUE)) {
  cat("✗ ProteomicsApp 未安装\n")
  cat("请先运行: source('install_package.R')\n\n")
  stop("包未安装")
}
cat("✓ ProteomicsApp 已安装\n\n")

# 测试加载
cat("测试 2/3: 加载包...\n")
tryCatch({
  library(ProteomicsApp)
  cat("✓ 包加载成功\n\n")
}, error = function(e) {
  cat("✗ 包加载失败:\n", e$message, "\n\n")
  stop("加载失败")
})

# 测试启动函数
cat("测试 3/3: 测试 run_app 函数...\n")
if (exists("run_app")) {
  cat("✓ run_app 函数可用\n\n")
} else {
  cat("✗ run_app 函数不存在\n\n")
  stop("函数缺失")
}

# 检查依赖
cat("检查核心依赖...\n")
core_deps <- c("shiny", "tidyverse", "DEP", "clusterProfiler")
missing_deps <- core_deps[!sapply(core_deps, requireNamespace, quietly = TRUE)]

if (length(missing_deps) > 0) {
  cat("⚠ 缺少依赖包:\n")
  cat(paste("  -", missing_deps, collapse = "\n"), "\n\n")
} else {
  cat("✓ 所有核心依赖已安装\n\n")
}

# 显示版本信息
cat("=============================================================================\n")
cat("包信息:\n")
cat("=============================================================================\n")
cat("  包名称:", "ProteomicsApp", "\n")
cat("  版本号:", as.character(packageVersion("ProteomicsApp")), "\n")
cat("  R版本:", R.version.string, "\n\n")

cat("=============================================================================\n")
cat("测试完成！\n")
cat("=============================================================================\n\n")
cat("使用方法:\n")
cat("  library(ProteomicsApp)\n")
cat("  run_app()\n\n")
cat("注意：\n")
cat("  - 首次启动可能需要较长时间加载依赖\n")
cat("  - 确保您有示例数据文件（Excel格式）\n")
cat("  - 查看用户指南: USER_GUIDE.md\n\n")

# 可选：自动启动应用
cat("是否现在启动应用？[y/N]: ")
if (interactive()) {
  response <- readline()
  if (tolower(response) %in% c("y", "yes")) {
    cat("\n启动应用...\n\n")
    run_app()
  }
}
