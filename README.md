# ProteomicsApp

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

一个全面的蛋白质组学数据分析Shiny应用，提供数据预处理、质量控制、差异表达分析和通路富集分析等功能。

## 功能特性

### 数据预处理 (Preprocessing)
- 📊 **数据导入**: 支持Excel文件导入（矩阵、元数据、注释）
- 🔍 **质量控制**: 可自定义缺失值阈值过滤
- 🧬 **数据插补**: 支持多种插补方法（missForest, MinProb, QRILC, KNN等）
- 📈 **异常值检测**: PCA自动检测和手动选择
- 💾 **缓存系统**: 智能缓存，加速重复分析

### 分析配置 (Analysis Config)
- 🔬 **多种分析模式**:
  - 组间比较 (Group Comparison)
  - 多组方差分析 (ANOVA)
  - 连续变量回归 (Continuous Regression)
- ⚙️ **灵活参数设置**: P值、LogFC、相关系数等阈值可调
- 📝 **结果管理**: 自动保存分析结果，支持历史记录加载
- 🧮 **GSVA分析**: 可选的基因集变异分析

### 可视化 (Visualization)
- 🌋 **火山图**: 交互式差异蛋白展示
- 🔥 **热图**: Top差异蛋白表达模式
- 🎯 **PCA图**: 样本分组可视化
- 📊 **单蛋白深度分析**: 指定蛋白的表达趋势
- 🧬 **富集分析**:
  - ORA (Over-Representation Analysis)
  - GSEA (Gene Set Enrichment Analysis)
  - GSVA可视化
- 📥 **报告导出**: 一键下载分析报告

## 安装

### 系统要求

- R >= 4.0.0
- Bioconductor >= 3.12

### 从GitHub安装

```r
# 安装devtools（如果尚未安装）
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# 从GitHub安装ProteomicsApp
devtools::install_github("your-username/ProteomicsApp")
```

### 本地安装

如果您有本地包文件：

```r
# 方法1：从源码安装
devtools::install_local("path/to/ProteomicsApp")

# 方法2：构建并安装
# 在终端中：
# cd /path/to/ProteomicsApp
# R CMD build .
# R CMD INSTALL ProteomicsApp_1.0.0.tar.gz
```

### 安装依赖

首次使用前，建议运行以下脚本安装所有依赖包：

```r
# CRAN包
cran_packages <- c(
  "shiny", "shinydashboard", "shinyWidgets", "shinyjs", "shinybusy", "shinythemes",
  "tidyverse", "readxl", "openxlsx", "DT", "here",
  "ggplot2", "ggrepel", "patchwork", "R6", "broom",
  "missForest", "doParallel", "glmnet", "caret", "splines",
  "pheatmap", "viridis", "digest", "jsonlite", "ggplotify"
)

install.packages(cran_packages)

# Bioconductor包
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "DEP", "SummarizedExperiment", "limma",
  "clusterProfiler", "enrichplot", "org.Hs.eg.db", "ReactomePA",
  "Mfuzz", "ComplexHeatmap", "GSVA", "GSEABase"
)

BiocManager::install(bioc_packages)
```

## 快速开始

### 启动应用

```r
library(ProteomicsApp)

# 启动应用
run_app()

# 在指定端口启动
run_app(port = 8080)

# 在浏览器中启动
run_app(launch.browser = TRUE)
```

### 数据格式要求

应用需要Excel文件（.xlsx），包含三个工作表：

1. **Sheet1 (Matrix)**: 蛋白表达矩阵
   - 第一列：Protein ID
   - 其余列：样本表达值

2. **Sheet2 (Metadata)**: 样本元数据
   - 必须包含样本名列（与Matrix列名对应）
   - 至少包含一个分组列（如condition、tissue等）

3. **Sheet3 (Annotation)**: 蛋白注释
   - Protein列：蛋白ID
   - name列：基因名称
   - 其他注释信息（可选）

### 基本工作流程

1. **数据导入**
   - 上传Excel文件
   - 选择过滤列和值（可选）
   - 设置缺失值阈值
   - 点击"Load Data"

2. **数据预处理**
   - 选择插补方法
   - 运行插补
   - 检查PCA图，必要时删除异常值

3. **配置分析**
   - 选择分析模式
   - 设置分组和协变量
   - 调整统计阈值
   - 运行分析

4. **可视化与导出**
   - 浏览各种可视化结果
   - 调整阈值查看不同筛选条件
   - 下载报告和结果表格

## 项目结构

```
ProteomicsApp/
├── DESCRIPTION              # 包描述文件
├── NAMESPACE               # 命名空间
├── LICENSE                 # MIT许可证
├── README.md              # 本文件
├── R/                     # R函数
│   └── run_app.R         # 启动函数
├── inst/                  # 安装资源
│   ├── app/              # Shiny应用
│   │   └── app.R        # 主应用文件
│   ├── scripts/          # 核心分析脚本
│   │   ├── ProteomicsAnalysis.R
│   │   ├── ProteomicsGSVA.R
│   │   └── run_proteomics_pipeline.R
│   └── extdata/          # 示例数据（可选）
└── man/                   # 文档
    └── run_app.Rd
```

## 主要功能模块

### ProteomicsDataManager (R6类)
- 数据加载与预处理
- 质量控制与过滤
- 多种插补方法
- PCA分析和异常值检测

### DiffExpAnalyst (R6类)
- 差异表达分析（limma）
- 组间比较、ANOVA、回归分析
- 火山图、MA图等可视化

### EnrichmentAnalyst (R6类)
- ORA富集分析
- GSEA分析
- 支持GO、KEGG、Reactome、WikiPathways等数据库

### ProteomicsGSVA (R6类)
- 基因集变异分析
- 通路水平的差异分析
- 多种可视化选项

## 常见问题

**Q: 如何处理大数据集？**
A: 应用支持最大200MB的上传文件。对于更大的数据集，建议预先过滤低质量蛋白。

**Q: 插补方法如何选择？**
A:
- missForest：推荐，精度高但较慢
- MinProb/QRILC：适合MNAR（缺失非随机）数据
- KNN：适合MAR（缺失随机）数据

**Q: 缓存数据存储在哪里？**
A: 缓存保存在`results_web_session/`目录下，按文件名和过滤条件组织。

**Q: 如何自定义分析参数？**
A: 在Analysis Config页面可调整所有统计阈值和分析选项。

## 更新日志

### v1.0.0 (2025-02-02)
- 初始版本发布
- 完整的数据预处理流程
- 三种分析模式支持
- 丰富的可视化功能
- 智能缓存系统

## 贡献

欢迎提交Issue和Pull Request！

## 许可证

MIT License - 详见[LICENSE](LICENSE)文件

## 联系方式

- 问题反馈: [GitHub Issues](https://github.com/your-username/ProteomicsApp/issues)
- 邮件: your.email@example.com

## 引用

如果您在研究中使用了ProteomicsApp，请引用：

```
Your Name (2025). ProteomicsApp: A Comprehensive Proteomics Analysis Platform.
R package version 1.0.0. https://github.com/your-username/ProteomicsApp
```

---

**提示**: 首次启动应用可能需要较长时间加载依赖包，请耐心等待。
