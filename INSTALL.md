# ProteomicsApp - GitHub安装指南

## 快速安装（推荐）

### 方法1: 直接从GitHub安装

```r
# 1. 安装devtools（如果尚未安装）
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# 2. 从GitHub安装ProteomicsApp
devtools::install_github("your-username/ProteomicsApp")

# 3. 启动应用
library(ProteomicsApp)
run_app()
```

### 方法2: 克隆仓库后安装

```bash
# 1. 克隆仓库
git clone https://github.com/your-username/ProteomicsApp.git
cd ProteomicsApp

# 2. 在R中运行安装脚本
Rscript install_package.R

# 或在R控制台中
# source("install_package.R")
```

## 详细安装步骤

### 1. 准备环境

确保您的系统已安装：
- R >= 4.0.0
- RStudio（推荐）

### 2. 安装依赖管理工具

```r
# 安装 devtools
install.packages("devtools")

# 安装 BiocManager
install.packages("BiocManager")
```

### 3. 安装依赖包

#### CRAN 包

```r
cran_packages <- c(
  "shiny", "shinydashboard", "shinyWidgets", "shinyjs", "shinybusy", "shinythemes",
  "tidyverse", "readxl", "openxlsx", "DT", "here",
  "ggplot2", "ggrepel", "patchwork", "R6", "broom",
  "missForest", "doParallel", "glmnet", "caret", "splines",
  "pheatmap", "viridis", "digest", "jsonlite", "ggplotify"
)

install.packages(cran_packages)
```

#### Bioconductor 包

```r
BiocManager::install(c(
  "DEP", "SummarizedExperiment", "limma",
  "clusterProfiler", "enrichplot", "org.Hs.eg.db", "ReactomePA",
  "Mfuzz", "ComplexHeatmap", "GSVA", "GSEABase"
))
```

### 4. 安装ProteomicsApp

#### 从GitHub（推荐）

```r
devtools::install_github("your-username/ProteomicsApp")
```

#### 从本地源码

如果您已克隆仓库：

```r
# 方法A: 使用devtools
devtools::install("path/to/ProteomicsApp")

# 方法B: 使用命令行
# R CMD INSTALL ProteomicsApp
```

## 验证安装

```r
# 加载包
library(ProteomicsApp)

# 查看帮助
?run_app

# 启动应用
run_app()
```

如果成功，浏览器会自动打开应用界面。

## 升级包

### 从GitHub升级

```r
devtools::install_github("your-username/ProteomicsApp", force = TRUE)
```

### 查看版本

```r
packageVersion("ProteomicsApp")
```

## 卸载

```r
remove.packages("ProteomicsApp")
```

## 常见问题

### Q1: 安装依赖时出现编译错误

**解决方案：**
- Windows: 安装 [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
- macOS: 安装 Xcode Command Line Tools: `xcode-select --install`
- Linux: 安装开发工具: `sudo apt-get install r-base-dev` (Ubuntu/Debian)

### Q2: Bioconductor包安装失败

**解决方案：**
```r
# 更新BiocManager
install.packages("BiocManager")

# 指定Bioconductor版本
BiocManager::install(version = "3.18")

# 重新安装失败的包
BiocManager::install("package_name", force = TRUE)
```

### Q3: 网络连接问题

**解决方案：**
```r
# 使用国内镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
```

### Q4: 依赖冲突

**解决方案：**
```r
# 更新所有包
update.packages(ask = FALSE)

# 重新安装ProteomicsApp
devtools::install_github("your-username/ProteomicsApp", force = TRUE)
```

## 开发者安装

如果您想修改代码：

```bash
# 1. Fork 并克隆仓库
git clone https://github.com/your-username/ProteomicsApp.git
cd ProteomicsApp

# 2. 在RStudio中打开ProteomicsApp.Rproj

# 3. 安装开发依赖
devtools::install_dev_deps()

# 4. 加载包（不安装）
devtools::load_all()

# 5. 测试
devtools::test()

# 6. 检查包
devtools::check()

# 7. 构建文档
devtools::document()
```

## GitHub发布流程

### 1. 准备发布

```bash
# 更新版本号（在DESCRIPTION中）
# 更新NEWS.md
# 提交所有更改
git add .
git commit -m "Prepare for v1.0.0 release"
git push
```

### 2. 创建Tag

```bash
git tag -a v1.0.0 -m "Release version 1.0.0"
git push origin v1.0.0
```

### 3. 在GitHub创建Release

1. 访问 `https://github.com/your-username/ProteomicsApp/releases`
2. 点击 "Draft a new release"
3. 选择刚创建的tag（v1.0.0）
4. 填写Release标题和说明
5. 可选：上传构建的.tar.gz文件
6. 点击 "Publish release"

## 持续集成（可选）

创建 `.github/workflows/R-CMD-check.yaml`:

```yaml
on:
  push:
    branches: [main, master]
  pull_request:
    branches: [main, master]

name: R-CMD-check

jobs:
  R-CMD-check:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - uses: r-lib/actions/setup-r@v2
      - uses: r-lib/actions/setup-pandoc@v2
      - name: Install dependencies
        run: |
          install.packages(c("remotes", "rcmdcheck"))
          remotes::install_deps(dependencies = TRUE)
        shell: Rscript {0}
      - name: Check
        run: rcmdcheck::rcmdcheck(args = "--no-manual", error_on = "error")
        shell: Rscript {0}
```

## 联系支持

- 问题报告: [GitHub Issues](https://github.com/your-username/ProteomicsApp/issues)
- 功能请求: [GitHub Discussions](https://github.com/your-username/ProteomicsApp/discussions)
- 邮件: your.email@example.com
