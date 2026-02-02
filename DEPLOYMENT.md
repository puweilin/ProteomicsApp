# ProteomicsApp 部署清单

## 📋 部署前检查清单

### ✅ 必需文件
- [x] DESCRIPTION - 包描述文件
- [x] NAMESPACE - 命名空间
- [x] LICENSE - MIT许可证
- [x] README.md - 项目说明
- [x] NEWS.md - 更新日志
- [x] R/run_app.R - 启动函数
- [x] inst/app/app.R - 主应用
- [x] inst/scripts/ - 核心分析脚本
  - [x] ProteomicsAnalysis.R
  - [x] ProteomicsGSVA.R
  - [x] run_proteomics_pipeline.R
- [x] man/run_app.Rd - 函数文档

### ✅ 配置文件
- [x] .Rbuildignore - 构建忽略列表
- [x] .gitignore - Git忽略列表

### ✅ 辅助脚本
- [x] install_package.R - 自动安装脚本
- [x] build_package.R - 构建脚本
- [x] test_package.R - 测试脚本
- [x] upload_to_github.sh - GitHub上传脚本

### ✅ 文档
- [x] USER_GUIDE.md - 使用指南
- [x] INSTALL.md - 安装说明

## 🚀 部署步骤

### 步骤1: 本地测试

```r
# 在R中运行
source("test_package.R")
```

检查输出，确保所有测试通过。

### 步骤2: 构建包

```r
# 在R中运行
source("build_package.R")
```

这会：
- 生成文档
- 运行R CMD check
- 构建.tar.gz源码包
- 构建二进制包（可选）

### 步骤3: 本地安装测试

```r
# 安装刚构建的包
install.packages("ProteomicsApp_1.0.0.tar.gz", repos = NULL, type = "source")

# 测试
library(ProteomicsApp)
run_app()
```

### 步骤4: 准备GitHub

1. **更新DESCRIPTION中的信息**
   ```r
   # 编辑ProteomicsApp/DESCRIPTION
   # 更新以下字段：
   # - Authors@R: 您的名字和邮箱
   # - URL: https://github.com/your-username/ProteomicsApp
   # - BugReports: https://github.com/your-username/ProteomicsApp/issues
   ```

2. **更新README.md**
   - 将所有 `your-username` 替换为您的GitHub用户名
   - 将所有 `your.email@example.com` 替换为您的邮箱

3. **检查.gitignore**
   - 确保敏感文件不会被提交

### 步骤5: 上传到GitHub

**方法A: 使用上传脚本（推荐）**

```bash
cd ProteomicsApp
./upload_to_github.sh
```

按提示操作：
1. 输入提交信息
2. 输入GitHub仓库URL
3. 等待上传完成

**方法B: 手动上传**

```bash
cd ProteomicsApp

# 初始化Git（如已初始化则跳过）
git init

# 添加文件
git add .

# 提交
git commit -m "Initial commit: ProteomicsApp v1.0.0"

# 添加远程仓库
git remote add origin https://github.com/your-username/ProteomicsApp.git

# 推送
git branch -M main
git push -u origin main
```

### 步骤6: 在GitHub上创建Release（可选但推荐）

1. 访问您的仓库页面
2. 点击 "Releases" → "Create a new release"
3. Tag version: `v1.0.0`
4. Release title: `ProteomicsApp v1.0.0`
5. 描述: 复制NEWS.md中的内容
6. 上传构建的.tar.gz文件
7. 点击 "Publish release"

### 步骤7: 团队成员安装测试

通知团队成员测试安装：

```r
# 从GitHub安装
devtools::install_github("your-username/ProteomicsApp")

# 测试
library(ProteomicsApp)
run_app()
```

## 📝 更新DESCRIPTION示例

```r
Package: ProteomicsApp
Title: Proteomics Data Analysis Platform
Version: 1.0.0
Authors@R:
    person("Your", "Name",
           email = "your.email@example.com",
           role = c("aut", "cre"),
           comment = c(ORCID = "YOUR-ORCID-ID"))
Description: A comprehensive Shiny application for proteomics data analysis,
    including data preprocessing, quality control, differential expression analysis,
    and pathway enrichment analysis.
License: MIT + file LICENSE
URL: https://github.com/your-username/ProteomicsApp
BugReports: https://github.com/your-username/ProteomicsApp/issues
Encoding: UTF-8
Roxygen: list(markdown = TRUE)
RoxygenNote: 7.3.3
Depends:
    R (>= 4.0.0)
Imports:
    shiny,
    shinydashboard,
    shinyWidgets,
    shinyjs,
    shinybusy,
    shinythemes,
    tidyverse,
    dplyr,
    tidyr,
    readxl,
    openxlsx,
    DT,
    here,
    ggplot2,
    ggrepel,
    patchwork,
    R6,
    broom,
    missForest,
    doParallel,
    glmnet,
    caret,
    splines,
    pheatmap,
    viridis,
    digest,
    jsonlite,
    ggplotify,
    grid,
    DEP,
    SummarizedExperiment,
    limma,
    clusterProfiler,
    enrichplot,
    org.Hs.eg.db,
    ReactomePA,
    Mfuzz,
    ComplexHeatmap,
    GSVA,
    GSEABase
Suggests:
    testthat (>= 3.0.0),
    knitr,
    rmarkdown
biocViews: Proteomics, DifferentialExpression, Visualization
```

## 🔍 质量检查

### 运行R CMD check

```r
devtools::check()
```

应该没有ERROR，最多有WARNING或NOTE。

### 检查文档

```r
# 生成并查看文档
devtools::document()
?run_app
```

### 测试所有功能

1. 启动应用
2. 上传测试数据
3. 完成一次完整分析流程
4. 检查所有可视化
5. 下载报告

## 📊 监控与维护

### 收集用户反馈

- 设置GitHub Issues模板
- 创建Discussion区
- 收集常见问题

### 版本更新

```bash
# 更新版本号（在DESCRIPTION中）
# 更新NEWS.md
# 提交更改
git add .
git commit -m "Version bump to 1.0.1"
git tag -a v1.0.1 -m "Release v1.0.1"
git push origin main --tags
```

### 持续改进

- 修复bug
- 添加新功能
- 优化性能
- 更新文档

## ✅ 最终检查清单

部署前最后检查：

- [ ] 所有文件已提交到Git
- [ ] DESCRIPTION信息已更新
- [ ] README中的链接已更新
- [ ] 包能成功构建（无ERROR）
- [ ] 本地安装测试通过
- [ ] 应用能正常启动
- [ ] 测试数据能正常分析
- [ ] 文档完整且准确
- [ ] LICENSE文件存在
- [ ] .gitignore配置正确

## 🎉 完成！

如果所有检查项都通过，您的包已经准备好分享给团队了！

---

**下一步**: 创建Wiki页面，添加更多示例和教程
