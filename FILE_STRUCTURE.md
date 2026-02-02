# ProteomicsApp 项目文件结构

```
ProteomicsApp/
├── .Rbuildignore                    # R包构建忽略文件
├── .gitignore                       # Git忽略文件
├── DESCRIPTION                      # 包描述文件（需更新作者信息）
├── NAMESPACE                        # 命名空间（自动生成）
├── LICENSE                          # MIT许可证
├── README.md                        # 项目说明（需更新GitHub链接）
├── NEWS.md                          # 更新日志
├── USER_GUIDE.md                    # 详细使用指南
├── INSTALL.md                       # 安装说明
├── DEPLOYMENT.md                    # 部署清单（本文件）
│
├── R/                               # R函数目录
│   └── run_app.R                   # 主启动函数
│
├── inst/                            # 安装资源目录
│   ├── app/                        # Shiny应用
│   │   └── app.R                   # 主应用文件（已适配R包环境）
│   ├── scripts/                    # 核心分析脚本
│   │   ├── ProteomicsAnalysis.R   # 数据管理和差异分析
│   │   ├── ProteomicsGSVA.R       # GSVA通路分析
│   │   └── run_proteomics_pipeline.R  # 完整分析流程
│   └── extdata/                    # 示例数据（可选）
│
├── man/                             # 文档目录
│   └── run_app.Rd                  # run_app函数文档（自动生成）
│
├── install_package.R                # 自动安装脚本
├── build_package.R                  # 构建脚本
├── test_package.R                   # 测试脚本
└── upload_to_github.sh              # GitHub上传脚本（可执行）

生成的文件（构建后）/
├── ProteomicsApp_1.0.0.tar.gz      # 源码包
└── ProteomicsApp_1.0.0.tgz         # 二进制包（macOS/Linux）
    或 ProteomicsApp_1.0.0.zip      # 二进制包（Windows）
```

## 文件说明

### 核心文件（必需）

1. **DESCRIPTION**
   - 包的元信息
   - **需要更新**: 作者、邮箱、GitHub URL

2. **NAMESPACE**
   - 定义导出的函数
   - 由roxygen2自动生成，不要手动编辑

3. **R/run_app.R**
   - 唯一导出的函数
   - 启动Shiny应用

4. **inst/app/app.R**
   - 完整的Shiny应用代码
   - 已修改为适配R包环境

5. **inst/scripts/**
   - 三个核心R6类和分析流程
   - 不需要修改

### 文档文件

1. **README.md**
   - 项目首页
   - **需要更新**: 将所有`your-username`替换为实际GitHub用户名

2. **USER_GUIDE.md**
   - 详细使用教程
   - 包含所有功能说明

3. **INSTALL.md**
   - 安装指南
   - 包含故障排除

4. **NEWS.md**
   - 版本更新日志

5. **DEPLOYMENT.md**
   - 部署检查清单

### 辅助脚本

1. **install_package.R**
   - 一键安装所有依赖和包本身
   - 推荐团队成员使用

2. **build_package.R**
   - 自动构建和检查包
   - 生成.tar.gz文件

3. **test_package.R**
   - 快速测试包是否正常工作

4. **upload_to_github.sh**
   - 交互式GitHub上传助手
   - 需要可执行权限

## 使用场景

### 开发者场景

```bash
# 1. 修改代码
# 编辑 inst/app/app.R 或其他文件

# 2. 重新生成文档
Rscript -e "devtools::document()"

# 3. 测试
Rscript test_package.R

# 4. 构建
Rscript build_package.R

# 5. 提交
git add .
git commit -m "Your changes"
git push
```

### 团队成员场景

```r
# 首次安装
devtools::install_github("your-username/ProteomicsApp")

# 更新到最新版
devtools::install_github("your-username/ProteomicsApp", force = TRUE)

# 使用
library(ProteomicsApp)
run_app()
```

### 离线分发场景

```bash
# 1. 构建包
Rscript build_package.R

# 2. 分享文件
# 将生成的 ProteomicsApp_1.0.0.tar.gz 发送给同事

# 3. 同事安装
# 在R中运行:
# install.packages("ProteomicsApp_1.0.0.tar.gz", repos = NULL, type = "source")
```

## 重要提示

### ⚠️ 部署前必须修改

1. **DESCRIPTION** (第5-6行)
   ```r
   Authors@R: person("Your", "Name", email = "your.email@example.com", ...)
   ```

2. **README.md** 中的所有链接
   - 搜索并替换 `your-username`
   - 搜索并替换 `your.email@example.com`

### ✅ 不需要修改

- NAMESPACE（自动生成）
- man/run_app.Rd（自动生成）
- inst/scripts/下的R文件（除非要改进算法）

### 📦 构建产物

构建后会生成：
- `ProteomicsApp_1.0.0.tar.gz`: 源码包，跨平台
- `ProteomicsApp_1.0.0.tgz` 或 `.zip`: 二进制包，特定平台

分享时推荐使用源码包（.tar.gz），因为它适用于所有平台。

---

**准备就绪！** 按照DEPLOYMENT.md的步骤进行部署。
