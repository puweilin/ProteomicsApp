# 🚀 ProteomicsApp 快速开始

## 一、立即开始（3步上传到GitHub）

### 步骤1: 更新包信息

编辑 `ProteomicsApp/DESCRIPTION` 文件第5-6行：

```r
Authors@R: person("您的名字", "姓氏",
                  email = "your.email@example.com",
                  role = c("aut", "cre"))
```

### 步骤2: 测试包

在R中运行：

```r
setwd("ProteomicsApp")
source("test_package.R")
```

确保所有测试通过。

### 步骤3: 上传到GitHub

在终端运行：

```bash
cd ProteomicsApp
./upload_to_github.sh
```

按提示操作即可！

---

## 二、团队成员安装（仅需1行代码）

团队成员安装您的包：

```r
devtools::install_github("your-username/ProteomicsApp")
```

使用：

```r
library(ProteomicsApp)
run_app()
```

---

## 三、完整流程（详细版）

### 3.1 本地准备

1. **检查包结构**
   ```bash
   cd ProteomicsApp
   ls -l
   ```

   应该看到：
   - DESCRIPTION, NAMESPACE, LICENSE
   - R/, inst/, man/ 目录
   - 各种 .md 文档
   - 辅助脚本 (.R 和 .sh)

2. **安装依赖**
   ```r
   source("install_package.R")
   ```

   这会自动安装所有CRAN和Bioconductor依赖包。

3. **测试包**
   ```r
   source("test_package.R")
   ```

   检查输出，应该看到所有测试通过（✓）。

### 3.2 构建包（可选）

如果想生成可分发的包文件：

```r
source("build_package.R")
```

这会生成：
- `ProteomicsApp_1.0.0.tar.gz` - 源码包
- 可能的二进制包（.tgz 或 .zip）

### 3.3 上传到GitHub

#### 方法A: 使用上传脚本（推荐，简单）

```bash
./upload_to_github.sh
```

脚本会引导您：
1. 初始化Git（如果需要）
2. 添加所有文件
3. 输入提交信息
4. 输入GitHub仓库URL
5. 自动推送

#### 方法B: 手动上传（需要Git经验）

```bash
# 1. 在GitHub创建新仓库（不要初始化README）
# 仓库名: ProteomicsApp

# 2. 初始化Git
git init
git add .
git commit -m "Initial commit: ProteomicsApp v1.0.0"

# 3. 关联远程仓库
git remote add origin https://github.com/your-username/ProteomicsApp.git

# 4. 推送
git branch -M main
git push -u origin main
```

### 3.4 创建Release（推荐）

在GitHub仓库页面：

1. 点击 "Releases" → "Create a new release"
2. Tag: `v1.0.0`
3. Title: `ProteomicsApp v1.0.0`
4. 描述: 复制粘贴 `NEWS.md` 的内容
5. （可选）上传 `ProteomicsApp_1.0.0.tar.gz`
6. 点击 "Publish release"

---

## 四、常见问题

### Q1: 如何更新包版本？

1. 编辑 `DESCRIPTION`，修改 Version 字段
2. 更新 `NEWS.md`，添加新版本说明
3. 提交并推送
4. 创建新的Git tag

### Q2: 如何修改应用代码？

主应用代码在：`inst/app/app.R`

修改后：
1. 测试：`source("test_package.R")`
2. 提交：`git commit -m "描述修改内容"`
3. 推送：`git push`

团队成员更新：
```r
devtools::install_github("your-username/ProteomicsApp", force = TRUE)
```

### Q3: 如何分享给没有GitHub访问权限的同事？

构建包并分享文件：

```r
source("build_package.R")
# 将生成的 ProteomicsApp_1.0.0.tar.gz 发送给同事
```

同事安装：
```r
install.packages("ProteomicsApp_1.0.0.tar.gz", repos = NULL, type = "source")
```

### Q4: 报错 "cannot find package"？

确保已安装 devtools：
```r
install.packages("devtools")
```

### Q5: GitHub上传失败？

检查：
1. GitHub仓库已创建
2. 仓库URL正确
3. 有推送权限（SSH密钥或Personal Access Token已配置）

---

## 五、下一步

### 完成后的文档位置

所有文档都在 GitHub 仓库中：

- **用户指南**: `USER_GUIDE.md` - 详细使用教程
- **安装说明**: `INSTALL.md` - 安装故障排除
- **部署清单**: `DEPLOYMENT.md` - 完整部署流程
- **文件结构**: `FILE_STRUCTURE.md` - 项目文件说明

### 分享给团队

发送给团队成员：

```
大家好！

我已经将蛋白质组学分析工具打包成了R包，上传到GitHub。

安装方法：
devtools::install_github("your-username/ProteomicsApp")

使用方法：
library(ProteomicsApp)
run_app()

详细使用指南：
https://github.com/your-username/ProteomicsApp/blob/main/USER_GUIDE.md

如有问题，请在 GitHub Issues 中反馈：
https://github.com/your-username/ProteomicsApp/issues
```

### 推荐的维护流程

1. **Bug修复**
   - 创建Issue
   - 修复代码
   - 更新版本号（如1.0.0 → 1.0.1）
   - 提交并创建Release

2. **新功能**
   - 创建分支
   - 开发新功能
   - 测试
   - 合并到main
   - 更新版本号（如1.0.0 → 1.1.0）
   - 创建Release

3. **重大更新**
   - 创建分支
   - 重构或大改
   - 充分测试
   - 更新文档
   - 合并到main
   - 更新版本号（如1.0.0 → 2.0.0）
   - 创建Release

---

## 📚 所有文档索引

1. **README.md** - 项目总览和功能介绍
2. **QUICKSTART.md** (本文件) - 快速开始指南
3. **USER_GUIDE.md** - 详细使用教程
4. **INSTALL.md** - 安装说明和故障排除
5. **DEPLOYMENT.md** - 完整部署检查清单
6. **FILE_STRUCTURE.md** - 项目文件结构说明
7. **NEWS.md** - 版本更新日志

---

## ✅ 检查清单

上传前最后确认：

- [ ] 已更新 DESCRIPTION 中的作者信息
- [ ] 已运行 `test_package.R` 且测试通过
- [ ] 在GitHub创建了新仓库
- [ ] 已运行 `upload_to_github.sh` 或手动推送成功
- [ ] 在GitHub上能看到所有文件
- [ ] 已创建Release（推荐）
- [ ] 已通知团队成员

---

**🎉 完成！您的R包已经可以供团队使用了！**

有问题？查看 `USER_GUIDE.md` 或创建 GitHub Issue。
