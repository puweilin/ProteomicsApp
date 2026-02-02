# ProteomicsApp 使用指南

## 快速开始

### 1. 安装包

```r
# 从GitHub安装
devtools::install_github("your-username/ProteomicsApp")

# 或从本地源码安装
source("install_package.R")
```

### 2. 启动应用

```r
library(ProteomicsApp)
run_app()
```

### 3. 准备数据

应用需要Excel文件（.xlsx格式），包含三个工作表：

#### Sheet1: Matrix (表达矩阵)
```
| Protein_ID  | Sample1 | Sample2 | Sample3 | ... |
|-------------|---------|---------|---------|-----|
| P12345      | 23.5    | 24.1    | NA      | ... |
| Q67890      | 18.2    | 17.9    | 18.5    | ... |
```

#### Sheet2: Metadata (样本信息)
```
| sample_id | condition | tissue | age | ... |
|-----------|-----------|--------|-----|-----|
| Sample1   | Control   | Liver  | 25  | ... |
| Sample2   | Treatment | Liver  | 27  | ... |
```

#### Sheet3: Annotation (蛋白注释)
```
| Protein   | name  | description          | ... |
|-----------|-------|----------------------|-----|
| P12345    | ACTB  | Actin beta          | ... |
| Q67890    | GAPDH | Glyceraldehyde-3... | ... |
```

## 详细使用流程

### 第一步：数据预处理 (Preprocessing)

#### 1.1 上传数据

1. 点击 "Browse" 上传Excel文件
2. 选择过滤条件（可选）
   - 例如：Filter by "tissue" = "Liver"
   - 或选择 "None" 使用所有样本
3. 设置缺失值阈值（推荐：70%）
4. 点击 "Load Data"

#### 1.2 查看数据概览

在 "Data Overview" 标签页查看：
- 样本数量和蛋白数量
- 过滤统计信息
- 缺失值模式图

#### 1.3 数据插补

1. 选择插补方法：
   - **missForest**: 推荐，精度高（较慢）
   - **MinProb**: 适合MNAR数据
   - **QRILC**: 适合MNAR数据
   - **KNN**: 适合MAR数据

2. 点击 "Run Imputation"

3. 在 "Imputation QC" 查看质量评估

#### 1.4 异常值检测与移除（可选）

1. 点击 "Auto-Detect (PCA)" 自动检测异常样本
2. 在 "PCA (Outlier Check)" 查看PCA图
3. 手动调整要移除的样本（可选）
4. 点击 "Confirm Update" 应用更改

#### 1.5 查看Clean Matrix

在 "Clean Matrix" 标签页：
- 查看完整的插补后矩阵
- 红色标记的值为插补值
- 可下载为Excel文件

### 第二步：分析配置 (Analysis Config)

#### 2.1 选择分析模式

三种模式可选：

**Group Comparison (组间比较)**
- 用于两组样本的差异分析
- 示例：Treatment vs Control
- 需要设置：Control组、Case组

**ANOVA (多组方差分析)**
- 用于三组及以上的比较
- 示例：Tissue A vs B vs C
- 自动进行成对比较

**Continuous Regression (连续变量回归)**
- 用于连续变量的相关性分析
- 示例：蛋白表达与年龄的关系
- 计算Spearman相关系数

#### 2.2 设置参数

1. **Condition Column**: 选择分组变量
2. **Control/Case**: 选择对照组和实验组（Group模式）
3. **Covariates**: 选择协变量（可选）
4. **阈值设置**:
   - P-val Cutoff: 0.05（推荐）
   - LogFC Cutoff: 0.58（推荐，约1.5倍）
   - Use Adjusted P-val: 勾选（推荐）
5. **GSVA**: 是否运行通路分析（耗时较长）

#### 2.3 运行分析

1. 检查所有参数设置
2. 点击 "Run New Analysis"
3. 等待分析完成（显示进度）
4. 自动跳转到可视化页面

#### 2.4 加载历史结果（可选）

如果之前运行过相同参数的分析：
1. 在 "Found History" 下拉菜单选择
2. 点击 "Load Selected"

### 第三步：可视化与结果导出 (Visualization)

#### 3.1 火山图与差异蛋白表格

**Volcano & Table** 标签页：

1. **交互式火山图**
   - X轴：LogFC（或相关系数）
   - Y轴：-log10(P-value)
   - 红色：上调；蓝色：下调
   - 使用鼠标框选感兴趣的点

2. **调整阈值**
   - 使用侧边栏滑块实时调整
   - 观察显著蛋白数量变化

3. **差异蛋白表格**
   - 完整的统计结果
   - 可排序、筛选、搜索
   - 导出为CSV或Excel

#### 3.2 热图

**Heatmap** 标签页：

1. 设置参数：
   - Top N Proteins: 显示前N个差异蛋白
   - Cluster Rows: 是否聚类行
   - Cluster Cols: 是否聚类列

2. 热图解读：
   - 颜色：蓝色（低表达）→ 白色 → 红色（高表达）
   - 行：蛋白质
   - 列：样本（按条件分组）

#### 3.3 PCA图

**PCA** 标签页：

- 查看样本在主成分空间的分布
- 颜色表示分组
- 评估样本分组效果
- 检测批次效应

#### 3.4 单蛋白深度分析

**Deep Dive** 标签页：

1. 在侧边栏选择感兴趣的蛋白
2. 查看该蛋白在不同条件下的表达
3. 包含统计检验结果

#### 3.5 富集分析

**Enrichment (ORA)** 标签页：

1. 选择数据库（GO_BP, KEGG, Reactome等）
2. 选择方向（上调或下调蛋白）
3. 查看富集的通路/功能
4. 点状图显示富集结果
5. 表格显示详细信息

**Enrichment (GSEA)** 标签页：

1. 选择数据库
2. 选择特定通路查看富集曲线
3. 点状图显示所有显著通路

**GSVA Analysis** 标签页（如果运行了GSVA）：

1. 选择数据库
2. 选择可视化类型：
   - Volcano: 通路水平的火山图
   - Heatmap: 通路表达热图
   - Barplot: Top通路柱状图

#### 3.6 下载报告

点击侧边栏的 "Download Report" 下载：
- 所有可视化图表（PDF格式）
- 差异蛋白表格（Excel格式）
- 富集分析结果（Excel格式）

## 高级技巧

### 1. 缓存管理

应用会自动缓存处理过的数据：
- 位置：`results_web_session/文件名/过滤条件/`
- 包含：预处理数据、分析结果
- 优势：快速重新加载，节省时间

清除缓存：
```r
# 删除特定文件的缓存
unlink("results_web_session/your_file/", recursive = TRUE)

# 清除所有缓存
unlink("results_web_session/", recursive = TRUE)
```

### 2. 自定义协变量

在limma模型中包含协变量：
1. 在Metadata中添加协变量列（如age, batch）
2. 在Analysis Config页面选择协变量
3. 模型会自动调整这些因素的影响

### 3. 批量导出

使用R脚本批量处理：
```r
library(ProteomicsApp)

# 加载缓存的结果
diff_tool <- readRDS("results_web_session/.../diff_tool.rds")
enrich_tool <- readRDS("results_web_session/.../enrich_tool.rds")

# 导出表格
write.csv(diff_tool$diff_results, "results.csv")
write.csv(enrich_tool$ora_up$GO_BP@result, "enrichment_GO_BP.csv")

# 生成自定义图表
library(ggplot2)
volcano <- diff_tool$plot_volcano(logfc_cutoff = 0.58, pval_cutoff = 0.05)
ggsave("volcano.pdf", volcano, width = 8, height = 6)
```

### 4. 多个比较分析

进行多个组间比较：
1. 第一次分析：Treatment1 vs Control
2. 保存结果（自动缓存）
3. 修改参数：Treatment2 vs Control
4. 运行第二次分析
5. 使用历史记录切换查看

### 5. 自定义基因集

使用自定义通路数据库：
```r
# 准备GMT格式的基因集文件
# 在脚本中修改enrichment数据库路径
# （需要修改源码，高级用户）
```

## 常见问题解答

### Q1: 为什么插补很慢？

A: missForest使用随机森林，计算量大。建议：
- 使用多核：脚本中设置`cores = 4`
- 或选择更快的方法：MinProb, QRILC

### Q2: GSVA一直在运行？

A: GSVA需要大量计算，特别是大数据集。建议：
- 首次分析时取消勾选 "Run GSVA"
- 先查看差异分析结果
- 需要时单独运行GSVA

### Q3: 如何解读火山图？

A:
- **右上角（红色）**: 显著上调蛋白
- **左上角（蓝色）**: 显著下调蛋白
- **中间（灰色）**: 非显著变化
- **阈值线**: 虚线表示LogFC和P-value阈值

### Q4: 富集分析没有结果？

A: 可能原因：
- 显著蛋白太少（<10个）
- 基因名注释不完整
- 选择的数据库不适合
- P-value阈值太严格

解决方案：
- 放宽阈值（P < 0.1）
- 检查基因名格式（应为Gene Symbol）
- 尝试不同数据库

### Q5: 如何处理批次效应？

A:
1. 在Metadata中添加batch列
2. 在Analysis Config中添加batch为协变量
3. limma会自动调整批次效应

或使用ComBat批次矫正（需要额外脚本）。

### Q6: 内存不足怎么办？

A:
- 减少样本数（使用过滤功能）
- 提高缺失值阈值（过滤更多蛋白）
- 关闭GSVA
- 增加系统内存

## 最佳实践

### 1. 数据准备

- **样本数**: 每组至少3个生物学重复
- **缺失值**: 尽量<30%
- **注释**: 使用标准Gene Symbol
- **元数据**: 包含所有相关变量

### 2. 参数选择

- **缺失值阈值**: 50-70%
- **插补方法**: missForest（小数据）或MinProb（大数据）
- **P-value**: 0.05，使用校正后P值
- **LogFC**: 0.58（1.5倍）或1（2倍）

### 3. 结果验证

- 检查PCA图，确保分组清晰
- 查看Top蛋白，确认生物学意义
- 与已知标志物对比
- Western Blot验证关键蛋白

### 4. 报告撰写

建议包含：
- 样本信息表
- 数据质控图（PCA, 缺失值分布）
- 火山图和热图
- Top差异蛋白表格
- 富集分析结果
- 方法学描述（引用使用的工具）


---

**祝使用愉快！**
