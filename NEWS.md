# News and Changes

## ProteomicsApp 1.0.0 (2025-02-02)

### 首次发布

#### 新功能

- **数据预处理模块**
  - 支持Excel文件导入（矩阵、元数据、注释）
  - 自定义缺失值阈值过滤
  - 多种插补方法支持：missForest, MinProb, QRILC, KNN, BPCA等
  - PCA异常值自动检测和手动移除
  - 智能缓存系统

- **分析配置模块**
  - 三种分析模式：
    - 组间比较 (Group Comparison)
    - 多组方差分析 (ANOVA)
    - 连续变量回归 (Continuous Regression)
  - 灵活的统计参数设置
  - 分析结果历史管理
  - 可选GSVA通路分析

- **可视化模块**
  - 交互式火山图
  - Top差异蛋白热图
  - PCA样本分组图
  - 单蛋白表达深度分析
  - ORA富集分析可视化（GO, KEGG, Reactome, WikiPathways）
  - GSEA富集分析可视化
  - GSVA通路分析可视化
  - 报告一键导出

- **用户体验**
  - 现代化专业界面设计
  - 实时进度反馈
  - 详细的分析日志
  - 交互式数据表格

#### 技术特性

- R6面向对象设计
- 模块化代码结构
- 完整的R包封装
- 详细的文档和示例
- 自动依赖管理

#### 依赖包

- Shiny生态系统：shiny, shinydashboard, shinyWidgets
- 数据处理：tidyverse, dplyr, tidyr
- 统计分析：limma, DEP
- 富集分析：clusterProfiler, enrichplot, ReactomePA
- 可视化：ggplot2, pheatmap, ComplexHeatmap
- 通路分析：GSVA, GSEABase

### 未来计划

- v1.1.0: 添加更多插补方法和质控指标
- v1.2.0: 支持蛋白-蛋白互作网络分析
- v1.3.0: 添加时序分析功能
- v2.0.0: 支持多组学整合分析
