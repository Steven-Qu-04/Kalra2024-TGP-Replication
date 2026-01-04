# Kalra2024-TGP-Replication
基于黑腹果蝇及其近缘种转录组数据的跨代可塑性（TGP）顺反调控机制研究。本仓库复现了 Kalra et al. (2024) 的全因子线性模型分析流程，核心创新在于引入了多维动态 Q 值选取策略，显著提升了对效应量微弱的环境交互基因的识别敏感度，揭示了反式调控在跨代遗传响应中的主导作用。
复现文章：Kalra, S., Lanno, S., Sanchez, G. et al. cis- and trans-regulatory contributions to a hierarchy of factors influencing gene expression variation. Commun Biol 7, 1563 (2024). https://doi.org/10.1038/s42003-024-07255-6

## Project Overview / 项目概览
本项目系统解析了遗传（物种）、发育阶段及环境压力对果蝇基因表达的层级贡献。This project systematically dissects the hierarchical contributions of genetics (species), developmental stage, and environmental stress to gene expression in Drosophila.Key 
## Methodological Highlight / 核心方法亮点:
### Dynamic Q-value Selection / 动态 Q 值选取: 
不同于单一的 FDR 阈值，本项目通过交叉比对 $Q_P$, $Q_H$, 和 $Q_C$ 值，在保持低假阳性的同时，找回了 24 个在常规统计中易被忽略的关键交互基因。
### is/Trans Decomposition / 顺反调控分解:
利用 F1 杂交种的等位基因特异性表达（ASE）量化各因子的分子调控逻辑。

## 📂 Repository Structure / 目录结构
```
├── scripts/                        # Core analysis scripts / 核心分析脚本
│   ├── Main_script.R               # Main pipeline execution / 主流程执行脚本
│   ├── Custom_functions.R          # Helper functions for Q-value & linear models / 自定义功能函数
│   ├── Updated_counts.R            # Data normalization & processing / 统计计数与预处理
│   ├── limma_ana.R                 # Limma-based differential expression analysis / 差异表达分析
│   └── draw2.R                     # Visualization & plotting automation / 自动化绘图脚本
│
├── data/                           # Processed data and metadata / 处理后的数据与元数据
│   ├── Final_combined_FPKM.csv     # Merged expression matrix / 整合后的表达矩阵
│   ├── Factors_and_data.csv        # Factor design matrix (Stage/Mat/Cur) / 因子设计矩阵
│   └── 2_orthologs.txt             # Orthologous gene information / 种间同源基因信息
│
├── results/                        # Analysis outputs / 分析结果导出
│   ├── 1_PCA/                      # Principal Component Analysis plots / PCA 聚类分析图
│   ├── 2_LinearModel/              # Statistical results (Q-values, Estimates) / 线性模型统计结果
│   │   ├── Parent/                 # Parental species divergence / 亲本种间分化结果
│   │   ├── Hybrid/                 # Hybrid ASE analysis / 杂交种等位基因特异性表达分析
│   │   └── Combined/               # Full-factorial interaction results / 全因子交互作用分析
│   ├── 3_Allele_Interactions/      # Cis/Trans mechanism classification / 顺反调控机制判定
│   └── 4_OddRatios/                # Enrichment & Fisher's test results / 费舍尔精确检验与优势比
│
├── images/                         # Generated visualizations / 论文与报告所用图表
│   ├── 0_Sample_Clustering.png     # Sample hierarchical clustering / 样本层级聚类图
│   ├── 1_DEGs_Summary_Barplot.png  # Summary of Differentially Expressed Genes / DEGs 统计图
│   ├── 2_Mechanism_Stacked_Barplot.png # Cis/Trans distribution / 顺反调控分布堆叠图
│   └── Venn_Diagram_Results.png    # Overlap of different factors / 不同因子影响基因的韦恩图
│
├── raw_counts/                     # Source RNA-seq FPKM files / 原始 FPKM 数据
│   ├── Adults/                     # Adult stage samples (sec/sim/hyb) / 成虫样本
│   └── Larvae/                     # Larval stage samples (sec/sim/hyb) / 幼虫样本
│
└── docs/
    └── 42003_2024_Article_7255.pdf # Reference paper (Kalra et al., 2024) / 参考原文文献
```
## 📝 File Descriptions / 文件详细说明
Main_script.R: The entry point of the project. It integrates data merging, linear modeling, and result export. (项目主程序，整合数据合并、模型构建与结果导出)
Custom_functions.R: Contains the logic for the Dynamic Q-value selection strategy, essential for identifying subtle interaction effects. (包含本项目核心的“动态 Q 值”选取逻辑)
2_LinearModel/: Holds all statistical parameters ($P, Q, T$ values) across three different model levels: Parent, Hybrid, and Combined. (存储 Parent、Hybrid 和 Combined 三个层级的全套统计参数)
images/: High-resolution plots ready for academic reporting. (存储符合学术报告要求的高分辨率图表)
