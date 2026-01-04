# 1. 加载必要库
if (!require("pheatmap")) install.packages("pheatmap")
library(limma)
library(ggplot2)
library(pheatmap)

# 2. 数据预处理
# 读取并过滤低表达基因 (至少在3个样本中FPKM > 1)，消除零方差警告
raw_data <- read.csv("Final_combined_FPKM.csv", row.names = 1)
keep <- rowSums(raw_data > 1) >= 3
expr_filtered <- log2(raw_data[keep, ] + 1)

# 3. 提取 F1 Hybrid 和 Sim 亲本数据 (以 Adult 为例)
# 筛选列名包含 "Adult" 且为 "Hybrid" 或 "Sim" 的样本
sel_cols <- grepl("Adult", colnames(expr_filtered)) & 
            (grepl("Hybrid", colnames(expr_filtered)) | (grepl("Sim", colnames(expr_filtered)) & !grepl("Hybrid", colnames(expr_filtered))))
expr_f1 <- expr_filtered[, sel_cols]

# 4. 构建分组与设计矩阵
samples <- colnames(expr_f1)
genotype <- ifelse(grepl("Hybrid", samples), "Hybrid", "Sim")
env <- gsub(".*_([A-Z]{2})\\d$", "\\1", samples) # 提取 CC, CO, OC, OO
group <- factor(paste(genotype, env, sep="_"))

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- lmFit(expr_f1, design)

# 5. 定义对比 (重点：非加性表达与环境特异性)
cont_matrix <- makeContrasts(
  # 杂种与亲本的差异 (非加性表达)
  F1_vs_P_CC = Hybrid_CC - Sim_CC,
  F1_vs_P_OO = Hybrid_OO - Sim_OO,
  # F1 特异性环境响应：(Hybrid_O - Hybrid_C) - (Sim_O - Sim_C)
  F1_Env_Specific = (Hybrid_OO - Hybrid_CC) - (Sim_OO - Sim_CC),
  levels = design
)

fit2 <- contrasts.fit(fit, cont_matrix)
fit2 <- eBayes(fit2)

# 6. 导出结果至 CSV
# 统计显著基因数量 (P < 0.05)
results <- decideTests(fit2, method="global", adjust.method="none", p.value=0.05)
write.csv(summary(results), "F1_Significant_Counts_Summary.csv")

# 导出 F1 环境特异性响应基因详细列表
f1_spec_list <- topTable(fit2, coef="F1_Env_Specific", n=Inf, p.value=0.05, adjust="none")
write.csv(f1_spec_list, "F1_Env_Specific_Genes_Detailed.csv")

# 7. 高清绘图输出 (4:3, 300 DPI, PNG格式)
img_w <- 8   # 英寸
img_h <- 6   # 英寸
img_res <- 300

# --- 图片 1: PCA 聚类图 ---
pca_data <- prcomp(t(expr_f1))
pca_df <- data.frame(pca_data$x, Genotype=genotype, Env=env)

png("F1_PCA_Analysis.png", width=img_w, height=img_h, units="in", res=img_res)
p1 <- ggplot(pca_df, aes(PC1, PC2, color=Genotype, shape=Env)) +
  geom_point(size=4, alpha=0.8) +
  theme_bw() +
  labs(title="PCA: F1 Hybrid vs Sim Parent", subtitle="Adult Samples") +
  theme(aspect.ratio=3/4)
print(p1)
dev.off()

# --- 图片 2: F1 特异性基因热图 ---
# 选取 F1_Env_Specific 中最显著的前 40 个基因
top_genes <- rownames(topTable(fit2, coef="F1_Env_Specific", n=40, adjust="none"))
plot_mat <- expr_f1[top_genes, ]



png("F1_Specific_Response_Heatmap.png", width=img_w, height=img_h, units="in", res=img_res)
pheatmap(plot_mat, 
         scale="row", 
         clustering_distance_cols="euclidean",
         show_colnames=TRUE, 
         main="Top 40 Genes: F1 Specific Environmental Response",
         fontsize=8)
dev.off()

cat("分析已完成！\n1. CSV 报告已生成。\n2. PNG 图片 (4:3, 300DPI) 已保存至您的工作目录。\n")