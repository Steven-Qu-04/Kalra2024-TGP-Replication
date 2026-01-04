# 1. 过滤低表达基因 (至少在 3 个样本中 FPKM > 1)
# 这一步能消除 Zero sample variances 警告，提高统计效能
keep <- rowSums(raw_data > 1) >= 3
expr_filtered <- expr_set[keep, ]

# 2. 重新拟合 (沿用之前的 design 和 cont.matrix)
fit_f <- lmFit(expr_filtered, design)
fit_f2 <- contrasts.fit(fit_f, cont.matrix)
fit_f2 <- eBayes(fit_f2)

# 3. 提取交互项显著的基因名 (P < 0.05)
inter_res <- topTable(fit_f2, coef="Inter_Mat_Cur", p.value=0.05, number=Inf, adjust.method="none")
inter_gene_names <- rownames(inter_res)

cat("过滤后交互项显著基因数量：", length(inter_gene_names), "\n")

# 4. 绘制热图查看这 24 个基因的特异性表现
if (!require("pheatmap")) install.packages("pheatmap")
library(pheatmap)

# 提取这些基因在不同环境组合下的平均表达量
# 重点对比 CC, CO, OC, OO
plot_data <- expr_filtered[inter_gene_names, ]



pheatmap(plot_data, 
         show_colnames = T, 
         show_rownames = F, 
         scale = "row", 
         cluster_cols = T,
         main = "Interaction Significant Genes (P < 0.05)")
        # 保存热图为 PNG 文件
        png("interaction_significant_genes.png", width = 800, height = 600)
        pheatmap(plot_data, 
                    show_colnames = T, 
                    show_rownames = F, 
                    scale = "row", 
                    cluster_cols = T,
                    main = "Interaction Significant Genes (P < 0.05)")
        dev.off()