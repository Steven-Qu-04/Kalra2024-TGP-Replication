library(limma)
basedir <- "E:/BaiduSyncdisk/试验设计/"
# 1. 准备数据
raw_data <- read.csv("Final_combined_FPKM.csv", row.names = 1)
expr_set <- log2(raw_data + 1)

# 2. 建立组合分组
samples <- colnames(expr_set)
groups <- factor(paste(
  ifelse(grepl("Adult", samples), "A", "L"),
  ifelse(grepl("Hybrid", samples), "H", "S"),
  substr(gsub(".*_([A-Z]{2})\\d$", "\\1", samples), 1, 1),
  substr(gsub(".*_([A-Z]{2})\\d$", "\\1", samples), 2, 2),
  sep = "_"
))

design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)
fit <- lmFit(expr_set, design)

# 3. 定义对比矩阵
cont.matrix <- makeContrasts(
  Stage_Effect = (L_S_C_C + L_S_O_O)/2 - (A_S_C_C + A_S_O_O)/2,
  Maternal_Effect = (A_S_O_C + A_S_O_O)/2 - (A_S_C_C + A_S_C_O)/2,
  Inter_Mat_Cur = (A_S_O_O - A_S_O_C) - (A_S_C_O - A_S_C_C),
  levels = design
)

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# 4. 统计结果并导出到 CSV
# 设定阈值
p_cutoff <- 0.05
lfc_cutoff <- 0

# 导出汇总统计表
results <- decideTests(fit2, method="global", adjust.method="none", p.value=p_cutoff, lfc=lfc_cutoff)
summary_table <- summary(results)
write.csv(summary_table, "DE_logic_summary.csv")
cat("--- 差异表达统计已保存至 DE_logic_summary.csv ---\n")
print(summary_table)

# 导出每个对比组的完整数据 (TopTable)
# 使用循环将 Stage_Effect, Maternal_Effect, Inter_Mat_Cur 分别存表
all_coefs <- colnames(cont.matrix)
lapply(all_coefs, function(x) {
  temp_tab <- topTable(fit2, coef = x, n = Inf, adjust.method = "none")
  write.csv(temp_tab, paste0("Full_Results_", x, ".csv"))
})

# 5. 绘制并保存韦恩图 (4:3, 300 DPI)
# 计算尺寸：4:3 比例下，若宽度为 8 英寸，高度则为 6 英寸
png(filename = "Venn_Diagram_Results.png", 
    width = 8, height = 6, units = "in", res = 300)

vennDiagram(results, 
            include = "both", 
            circle.col = c("#E41A1C", "#377EB8", "#4DAF4A"),
            main = "Overlap of Main Effects and Interactions\n(P < 0.05, Unadjusted)")

dev.off() # 关闭设备，保存文件
cat("--- 韦恩图已保存为 Venn_Diagram_Results.png (300 DPI) ---\n")

# 6. 提取交互项最显著的前 100 个基因
top_inter <- topTable(fit2, coef="Inter_Mat_Cur", number=100, sort.by="P", adjust.method = "none")
write.csv(top_inter, "Top_100_Interaction_Genes_Raw.csv")