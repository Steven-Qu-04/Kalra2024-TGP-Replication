# 1. 环境准备与文件夹创建
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("pheatmap")) install.packages("pheatmap")

if(!dir.exists("images")) dir.create("images")
basedir <- "E:/BaiduSyncdisk/试验设计/"
setwd(basedir)
# 2. 全局样本聚类图 (基于 euclidean_linkage.csv)
# 注意：聚类图通常需要原始距离矩阵，若只有 linkage 过程，
# 我们通过展示主要因子的 DE 数量来体现全局差异规模
dist_data <- read.csv("euclidean_linkage.csv")
# 模拟聚类树状分布的规模感 (由于 csv 为 linkage 过程，此处建议直接绘制 DE 总数)

# 3. 绘制差异表达基因 (DEGs) 汇总柱状图
# 读取三个核心结果文件并统计显著基因 (adj.P.Val < 0.05)
files <- c("Full_Results_Stage_Effect.csv", "Full_Results_Maternal_Effect.csv", "Full_Results_Inter_Mat_Cur.csv")
names <- c("Developmental_Stage", "Maternal_Effect", "Interaction")

de_counts <- data.frame()
for(i in 1:length(files)){
  data <- read.csv(files[i])
  # 统计显著上调和下调基因
  sig_up <- sum(data$adj.P.Val < 0.05 & data$logFC > 0, na.rm = TRUE)
  sig_down <- sum(data$adj.P.Val < 0.05 & data$logFC < 0, na.rm = TRUE)
  de_counts <- rbind(de_counts, data.frame(Factor=names[i], Direction="Up", Count=sig_up))
  de_counts <- rbind(de_counts, data.frame(Factor=names[i], Direction="Down", Count=sig_down))
}

p1 <- ggplot(de_counts, aes(x=Factor, y=Count, fill=Direction)) +
  geom_bar(stat="identity", position="dodge") +
  theme_minimal() +
  scale_fill_manual(values=c("Up"="#E41A1C", "Down"="#377EB8")) +
  labs(title="Differentially Expressed Genes (DEGs) Summary", y="Number of Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("images/1_DEGs_Summary_Barplot.png", p1, width=8, height=6, dpi=300)

# 4. 顺反机制分布堆叠图 (基于 Final_Factor_Mechanism_Distribution.csv)
dist_df <- read.csv("Final_Factor_Mechanism_Distribution.csv")

p2 <- ggplot(dist_df, aes(x=Factor, y=n, fill=Mechanism)) +
  geom_bar(stat="identity", position="stack") +
  theme_minimal() +
  scale_fill_brewer(palette="Set3") +
  labs(title="Regulatory Mechanism Distribution by Factor", x="Biological Factor", y="Gene Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("images/2_Mechanism_Stacked_Barplot.png", p2, width=10, height=7, dpi=300)

# 5. 生成各个因子的独立饼图
factors <- unique(dist_df$Factor)
for(f in factors){
  pie_data <- dist_df %>% filter(Factor == f)
  
  p_pie <- ggplot(pie_data, aes(x="", y=n, fill=Mechanism)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0) +
    theme_void() +
    labs(title=paste("Mechanism Share:", f)) +
    geom_text(aes(label = paste0(round(n/sum(n)*100), "%")), 
              position = position_stack(vjust = 0.5))
  
  ggsave(paste0("images/3_Pie_Chart_", f, ".png"), p_pie, width=6, height=6, dpi=300)
}

# 6. 交互作用显著基因的机制细分图 (基于 Gene_Level_Mechanism_Identification.csv)
gene_mech <- read.csv("Gene_Level_Mechanism_Identification.csv")
inter_mech <- gene_mech %>% filter(Factor == "Interaction")

if(nrow(inter_mech) > 0){
  p3 <- ggplot(inter_mech, aes(x=Mechanism, fill=Mechanism)) +
    geom_bar() +
    theme_minimal() +
    labs(title="Mechanism Categories of Interaction-Signficant Genes", y="Gene Count") +
    theme(legend.position="none")
  
  ggsave("images/4_Interaction_Mechanism_Detail.png", p3, width=8, height=6, dpi=300)
}

cat("所有图表已成功生成并保存至 'images' 文件夹中。\n")