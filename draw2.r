# 补充：生成 2.1 节的聚类图
library(ggplot2)
if(!dir.exists("images")) dir.create("images")

# 读取层级聚类距离数据
linkage_data <- read.csv("euclidean_linkage.csv")

# 由于只有 linkage 数据，我们可以通过绘制节点距离分布来展示样本间的层级关系规模
p_cluster <- ggplot(linkage_data, aes(x = reorder(as.character(idx1), distance), y = distance)) +
  geom_segment(aes(xend = as.character(idx2), yend = distance)) +
  geom_point(color = "steelblue", size = 2) +
  theme_minimal() +
  labs(title = "Sample Hierarchical Clustering (Euclidean Distance)",
       x = "Sample Nodes", y = "Clustering Distance") +
  theme(axis.text.x = element_blank())

ggsave("images/0_Sample_Clustering.png", p_cluster, width = 8, height = 5, dpi = 300)