import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
import pandas as pd

np.random.seed(42)

# ==============================
# 1. 生成适合马氏距离的数据
# ==============================
n_samples = 20
n_features = 10

# 构造相关性强的协方差矩阵
# 对角线为 1，其余为 0.7 → 特征强相关
corr = np.full((n_features, n_features), 0.7)
np.fill_diagonal(corr, 1.0)

# Cholesky 分解
L = np.linalg.cholesky(corr)

# 生成标准正态数据
X_raw = np.random.randn(n_samples, n_features)

# 产生相关数据
X = X_raw @ L.T

print("Data shape:", X.shape)

# ==============================
# 2. 计算距离矩阵
# ==============================
# 欧氏距离
D_euclid = pdist(X, metric='euclidean')

# 马氏距离
cov = np.cov(X, rowvar=False)
cov_inv = np.linalg.pinv(cov)  # 伪逆即可
D_mahal = pdist(X, metric='mahalanobis', VI=cov_inv)

# ==============================
# 3. 聚类
# ==============================
Z_euclid = linkage(D_euclid, method='single')
Z_mahal = linkage(D_mahal, method='single')

# ==============================
# 4. 保存聚类树到文件
# ==============================
pd.DataFrame(Z_euclid, columns=['idx1', 'idx2', 'distance', 'sample_count']).to_csv("euclidean_linkage.csv", index=False)
pd.DataFrame(Z_mahal, columns=['idx1', 'idx2', 'distance', 'sample_count']).to_csv("mahalanobis_linkage.csv", index=False)
print("Linkage matrices saved.")

# ==============================
# 5. 绘图对比
# ==============================
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

dendrogram(Z_euclid, ax=axes[0], labels=np.arange(n_samples))
axes[0].set_title("Single Linkage - Euclidean Distance")
axes[0].set_ylabel("Distance")

dendrogram(Z_mahal, ax=axes[1], labels=np.arange(n_samples))
axes[1].set_title("Single Linkage - Mahalanobis Distance")
axes[1].set_ylabel("Distance")

plt.tight_layout()
plt.show()
