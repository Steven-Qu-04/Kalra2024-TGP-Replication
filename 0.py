import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform

# =========================================================
# PARAMETER INTERFACE
# =========================================================
LINKAGE_METHOD = "complete"   # "single" or "complete"
OUTPUT_FILE = f"{LINKAGE_METHOD}_linkage_demo.mp4"

# =========================================================
# 1. Generate data
# =========================================================
np.random.seed(42)
X = np.random.randn(10, 2)
n = len(X)

# =========================================================
# 2. Hierarchical clustering
# =========================================================
Z = linkage(X, method=LINKAGE_METHOD)
D = squareform(pdist(X))

# =========================================================
# 3. Figure
# =========================================================
fig = plt.figure(figsize=(12, 5))
ax_scatter = fig.add_subplot(121)
ax_dendro = fig.add_subplot(122)

# =========================================================
# 4. Scatter base
# =========================================================
ax_scatter.scatter(X[:, 0], X[:, 1], c="lightgray", s=60, zorder=1)

for i, p in enumerate(X):
    ax_scatter.text(p[0] + 0.03, p[1] + 0.03, str(i), fontsize=9)

ax_scatter.set_title(f"{LINKAGE_METHOD.capitalize()} Linkage: Active Merge")

# dynamic artists
line_artist, = ax_scatter.plot([], [], "r-", lw=2, zorder=3)
clusterA_artist = ax_scatter.scatter([], [], c="dodgerblue", s=100, zorder=2)
clusterB_artist = ax_scatter.scatter([], [], c="limegreen", s=100, zorder=2)

# =========================================================
# 5. Helper functions
# =========================================================
def get_clusters(step):
    clusters = {i: [i] for i in range(n)}
    for k in range(step):
        a, b = int(Z[k, 0]), int(Z[k, 1])
        new_id = n + k
        clusters[new_id] = clusters[a] + clusters[b]
        clusters.pop(a)
        clusters.pop(b)
    return clusters


def extreme_pair(c1, c2):
    """
    Return nearest or farthest pair depending on linkage
    """
    best_val = np.inf if LINKAGE_METHOD == "single" else -np.inf
    pair = (None, None)

    for i in c1:
        for j in c2:
            d = D[i, j]
            if (LINKAGE_METHOD == "single" and d < best_val) or \
               (LINKAGE_METHOD == "complete" and d > best_val):
                best_val = d
                pair = (i, j)

    return pair

# =========================================================
# 6. Animation
# =========================================================
def animate(step):
    # ---------- dendrogram ----------
    ax_dendro.clear()
    dendrogram(
        Z,
        labels=[str(i) for i in range(n)],
        color_threshold=Z[step, 2],
        ax=ax_dendro
    )
    ax_dendro.set_title(
        f"Hierarchical Clustering ({LINKAGE_METHOD.capitalize()} Linkage)"
    )
    ax_dendro.set_ylabel("Distance")

    # ---------- cluster reconstruction ----------
    clusters = get_clusters(step)
    a, b = int(Z[step, 0]), int(Z[step, 1])

    clusterA = clusters[a]
    clusterB = clusters[b]

    # ---------- highlight clusters ----------
    clusterA_artist.set_offsets(X[clusterA])
    clusterB_artist.set_offsets(X[clusterB])

    # ---------- connection line ----------
    i, j = extreme_pair(clusterA, clusterB)
    line_artist.set_data(
        [X[i, 0], X[j, 0]],
        [X[i, 1], X[j, 1]]
    )

    return line_artist, clusterA_artist, clusterB_artist

# =========================================================
# 7. Create animation
# =========================================================
ani = FuncAnimation(
    fig,
    animate,
    frames=len(Z),
    interval=1200,
    repeat=False
)

# =========================================================
# 8. Save MP4 (requires ffmpeg)
# =========================================================
ani.save(
    OUTPUT_FILE,
    writer="ffmpeg",
    fps=1
)

plt.tight_layout()
plt.show()
