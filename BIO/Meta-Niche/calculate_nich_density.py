
import numpy as np
import pandas as pd
from scipy import stats

import matplotlib.pyplot as plt

niche = pd.read_csv('scaled_coords_2d.r15.txt', 
        sep='\t', header=0, index_col=0)

niche = niche.groupby('hood_inte.seurat_clusters').filter(lambda x: len(x) >= 300)

xmin = np.min(niche.x)
xmax = np.max(niche.x)
ymin = np.min(niche.y)
ymax = np.max(niche.y)
X, Y = np.mgrid[xmin:xmax:0.1, ymin:ymax:0.1]

sort_info_cluster = []
sort_info_count = []
sort_info_max = []
density_sort = []
for index, cluster in enumerate(niche['hood_inte.seurat_clusters'].unique()):
    data = niche[niche['hood_inte.seurat_clusters'] == cluster]
    count_ = data.shape[0]

    x = data.x.values
    y = data.y.values
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([x, y])
    kernel = stats.gaussian_kde(values)
    Z = np.reshape(kernel(positions).T, X.shape)
    max_ = Z.max()
    #border_values = Z[]
    density_sort.append((cluster, Z, count_, max_))
    sort_info_count.append(count_)
    sort_info_max.append(max_)
    sort_info_cluster.append(cluster)

density_sort = sorted(density_sort, key=lambda x: (x[2], x[3]), reverse=True)

sort_info_dis = []
min_, max_ = 75, 225
plt.style.use('dark_background')
fig, axes = plt.subplots(4, 7, figsize=(21, 12), sharex=True, sharey=True)
for index, (cluster, Z, count_, density_) in enumerate(density_sort):
    row = index // 7
    col = index % 7
    ax = axes[row, col]
    
    contours = ax.contour(
            X[min_:max_, min_:max_], 
            Y[min_:max_, min_:max_], 
            Z[min_:max_, min_:max_], 
            levels=[0.001], 
            colors='r', 
            linestyles='dashed'
            )

    _, _, _, _, _, d2 = contours.find_nearest_contour(
            0, 0, pixel=False
            )
    sort_info_dis.append(d2)
    
    ax.plot(0, 0, '+b')
    circle = plt.Circle((0, 0), 5, color='black', ls='--', 
            lw=0.5, fill=False)
    ax.add_patch(circle)
    ax.imshow(Z, extent=[-15, 15, -15, 15], origin='lower', 
            cmap='RdGy_r', alpha=1, vmin=0, vmax=0.0025)
            #cmap='RdGy_r', alpha=1)
    ax.text(10, 10, cluster)
    ax.text(5, 5, d2)

fig.savefig('niche_density.pdf')


data = pd.DataFrame({ 
    'cluster': sort_info_cluster,
    'count': sort_info_count, 
    'max_density': sort_info_max,
    'distance': sort_info_dis,
    })
data = data.sort_values(['count', 'distance', 'max_density'], ascending=[False, True, False])
data.to_csv('niche_density.txt', sep='\t')


