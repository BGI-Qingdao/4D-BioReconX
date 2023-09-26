
from collections import Counter
import pandas as pd

meta = pd.read_csv('NsC_seurat_meta.deconv.txt', sep='\t', header=0)
meta['index'] = meta['label'].astype(int).astype(str) + '.' + meta['Batch']
meta = meta.set_index('index')
meta = meta[meta['timepoint'] == 'WT']

coords = pd.read_csv('WT_stat.txt', sep='\t', header=0)
coords['index'] = coords['label'].astype(int).astype(str) + '.' + coords['Batch']
coords = coords.set_index('index')
coords = coords.rename(columns={'cx': 'x', 'cy': 'y'})
coords = coords[['x', 'y']]

meta = meta.merge(coords, how='left', left_index=True, right_index=True)
meta['hood_inte.seurat_clusters'] = 'c' + meta['hood_inte.seurat_clusters'].astype(int).astype(str)

from sklearn.neighbors import NearestNeighbors

all_types = []
all_coords = []
all_dis = []
for section in meta['Batch'].unique():
    data = meta[meta['Batch'] == section]

    nb2_meta = data[data['hood_inte.seurat_clusters'] == 'c3']
    nb2_arr = nb2_meta[['x', 'y']].values
    other_meta = data[data['hood_inte.seurat_clusters'] != 'c3']
    other_arr = other_meta[['x', 'y']].values

    #stem cell distribution
    nbrs = NearestNeighbors(n_neighbors=100, algorithm='ball_tree').fit(other_arr)
    distances, indices = nbrs.radius_neighbors(nb2_arr, radius = 15)

    for index, (dis, ind) in enumerate(zip(distances, indices)):
        if len(ind) == 0 and len(dis) == 0:
            continue
        x, y = nb2_arr[index]
        nbrs_array = other_arr[ind] - [x, y]
        all_coords.extend(list(nbrs_array))

        all_dis.extend(list(dis))
        
        types = other_meta.iloc[ind]['hood_inte.seurat_clusters'].values
        all_types.extend(list(types))

niche = pd.DataFrame(all_coords, columns=['x', 'y'])
niche['dist'] = all_dis
niche['hood_inte.seurat_clusters'] = all_types
niche.to_csv('scaled_coords_2d.r15.txt', sep='\t')


import matplotlib.pyplot as plt
import sys
sys.path.append('/dellfsqd2/ST_OCEAN/USER/hankai/Project/10.smed/code/site-packages')
from WBRtools import dark_colors

color_map = dict((f'c{i}', c) for i, c in enumerate(dark_colors))

plt.style.use('dark_background')
fig, ax = plt.subplots(figsize=(10, 10))
x = niche.x.values
y = niche.y.values
#s = niche['dist'].max() - niche['dist']
#s = niche['dist']
c = [color_map[i] for i in niche['hood_inte.seurat_clusters'].values]

ax.scatter(x, y, s=5, c=c, 
        alpha=0.5,
        linewidths=0, edgecolors=None,
        rasterized=True,
        )

ax.scatter(0, 0, s=100, c='w', marker='+',
        linewidths=0, edgecolors=None,
        )

ax.set_yticks([])
ax.set_ylabel('')
_all = ['left', 'right', 'top']
for spine in _all:
    ax.spines[spine].set_visible(False)

fig.savefig('niche.pdf', dpi=600)


