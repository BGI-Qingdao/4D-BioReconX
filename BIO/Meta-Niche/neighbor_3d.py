
from collections import Counter
import pandas as pd

meta = pd.read_csv('NsC_seurat_meta.deconv.txt', sep='\t', header=0)
meta['index'] = meta['label'].astype(int).astype(str) + '.' + meta['Batch']
meta = meta.set_index('index')

coords = pd.read_csv('rotation_pos.txt', sep=',', header=None, names=['Batch', 'label', 'x', 'y', 'z'])
coords['index'] = coords['label'].astype(int).astype(str) + '.' + coords['Batch']
coords = coords.set_index('index')

meta = meta.merge(coords, how='left', left_index=True, right_index=True)
meta = meta[meta['timepoint'] == 'WT']
meta['hood_inte.seurat_clusters'] = 'c' + meta['hood_inte.seurat_clusters'].astype(int).astype(str)

nb2_meta = meta[meta['hood_inte.seurat_clusters'] == 'c3']
nb2_arr = nb2_meta[['x', 'y', 'z']].values
other_meta = meta[meta['hood_inte.seurat_clusters'] != 'c3']
other_arr = other_meta[['x', 'y', 'z']].values

from sklearn.neighbors import NearestNeighbors
#stem cell distribution
nbrs = NearestNeighbors(n_neighbors=100, algorithm='ball_tree').fit(other_arr)
distances, indices = nbrs.radius_neighbors(nb2_arr, radius = 15)

all_types = []
all_coords = []
for index, (dis, ind) in enumerate(zip(distances, indices)):
    if len(ind) == 0 and len(dis) == 0:
        continue
    x, y, z = nb2_arr[index]
    nbrs_array = other_arr[ind] - [x, y, z]
    all_coords.extend(list(nbrs_array))
    
    types = other_meta.iloc[ind]['hood_inte.seurat_clusters'].values
    all_types.extend(list(types))

print(len(all_types))
print(len(all_coords))
niche = pd.DataFrame(all_coords, columns=['x', 'y', 'z'])
niche['hood_inte.seurat_clusters'] = all_types
niche.to_csv('test.txt', sep='\t')


