
from collections import defaultdict
import numpy as np
import pandas as pd
import anndata
from scipy import sparse

def binning(adata, bin_size=100, dim='3d', use_value='sum', celltype=None):
    obs = adata.obs.copy()
    obs['x_bin'] = (obs.x / bin_size).astype(int)
    obs['y_bin'] = (obs.y / bin_size).astype(int)

    cols = ['x_bin', 'y_bin']
    if dim == '2d':
        cols.append('z')
    if celltype:
        if isinstance(celltype, str):
            celltype = [celltype]
        cols.extend(celltype)
    grouped_df = obs.groupby(cols)
    
    mtx = []
    obs_names = []
    metadata = {'cell_number': [], 'x': [], 'y': []}
    if dim == '2d':
        metadata['z'] = []
    if celltype:
        metadata[celltype] = []
    for name, group in grouped_df:
        cell_name = '_'.join([f'{x}' for x in name])
        obs_names.append(cell_name)
        
        metadata['cell_number'].append(group.shape[0])
        metadata['x'].append(name[0])
        metadata['y'].append(name[1])
        if dim == '2d':
            metadata['z'].append(name[2])
        if celltype:
            if dim == '3d':
                metadata[celltype].append(name[2])
            elif dim == '2d':
                metadata[celltype].append(name[3])

        if use_value == 'mean':
            data = adata[group.index].X.mean(axis=0, dtype=np.float64)
        elif use_value == 'sum':
            data = adata[group.index].X.sum(axis=0, dtype=np.int32)
        mtx.append(data)
    
    mtx = np.vstack(mtx)
    obs = pd.DataFrame(data=metadata, index=obs_names)
    var = adata.var.copy(deep=True)
    grouped_adata = anndata.AnnData(
            X=sparse.csr_matrix(mtx),
            obs=obs,
            var=var
            )
    
    for layer in adata.layers.keys():
        mtx = []
        for _, group in grouped_df:
            if use_value == 'mean':
                data = adata[group.index].layers[layer].mean(axis=0, dtype=np.float64)
            elif use_value == 'sum':
                data = adata[group.index].layers[layer].sum(axis=0, dtype=np.int32)
            mtx.append(data)
        mtx = np.vstack(mtx)
        grouped_adata.layers[layer] = sparse.csr_matrix(mtx)

    return grouped_adata

def _distance_calculate(orig_point, points):
    from sklearn.neighbors import NearestNeighbors
    neigh = NearestNeighbors(n_neighbors=1)
    neigh.fit(orig_point)

    neigh_dist, neigh_ind = neigh.kneighbors(points, 
            n_neighbors=1, return_distance=True)

    return neigh_dist.min(), neigh_dist.max()

def ncseg(adata, celltype=None, meta_nCell=20, min_nCell=5, use_value='sum', 
        n_neighbors=None, radius=None,):
    """
    Description
    -----------------------------------------------------
    Non-Contiguous segmentation for spatial coordinate based on pre-
    defined cell attributes and expected group member
    
    Parameters
    -----------------------------------------------------
    adata: anndata.AnnData
        anndata object to group, .X must be raw counts
    celltype: str
        column name in adata.obs which should be detect as 1st 
    meta_nCell: int
        number of cells to aggregate as a single SPC cell 
    min_nCell: int
        minimal number of cells in a SPC cell
    use_value: {'sum', 'mean'}, default sum 
        'sum' or 'mean' of counts should be used as SPC raw counts

    n_neighbors: deprecated
    radius: deprecated
    """

    def call_neighbor(coords_df, n_neighbors=None, radius=None, 
            meta_nCell=20, min_nCell=5, meta_prefix=None):

        from shapely.geometry import MultiPoint

        assert 'x' in coords_df.columns
        assert 'y' in coords_df.columns
        cols = ['x', 'y']
        if 'z' in coords_df.columns:
            cols.append('z')

        array = coords_df[cols].to_numpy()
        n_clusters = int(coords_df.shape[0] / meta_nCell)
        
        if radius or n_neighbors:
            # *** do NOT use this, raius not right
            # bug not fixed
            from sklearn.neighbors import kneighbors_graph, radius_neighbors_graph
            from sklearn.cluster import AgglomerativeClustering, DBSCAN
            if radius is None:
                graph = kneighbors_graph(array, n_neighbors=n_neighbors, 
                        mode='connectivity', metric='euclidean')
                distance = kneighbors_graph(array, n_neighbors=n_neighbors, 
                        mode='distance', metric='euclidean').toarray()
            else:
                graph = radius_neighbors_graph(array, radius=radius, 
                        mode='connectivity', metric='euclidean')
                distance = radius_neighbors_graph(array, radius=radius, 
                        mode='distance', metric='euclidean').toarray()
            clustering = AgglomerativeClustering(n_clusters=n_clusters, affinity='precomputed',
                    connectivity=graph, linkage='complete').fit(distance).labels_
            #clustering = DBSCAN(eps=3, metric='precomputed').fit(distance).labels_
            #unique, counts = np.unique(clustering.labels_, return_counts=True)
        else:        
            from scipy.cluster.hierarchy import linkage, fcluster
            from scipy.spatial.distance import pdist
            Y = pdist(array, metric='euclidean')
            Z = linkage(Y, 'complete')
            clustering = fcluster(Z, t=n_clusters, criterion='maxclust')

        unique, counts = np.unique(clustering, return_counts=True)

        indices_array = coords_df.index.values
        cls2labels = dict()
        for cls, number in zip(unique, counts):
            if number < min_nCell:
                continue
            
            if number < meta_nCell:
                adapted_meta_copies = 0
                adapted_meta_nCell = meta_nCell
            else:
                adapted_meta_copies = number // meta_nCell
                left_nCell = number % meta_nCell
                if left_nCell > 0:
                    left_nCell_each_copy = left_nCell // adapted_meta_copies + 1
                    if left_nCell_each_copy > int(meta_nCell * 0.5):
                        # (nCell of last copy < meta_nCell) 
                        # & (nCell of last copy > meta_nCell / 2)
                        adapted_meta_copies += 1
                        adapted_meta_nCell = meta_nCell
                    else:
                        # nCell of each copy > meta_nCell
                        # & nCell of last copy > nCell of previous copies
                        adapted_meta_copies = adapted_meta_copies
                        adapted_meta_nCell = meta_nCell + left_nCell_each_copy
                else:
                    adapted_meta_copies = adapted_meta_copies
                    adapted_meta_nCell = meta_nCell
            
            indices = indices_array[[i for i, c in enumerate(clustering) if c == cls]]
            assert len(indices) == number

            if adapted_meta_copies > 1:
                indices = [indices[i:i+adapted_meta_nCell] for i in 
                        range(0, len(indices), adapted_meta_nCell)]
                for copy_index, meta_indices in enumerate(indices):
                    if len(meta_indices) < min_nCell:
                        continue
                    meta_coords = coords_df.loc[meta_indices][cols].values
                    meta_centorid = [np.median(i) for i in zip(*meta_coords)]
                    min_dis, max_dis = _distance_calculate([meta_centorid], meta_coords)
                    cls2labels[f'{meta_prefix}_{cls}_{copy_index}'] = [meta_indices, meta_centorid, min_dis, max_dis]
            else:
                meta_coords = coords_df.loc[indices][cols].values
                meta_centorid = [np.median(i) for i in zip(*meta_coords)]
                min_dis, max_dis = _distance_calculate([meta_centorid], meta_coords)
                cls2labels[f'{meta_prefix}_{cls}'] = [indices, meta_centorid, min_dis, max_dis]
        return cls2labels

    obs = adata.obs.copy(deep=True)
    grouped_df = obs.groupby([celltype])
    
    neighborhood = {}
    for name, group in grouped_df:
        cls2labels = call_neighbor(group, n_neighbors=n_neighbors, radius=radius, 
                meta_nCell=meta_nCell, min_nCell=min_nCell, meta_prefix=name)
        filtered_nCell = (group.shape[0] - sum([len(x[0]) for x in cls2labels.values()])) / group.shape[0]
        print(f'\t{filtered_nCell} cells filtered for {name}', flush=True)

        neighborhood.update(cls2labels)
    print(neighborhood)

    mtx = []
    obs_names = []
    metadata = defaultdict(list)
    for hood_name, (labels, coords, min_rad, max_rad) in neighborhood.items():

        obs_names.append(hood_name)
            
        metadata[celltype].append(hood_name.split('_')[0])
        metadata['cell_number'].append(len(labels))
        metadata['hood'].append(','.join(labels))

        metadata['x'].append(coords[0])
        metadata['y'].append(coords[1])
        if len(coords) == 3:
            metadata['z'].append(coords[2])
        metadata['min_radius'].append(min_rad)
        metadata['max_radius'].append(max_rad)

        if use_value == 'mean':
            data = adata[labels].X.mean(axis=0, dtype=np.float64)
        elif use_value == 'sum':
            data = adata[labels].X.sum(axis=0, dtype=np.int32)
        mtx.append(data)

    mtx = np.vstack(mtx)
    obs = pd.DataFrame(data=metadata, index=obs_names)
    var = adata.var.copy(deep=True)
    grouped_adata = anndata.AnnData(
            X=sparse.csr_matrix(mtx),
            obs=obs,
            var=var
            )
    
    for layer in adata.layers.keys():
        mtx = []
        for _, (labels, _, _, _) in neighborhood.items():
            if use_value == 'mean':
                data = adata[labels].layers[layer].mean(axis=0, dtype=np.float64)
            elif use_value == 'sum':
                data = adata[labels].layers[layer].sum(axis=0, dtype=np.int32)
            mtx.append(data)
        mtx = np.vstack(mtx)
        grouped_adata.layers[layer] = sparse.csr_matrix(mtx)

    return grouped_adata

def ncseg_by_ref(adata, by=None, meta=None, use_value='sum'):
    """Non-Contiguous segmentation by a pre-defined group category
    
    adata: anndata object to group
    by: cellname to labels dict, labels must in obs_names of adata
    meta: additional meta to add to the result, must be 
          dataframe with keys in 'by' as index
    use_value: 'sum' or 'mean' of counts should be used as SPC raw counts
    """

    mtx = []
    obs_names = []
    metadata = defaultdict(list)
    for hood_name, labels in by.items():

        obs_names.append(hood_name)
        
        if use_value == 'mean':
            data = adata[labels].X.mean(axis=0, dtype=np.float64)
        elif use_value == 'sum':
            data = adata[labels].X.sum(axis=0, dtype=np.int32)
        mtx.append(data)

    mtx = np.vstack(mtx)
    obs = pd.DataFrame(data=metadata, index=obs_names)
    if meta is not None:
        obs = obs.merge(meta, how='left', left_index=True, right_index=True)
    var = adata.var.copy(deep=True)
    grouped_adata = anndata.AnnData(
            X=sparse.csr_matrix(mtx),
            obs=obs,
            var=var
            )
    
    for layer in adata.layers.keys():
        mtx = []
        for _, labels in by.items():
            if use_value == 'mean':
                data = adata[labels].layers[layer].mean(axis=0, dtype=np.float64)
            elif use_value == 'sum':
                data = adata[labels].layers[layer].sum(axis=0, dtype=np.int32)
            mtx.append(data)
        mtx = np.vstack(mtx)
        grouped_adata.layers[layer] = sparse.csr_matrix(mtx)

    return grouped_adata


