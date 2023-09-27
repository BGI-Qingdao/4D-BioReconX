
import pandas as pd
from sklearn.neighbors import NearestNeighbors

def cal_niche(df, group, group_col='celltype', radius=15):
    """calculate niche component

    Parameters
    ---------------------------
    df: pd.DataFrame
        dataframe contain coordinates and group information
    group: str
        group used as target
    group_col: str
        column name used as group information
    radius: int
        radius for neighbor calculation
    
    return: pd.DataFrame
        dataframe contain reformed niche information
    """

    coords_cols = ['x', 'y']
    if 'z' in df.columns:
        coords_cols += ['z']

    meta_cols = [col for col in df.columns if col not in coords_cols]

    target_df = df[df[group_col] == group]
    target_arr = target_df[coords_cols].values
    other_df = df[df[group_col] != group]
    other_arr = other_df[coords_cols].values

    nbrs = NearestNeighbors(n_neighbors=100, algorithm='ball_tree').fit(other_arr)
    distances, indices = nbrs.radius_neighbors(target_arr, radius=radius)

    nbrs_coords = []
    nbrs_values = []
    for index, (dis, ind) in enumerate(zip(distances, indices)):
        if len(ind) == 0 and len(dis) == 0:
            continue
        if 'z' in coords_cols:
            x, y, z = target_arr[index]
            nbrs_arr = other_arr[ind] - [x, y, z]
        else:
            x, y = target_arr[index]
            nbrs_arr = other_arr[ind] - [x, y]
        nbrs_coords.extend(list(nbrs_arr))
    
        meta_vals = other_df.iloc[ind][meta_cols].values
        nbrs_values.extend(list(meta_vals))
    
    niche = pd.DataFrame(nbrs_coords, columns=coords_cols)
    niche[meta_cols] = nbrs_values
    return niche

def cal_component(niche, group_col):
    df = niche.groupby(group_col).size().to_frame('counts')
    df['percentage'] = df['counts'] / df['counts'].sum()
    return df

def cal_2d_density(df, group_col, groups=None, min_num=50, save_figure='out.pdf'):
    """calculate density for each group

    Parameters
    ---------------------------
    df: pd.DataFrame
        dataframe contain reformed niche information
    group_col: str
        column name used as group information
    min_num: int
        minimal number used for calculcation
    save_figure: str
        figure file name
    """
    import numpy as np
    from scipy import stats
    import matplotlib.pyplot as plt

    df = df.groupby(group_col).filter(lambda x: len(x) >= min_num)

    xmin, xmax = df.x.min(), df.x.max()
    ymin, ymax = df.y.min(), df.y.max()
    X, Y = np.mgrid[xmin:xmax:0.1, ymin:ymax:0.1]

    
    data = []
    for name, group in df.groupby(group_col):

        _count = group.shape[0]
        x, y = group.x.values, group.y.values

        pos = np.vstack([X.ravel(), Y.ravel()])
        val = np.vstack([x, y])
        kernel = stats.gaussian_kde(val)
        Z = np.reshape(kernel(pos).T, X.shape)
        _max = Z.max()

        data.append([name, Z, _count, _max])
    data = sorted(data, key=lambda x: (x[2], x[3]), reverse=True)
    
    if groups is not None:
        data = [x for x in data if x[0] in groups]
    
    _min_lim, _max_lim = 75, 225
    ncols = 6
    nrows = len(data) / ncols

    plt.style.use('dark_background')
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)
    axes = axes.flatten()

    results = []
    for ax, (name, Z, _count, _max) in zip(axes, data):
        contours = ax.contour(
            X[_min_lim:_max_lim, _min_lim:_max_lim],
            Y[_min_lim:_max_lim, _min_lim:_max_lim],
            Z[_min_lim:_max_lim, _min_lim:_max_lim],
            levels=[0.001],
            colors='r',
            linestyles='dashed'
        )

        *d2 = contours.find_nearest_contour(0, 0, pixel=False)

        ax.plot(0, 0, '+b')
        circle = plt.Circle((0, 0), 5, color='k', ls='--', lw=0.5, fill=False)
        ax.add_patch(circle)
        ax.imshow(Z, entent=[-15, 15, -15, 15], origin='lower', 
                  cmap='RdGy_r', alpha=1, vmin=0, vmax=0.0025)
        ax.text(10, 10, name)
        ax.text(5, 5, ds)

        results.append([name, _count, _max, d2])
    fig.savefig(save_figure)

    colnames = ['group', 'count', 'max_dens', 'nearest_dist']
    results = pd.DataFrame(results, columns=colnames)
    results = results.sort_values(['count', 'nearest_dist', 'max_dens'], 
                    ascending=[False, True, False])
    return results

