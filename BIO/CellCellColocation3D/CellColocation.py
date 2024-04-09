import sys
import getopt
import anndata as ad
import numpy as np
import pandas as pd
import igraph as ig
from scipy import stats
import matplotlib.pyplot as plt
import random

# keep result Repeatable
random.seed(42)


###########################################################
# grid functions 
###########################################################
#
# get grid points from one slice ( in sample area ) for probability sample
# 
def grid_in_slice(x, y, binsize):
    pos_data_x = x.copy()
    pos_data_y = y.copy()
    xmin = np.min(pos_data_x) - 2 * binsize
    xmax = np.max(pos_data_x) + 2 * binsize
    ymin = np.min(pos_data_y) - 2 * binsize
    ymax = np.max(pos_data_y) + 2 * binsize
    X, Y = np.mgrid[xmin:xmax:binsize, ymin:ymax:binsize]
    body_mask = np.zeros(X.shape)
    pos_data_x = pos_data_x - xmin
    pos_data_x = pos_data_x / binsize
    pos_data_x = pos_data_x.astype(int)
    pos_data_y = pos_data_y - ymin
    pos_data_y = pos_data_y / binsize
    pos_data_y = pos_data_y.astype(int)
    body_mask[pos_data_x, pos_data_y] = 1
    X_in_slice = X[body_mask == 1]
    Y_in_slice = Y[body_mask == 1]
    return X_in_slice, Y_in_slice


#
# get grid points from all slices ( in sample area ) for probability sample
# 
def get_3D_grid_sample_points(sample_raw_data):
    positions_for_density = []
    for sid in np.unique(sample_raw_data['z']):
        slice_data = sample_raw_data[sample_raw_data['z'] == sid]
        x_in_slice = slice_data['x'].to_numpy()
        y_in_slice = slice_data['y'].to_numpy()
        x_in_grid, y_in_grid = grid_in_slice(x_in_slice, y_in_slice, binsize)
        positions_for_density.append(np.vstack([x_in_grid, y_in_grid, (np.ones(x_in_grid.shape) * sid)]).T)
    positions = np.vstack(positions_for_density)
    return positions


###########################################################
# visualise functions 

###########################################################
# KDE, KL and MST
###########################################################

#
# get sampled density for all grid points
#
def kde3d_all(bootstraped_sample_data, positions, mincell):
    used_ct = []
    Pb_of_allct = []
    for ct in np.unique(bootstraped_sample_data['label']):
        ct_data = bootstraped_sample_data[bootstraped_sample_data['label'] == ct]
        cellnum = len(ct_data)
        if cellnum < mincell:
            continue
        used_ct.append(ct)
        all_points = np.vstack([ct_data['x'].to_numpy(), ct_data['y'].to_numpy(), ct_data['z'].to_numpy()])
        kde_ret = stats.gaussian_kde(all_points)
        Pb = kde_ret(positions.T)
        Pb = fill_zero(Pb)
        Pb_of_allct.append(Pb)
    return used_ct, Pb_of_allct


#
# fill 0 probability to avoid inf KL value
#
def fill_zero(Pb, esp=1e-20):
    Pb = Pb / np.sum(Pb)
    Pb[Pb < esp] = esp
    return Pb


#
# get KL diversity matrix
# 
def KL(used_ct, Pb_of_allct):
    kl_array = np.zeros((len(used_ct), len(used_ct)))
    for i, ct1 in enumerate(used_ct):
        for j, ct2 in enumerate(used_ct):
            if ct1 == ct2:
                continue
            # one row for one ref celltype
            kl_array[i, j] = stats.entropy(pk=Pb_of_allct[i], qk=Pb_of_allct[j], base=2)
    return kl_array


#
# get MST from KL matrix
#
def MSTdf_from_KL(used_ct, kl_array):
    ######################################################
    # get mst
    ######################################################
    yidx, xidx = np.nonzero(kl_array)
    edges = []
    weights = []
    for yid, xid in zip(yidx, xidx):
        edges.append([yid, xid])  # from KL ref to KL query
        weights.append(kl_array[yid, xid])
    ag = ig.Graph(n=len(used_ct), edges=edges, edge_attrs={'weight': weights}, directed=True)
    mst = ag.spanning_tree(weights=ag.es["weight"])
    mst_maskarray = np.array(mst.get_adjacency().data)
    ######################################################
    # gen result
    ######################################################
    yidx, xidx = np.nonzero(mst_maskarray)
    used_ct = np.array(used_ct, dtype='str')
    round_pd = pd.DataFrame()
    round_pd['from'] = used_ct[yidx.tolist()]
    round_pd['to'] = used_ct[xidx.tolist()]
    return round_pd


###########################################################
# cell type functions
###########################################################
def loading_data(in_h5ad_file, spatial_key, annotation_key):
    adata = ad.read_h5ad(in_h5ad_file)
    sample_raw_data = pd.DataFrame()
    sample_raw_data['x'] = adata.obsm[spatial_key][:, 0]
    sample_raw_data['y'] = adata.obsm[spatial_key][:, 1]
    sample_raw_data['z'] = adata.obsm[spatial_key][:, 2]
    sample_raw_data['label'] = adata.obs[annotation_key].to_list()

    total = len(sample_raw_data)
    cids = []
    proportions = []
    for ct in np.unique(sample_raw_data['label']):
        cids.append(ct)
        proportions.append(np.sum(sample_raw_data['label'] == ct))
    draw_points = pd.DataFrame()
    draw_points['cid'] = cids
    # draw_points['cid']= draw_points['cid'].astype(int)
    draw_points['fraction'] = proportions
    draw_points['fraction'] = draw_points['fraction'].astype(float)
    return sample_raw_data, draw_points


# main logic
###########################################################

def mainpipe(in_h5ad_file,
             prefix='prefix',
             bootnum=100,
             binsize=50,
             sample_fac=0.8,
             mincell=100,
             spatial_key='spatial',
             annotation_key='annotation'):
    # loading meta data
    sample_raw_data, drawpoints = loading_data(in_h5ad_file, spatial_key, annotation_key)
    # create points-in-sample for each slice
    positions = get_3D_grid_sample_points(sample_raw_data)
    # bootstrap
    temp_adjacent_list = []
    for rid in range(bootnum):
        # sample data
        bootstraped_sample_data = sample_raw_data.sample(frac=sample_fac)
        # kde density estimation
        used_ct, Pb_of_allct = kde3d_all(bootstraped_sample_data, positions, mincell)
        if len(used_ct) < 2:
            continue
        # KL divergency
        kl_array = KL(used_ct, Pb_of_allct)
        # save KL matrix
        kl_pd = pd.DataFrame(kl_array, index=used_ct, columns=used_ct)
        kl_pd.to_csv(f'{prefix}.KL.{rid}.txt', sep='\t')
        # MST
        round_pd = MSTdf_from_KL(used_ct, kl_array)
        round_pd['weight'] = [1.0 / float(bootnum)] * len(round_pd)
        temp_adjacent_list.append(round_pd)
    # save result
    final_net = pd.concat(temp_adjacent_list, ignore_index=True)
    final_net.to_csv(f'{prefix}.mst.all.csv', sep='\t', header=True, index=False)
    final_net = final_net.groupby(['from', 'to']).agg('sum').reset_index()
    final_net.to_csv(f'{prefix}.mst.merge.csv', sep='\t', header=True, index=False)


###########################################################
# entry, usage and parameter
###########################################################
def KL_usage():
    print("""
Usage: CellColocation.py -i <in.h5ad>
                         -o <prefix>
                         -b [binsize (for each slice), default 50 (refer to 50um if 1 unit=1um).]
                            Please set -b based on your coordinate system !!!
                            Notice: -b must be integer
                         -m [min cell number for a cell type, default 100.]
                         -f [sample fraction, default 0.8. The bootstrap fraction.]
                         -l [loop number, default 100. The number of iterations for bootstrap.]
                         -s [spatial key in obsm, default 'spatial']
                         -a [annotation key in obs, default 'annotation']
    """)


if __name__ == '__main__':
    # usage
    if len(sys.argv) == 1 or (len(sys.argv) == 2 and sys.argv[1] in ("-h", "--help")):
        KL_usage()
        exit(0)

    # default parameters
    mincell = 100  # ignore one celltype if total cell number less than 100
    bootnum = 100
    sample_fac = 0.8
    binsize = 50
    in_h5ad_file = ''
    prefix = ''
    spatial_key = 'spatial'
    anno_key = 'annotation'
    # parse inputs
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:o:b:m:f:l:a:s:", ["help"])
    except getopt.GetoptError:
        KL_usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            KL_usage()
            sys.exit(2)
        elif opt in ("-i"):
            in_h5ad_file = arg
        elif opt in ("-o"):
            prefix = arg
        elif opt in ("-b"):
            binsize = int(arg)
        elif opt in ("-m"):
            mincell = int(arg)
        elif opt in ('-f'):
            sample_fac = float(arg)
        elif opt in ('-l'):
            bootnum = int(arg)
        elif opt in ('-s'):
            spatial_key = arg
        elif opt in ('-a'):
            anno_key = arg

    # sanity check
    if in_h5ad_file == '' or prefix == '' or \
            binsize < 1 or mincell < 0 or bootnum < 1 or \
            sample_fac < 0.001 or sample_fac > 0.999:
        print('ERROR: invalid parameter!!!')
        KL_usage()
        sys.exit(1)
    # call main
    mainpipe(in_h5ad_file,
             prefix=prefix,
             bootnum=bootnum,
             binsize=binsize,
             sample_fac=sample_fac,
             mincell=mincell,
             spatial_key=spatial_key,
             annotation_key=anno_key)