#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date: Created on 28 Aug 2023 09:19
# @Author: Yao LI
# @File: clustering.py

import pandas as pd
from sklearn.linear_model import LogisticRegression
import hdbscan
import sys

from pcg_pattern import *
from scipy.ndimage import gaussian_filter1d
import numpy as np
from plot_heatmap import plot_heatmap
import scanpy as sc
from hdbscan import flat


'''
python clustering.py scaled_cell_norm.h5ad wt.h5ad
'''
# load total genes
scaled_fn = sys.argv[1]
# data = pd.read_csv(scaled_fn)
adata = sc.read_h5ad(scaled_fn)
data = adata.to_df()
data = sort_bin(data)

print(data)
print(data.shape)
gene_remains = filter_genes(data, 5)
sub_data = data[gene_remains]
print(sub_data.shape)
#sub_data = smooth(sub_data, 3)
new_norm1 = norm_scale(sub_data)

# clustering based on WT data
clusterer = hdbscan.HDBSCAN()
clusterer.fit(new_norm1.T.to_numpy())
labels = clusterer.labels_
labels = labels.astype(int)
used_labels = labels.copy()

# original WT labels
np.savetxt('1.wt_labels.txt', labels, fmt='%i', delimiter=',')
# cluster--gene
cluster_genes = cal_cluster_genes(new_norm1.T, labels)
# non_cluster genes
l = cluster_genes[-1]
with open('1.wt_no_cluster_genes.txt', 'w') as f:
    f.writelines('\n'.join(l))

# sort clusters
corder, non_cluster = sort_cluster(cluster_genes, sub_data)
# update original WT labels
labels = update_sorted_labels(corder, labels)
np.savetxt('1.wt_labels_sorted.txt', labels, fmt='%i', delimiter=',')
cg = cal_cluster_genes(new_norm1.T, labels)
save_cluster_gene(cg, '1.wt_cg_sorted.csv', '.')
new_gc = cal_gene_cluster(list(new_norm1.columns), labels)
save_gene_cluster(new_gc, '1.wt_gc_sorted.txt', '.')

# vi
gene_mat = new_norm1.T.to_numpy()
# draw heatmap
new_order = list(range(max(labels) + 1))
# draw_heatmap(gene_mat, labels, '1.wt', '.', new_order)
plot_heatmap(gene_mat, labels, '1.wt', new_order)
# draw umap
draw_umap(gene_mat, labels, non_cluster, '.', '1.wt_sorted', new_order)


# 2.
#load total data: contains all samples
save = '.'
wt_data = new_norm1.T
non_data = wt_data[used_labels==-1]
rest_data = wt_data[used_labels!=-1]
non_genes = list(non_data.index)
print(f'number of non genes {len(non_genes)}')


svc = LogisticRegression()
svc.fit(rest_data.to_numpy(),used_labels[used_labels!=-1])
recall_non_labels = svc.predict(non_data.to_numpy())

results = pd.DataFrame({'gene':non_genes, 'cluster':recall_non_labels.tolist()})
results.to_csv('2.labels_lr.txt', index=False, sep='\t')


# 4. projection
WT_data = new_norm1
wt_genes = list(WT_data.columns)

data_folder = sys.argv[2]
samples = ["0hpa1", "0hpa2","12hpa1", "12hpa2", "36hpa1", "36hpa2", "3dpa1", "3dpa2", 
"5dpa1", "5dpa2", "7dpa1", "7dpa2", "10dpa1", "10dpa2", "14dpa1", "14dpa2", "WT"]
sample_header = {}

# orignal WT labels
labels = np.loadtxt(sys.argv[3])
# update labels based on recall output
recall = pd.read_csv(sys.argv[4], delimiter='\t')
labels = update_recall_labels(labels, recall, wt_genes)
labels = labels.astype(int)
new_order = list(range(max(labels)+1))
draw_heatmap(WT_data.T.to_numpy(), labels, '3.recalled_WT', save, new_order)


# plot only WT data points
cluster_gene = cal_cluster_genes(WT_data.T, labels)
corder,non_cluster = sort_cluster(cluster_gene, WT_data)
labels = update_sorted_labels(corder, labels)
print(f'WT label size is {labels.shape}')
print(labels.max(), labels.min())
np.savetxt('3.labels_sorted2_wt.txt', labels, fmt='%i', delimiter=',')

new_order = list(range(max(labels)+1))
draw_heatmap(WT_data.T.to_numpy(), labels, '3.projected_WT', save, new_order)
draw_umap(WT_data.T.to_numpy(), labels, -1, save, '3.sorted_WT')


#--------------------------------#
#    ADD the rest of datasets    #
#--------------------------------#
# save all sample-label data in sample_header
n_wt_genes = [f'{x}_WT' for x in wt_genes]
sample_header['WT'] = pd.DataFrame({'gene':n_wt_genes, 'labels':labels.tolist()})
total_gene_df = []
total_gene_df.append(WT_data)

# projection
svc = LogisticRegression()
svc.fit(WT_data.T.to_numpy(),labels)
for sp in samples:
    if sp == 'WT':
        continue
    sample_fn = os.path.join(data_folder, f'{sp}.HVG5000.log_mean.h5ad')
    sample_data = sc.read_h5ad(sample_fn)
    data = sample_data.to_df()
    gene_remains = filter_genes(data, 5)
    sub_data = data[gene_remains]
    new_norm1 = norm_scale(sub_data)
    data_df = new_norm1
    total_gene_df.append(data_df)
    genes = list(data_df.columns)
    n_genes = [f'{x}_{sp}' for x in genes]
    print(f'sample {sp}: genes {len(genes)}')
    labels = svc.predict(data_df.T.to_numpy())
    sample_header[sp] = pd.DataFrame({'gene':n_genes, 'labels':labels.tolist()})
 
# save final results
results = []
for sp in samples:
    results.append(sample_header[sp])

current_result = pd.concat(results, ignore_index=True, sort=False)  ## !!!! sort=False, highly crucial
current_result.to_csv('3.projected_gene_cluster.csv',index=None,sep=',')


#--------------------------------------------#
#      Update latest clustering result       #
#--------------------------------------------#
# get the latest label
labels = current_result['labels'].to_numpy()
np.savetxt('3.labels_projected_final_all.txt', labels, fmt='%i', delimiter=',')
print(len(labels))
print(len(list(data.columns)))

# draw·
new_order = list(range(max(labels)+1))
# create total gene
total_gene_data = pd.concat(total_gene_df, axis=1)
gene_mat = total_gene_data.T.to_numpy()
#draw_heatmap(gene_mat, labels, '3.projected_all', save, new_order)
plot_heatmap(gene_mat, labels, '3.projected_all', new_order)
draw_umap(gene_mat, labels, -1, save, '3.sorted_all')

