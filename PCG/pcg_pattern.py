#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date: Created on 04 Apr 2022 11:02
# @Author: Yao LI


import argparse
import sys
#2023-07-21:
sys.path.append('/dellfsqd2/ST_OCEAN/USER/liyao1/03.planarian/pcg_ap_pattern')
#sys.path.append('/dellfsqd2/ST_OCEAN/USER/liyao1/scripts/pcg_ap_pattern')
import pandas as pd
import numpy as np
import hdbscan
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import umap
import umap.plot
import seaborn as sns
import os
from itertools import chain
from hdbscan import flat
from stereo_pcg.util import read_list, read_expr
from scipy.ndimage import gaussian_filter1d
from sklearn.decomposition import PCA
import scipy.stats
import time


def logger(smt :str):
    print("LOG : [{}] at time [{}]".format(smt,time.asctime(time.localtime(time.time()))),flush = True)


#-------------------------------------------------------
#          Pre-process data
#-------------------------------------------------------
def slide_window(gene, window_size:int) -> bool:
    for i in range(len(gene)-window_size+1):
        profile = gene[i:i+window_size]
        if not profile.isin([0]).any().any():
            return True


def filter_genes(data:pd.DataFrame, window_size: int=5) -> list:
    filter_genes = []
    for gene_id in data.columns:
        gene = data[gene_id]
        keep_flag = slide_window(gene, window_size)
        if keep_flag:
            filter_genes.append(gene_id)
    return filter_genes


def norm(sub_data:pd.DataFrame) -> pd.DataFrame:
    d = {}
    for g in sub_data.columns:
        norm_g = (sub_data[g]-np.mean(sub_data[g]))/np.std(sub_data[g])
        d[g] = norm_g
    new_norm = pd.DataFrame(d)
    return new_norm


def norm_scale(sub_data: pd.DataFrame, new_max: int=3, new_min: int=-3) -> pd.DataFrame:
    d = {}
    for g in sub_data.columns:
        target_range = new_max - new_min
        data_range = sub_data[g].max() - sub_data[g].min()
        scale_factor = target_range / data_range
        new_data = (sub_data[g] - sub_data[g].min()) * scale_factor + new_min
        d[g] = new_data
    new_norm = pd.DataFrame(d)
    return new_norm


def smooth(data: pd.DataFrame, sig:int=3):
    new_data = data.copy()
    for gene in data.columns:
        new = gaussian_filter1d(data[gene],sigma=sig)
        new_data[gene] = new
    return new_data


def pca_2d(new_norm):
    pca = PCA(n_components=2)
    p_data = pca.fit_transform(new_norm.to_numpy().T)
    #eigen_value = pca.explained_variance_
    #eigen_vector = pca.components_
    ratio = pca.explained_variance_ratio_
    print(f'first three are {ratio[0:3]}')
    return p_data


#--------------------------------------------
#       Clustering
#--------------------------------------------
def doUMAP(X:np.ndarray, labels:np.ndarray) -> np.ndarray :
    #X = X[labels!=-1]
    mapper = umap.UMAP().fit(X)
    #umap.plot.points(mapper, labels=labels[labels!=-1])
    X2 = mapper.transform(X)
    return X2
    #return X2, mapper


def doHDBSCAN(X :np.ndarray, min_num: int, epsilon, minsample:int) -> np.array:
    if min_num==0 and epsilon==0 and minsample==0:
        print('default param')
        clusterer = hdbscan.HDBSCAN()
    elif min_num > 0 and epsilon==0 and minsample==0:
        clusterer = hdbscan.HDBSCAN(min_cluster_size=min_num)
    else:
        clusterer = hdbscan.HDBSCAN(min_cluster_size=min_num, cluster_selection_epsilon=epsilon, min_samples=minsample)
    clusterer.fit(X)
    return clusterer


#--------------------------
#
#------------------
def draw_cluster(data_2d, clusterer, sample: str, path: str):
    color_palette = sns.color_palette('Paired', clusterer.labels_.max()+1)
    cluster_colors = [color_palette[x] if x >= 0
                      else (0.5, 0.5, 0.5)
                      for x in clusterer.labels_]
    cluster_member_colors = [sns.desaturate(x, p) for x, p in
                             zip(cluster_colors, clusterer.probabilities_)]
    plt.scatter(*data_2d.T, s=50, linewidth=0, c=cluster_member_colors, alpha=0.5)
    plt.savefig(os.path.join(path,f"{sample}_{clusterer.labels_.max()+1}clusters_umap.png"))
    plt.close()


def drawHDBSCAN_by_umapplot(mapper,labels,prefix):
    umap.plot.points(mapper, labels=labels)
    plt.legend(loc=(1.04,0))
    plt.savefig("{}.show_hdbscan_by_umapplot.png".format(prefix))
    plt.close()


def draw_heatmap(gene_mat, labels:np.ndarray, sample: str, path, orders: list=None):
    vmin = np.min(gene_mat)
    vmax = np.max(gene_mat)
    cluster_num = labels.max()+1
    if orders is None:
        fig, axs = plt.subplots(nrows=cluster_num+1, ncols=1,figsize=(6, 12))
        
        # don't draw -1 cluster in the heatmap
        #for i in range(cluster_num):
        for i in range(cluster_num-1):
            sns.heatmap(gene_mat[np.where(labels==i)], ax=axs[i], vmin=vmin, vmax=vmax, cmap="twilight_shifted", cbar=False)
            axs[i].set_xticks([])
            axs[i].set_yticks([])
        #sns.heatmap(gene_mat[np.where(labels==-1)], ax=axs[-1], vmin=vmin, vmax=vmax, cmap="twilight_shifted", cbar=False)
        #axs[-1].set_yticks([])
        axs[-1].set_yticks([])
        plt.savefig(os.path.join(path,f"{sample}_{cluster_num}clusters_heatmap.png"))
        plt.close()
    else:
        fig, axs = plt.subplots(nrows=len(orders)+1, ncols=1,figsize=(6, 12))
        axs[0].arrow(-1.5,5,12,0, width=0.05, length_includes_head=True, head_length=1, shape='right', color='black')
        axs[0].text(-1.5, 5.03, 'A', size=20, weight='bold', ha='left', va='baseline')
        axs[0].text(8.5, 5.03, 'P', size=20, weight='bold', ha='left', va='baseline')
        axs[0].spines['top'].set_visible(False)
        axs[0].spines['right'].set_visible(False)
        axs[0].spines['bottom'].set_visible(False)
        axs[0].spines['left'].set_visible(False)
        axs[0].set_xticks([])
        axs[0].set_yticks([])
        for i in range(len(orders)-1):
        #for i in range(len(orders)):
            sns.heatmap(gene_mat[np.where(labels==orders[i])], ax=axs[i+1], vmin=vmin, vmax=vmax, cmap="twilight_shifted", cbar=False)
            axs[i+1].set_xticks([])
            axs[i+1].set_yticks([])
        sns.heatmap(gene_mat[np.where(labels==orders[-1])], ax=axs[-1], vmin=vmin, vmax=vmax, cmap="twilight_shifted", cbar=False)
        axs[-1].set_yticks([])
        plt.savefig(os.path.join(path,f"{sample}_{len(orders)-1}clusters_ordered_heatmap.png"))
        plt.close()


def cal_gene_cluster(genes: list, labels: np.ndarray) -> dict:
    gene_label = {}
    for i in range(len(genes)):
        if genes[i] not in gene_label:
            gene_label[genes[i]] = int(labels[i])
        else:
            print("! this gene occured before !")
    return gene_label


def cal_cluster_genes(gene_mat_df: pd.DataFrame, labels: np.ndarray) -> dict:
    clusters = list(set(sorted(labels.tolist())))
    cluster_genes = {}
    for cid in clusters:
        cid = int(cid)
        cluster_genes[cid] = gene_mat_df.iloc[np.where(labels == cid)[0].tolist()].index.tolist()
    return cluster_genes


def get_gene(cid: str, cluster_genes: dict):
    try:
        return cluster_genes[cid]
    except KeyError:
        return "not available"


def get_cluster(gene: str, gene_label: dict):
    try :
        return gene_label[gene]
    except KeyError:
        return "not available"


#---------------------------------------------------#
#        Re-clustering of unassigned genes          #
#---------------------------------------------------#
def is_cluster_pearson(cid: int, nog: str, gene_mat_df, cluster_genes:dict, cutoff) -> bool:
    """Decide if a non-assigned gene has pcc higher than 0.8 with at least one gene in the cluster"""
    c_gene = get_gene(cid, cluster_genes)
    #is_cluster = False
    pccs = []
    for g in c_gene:
        pcc = np.corrcoef(gene_mat_df.loc[[g]], gene_mat_df.loc[[nog]])[0,1]
        pccs.append(pcc)
    if max(pccs) >= cutoff:
        return True, max(pccs)
    else:
        return False, max(pccs)
        #if pcc >= cutoff:
        #    is_cluster = True
        #    return is_cluster, pcc
    #return is_cluster, 0


def is_cluster_spearman(cid: int, nog: str, gene_mat_df: pd.DataFrame, cluster_genes: dict, cutoff) -> bool:
    """Decide if a non-assigned gene has spearman correlation coefficient value higher than cutoff with at least one gene in the cluster"""
    c_gene = get_gene(cid, cluster_genes)
    sccs = []
    for g in c_gene:
        scc = scipy.stats.spearmanr(gene_mat_df.loc[[g]], gene_mat_df.loc[[nog]]).correlation
        sccs.append(scc)
    if max(sccs) >= cutoff:
        return True, max(sccs), sum(1 for i in sccs if i >= cutoff)
    else:
        return False, max(sccs), 0


def re_clustering(gene_mat_df:pd.DataFrame, labels: np.ndarray, cluster_genes:dict, cmethod:str, cutoff) ->np.ndarray:
    """Hard clustering
    reassign gene to the cluster with higher mean expression level
    """
    labels = labels.copy()
    clusters = list(set(sorted(labels.tolist())))
    no_cluster_genes = get_gene(-1, cluster_genes)
    n = 0
    for nog in no_cluster_genes:
        if (n%100) == 0:
            print(f'reclustering {no_cluster_genes.index(nog)+1}th/{len(no_cluster_genes)} gene, {round((no_cluster_genes.index(nog)+1)/len(no_cluster_genes), 2)*100}%...')
        kept_genes = {}
        for cid in clusters[0:-1]:
            if cmethod =='pearson':
                isc, max_cor = is_cluster_pearson(cid, nog, gene_mat_df, cluster_genes, cutoff)
                # if the gene might re-assign to this cluster
                if isc:
                    kept_genes[max_cor] = cid
            #elif cmethod=='spearman':
                #isc, max_cor, num = is_cluster_spearman(cid, nog, gene_mat_df, cluster_genes, cutoff)

        # if any cluster
        if len(kept_genes)>0:
            # get gene index among all the genes
            nog_index = gene_mat_df.index.tolist().index(nog)
            # find the gene in the original label and change its cluster id to the cluster has the highest pcc value
            labels[nog_index] = kept_genes[max(kept_genes.keys())]  
        n += 1
    return labels


def sort_cluster(cg:dict, expr:pd.DataFrame)->list:
    """Sort clusters by the position of their peak value, from left to right"""
    order = {}
    for cid in cg.keys():
        cdata = expr[cg[cid]]
        new = cdata.mean(axis=1)
        order[cid] = new.idxmax()
    # sort clusters from left to right
    corder = [int(k) for k, v in sorted(order.items(), key=lambda item: item[1])]
    print(f'new order is {corder}')
    if -1 in corder:
        non_cluster = corder.index(-1)
        print(f'-1 cluster is now labeled as {non_cluster}')
    else:
        non_cluster = -1
    return corder, non_cluster


def sort_labels_(cg, expr, new_labels:np.ndarray, save)->np.ndarray:
    """Sort clusters by the position of their peak value, from left to right"""
    corder, non_cluster = sort_cluster(cg, expr)
    #non_cluster = corder.index(-1)
    #print(f'-1 cluster is now labeled as {non_cluster}')
    # update labels
    sorted_labels = new_labels.copy()
    for i, l in enumerate(corder):
        sorted_labels[sorted_labels==l] = i+999
    sorted_labels = sorted_labels - 999 
    np.savetxt(os.path.join(save,"labels.txt"), sorted_labels, fmt='%i', delimiter=',')
    print(f'max cluster number is {sorted_labels.max()}')
    print(f'min cluster number is {sorted_labels.min()}')
    return sorted_labels


def sort_labels(corder:list, new_labels:np.ndarray, save)->np.ndarray:
    """Sort clusters by the position of their peak value, from left to right"""
    sorted_labels = new_labels.copy()
    for i, l in enumerate(corder):
        sorted_labels[sorted_labels==l] = i+999
    sorted_labels = sorted_labels - 999 
    np.savetxt(os.path.join(save,"labels.txt"), sorted_labels, fmt='%i', delimiter=',')
    print(f'max cluster number is {sorted_labels.max()}')
    print(f'min cluster number is {sorted_labels.min()}')
    return sorted_labels


#--------------------------------
#         PCG analysis
#---------------------------------
def load_pcgs(fn: str) ->list:
    with open(fn) as f:
        pcgs = f.read().splitlines()
    return pcgs


def get_pcg_in_cluster(cluster_id: int, data: pd.DataFrame, labels: np.ndarray, pcgs: list) -> list:
    interg = []
    for g in pcgs:
        if g in data.iloc[np.where(labels==cluster_id)[0].tolist()].index.tolist():
            interg.append(g)
    return interg


def get_pcg_in_smes(cluster_id: int, data: pd.DataFrame, labels: np.ndarray, pcgs: list) -> list:
    interg = []
    for g in pcgs:
        gs = [i for i in data.iloc[np.where(labels==cluster_id)[0].tolist()].index.tolist() if i.startswith(g)]
        if len(gs) > 0:
            interg.append(gs)
    interg = list(chain.from_iterable(interg))
    return interg


def load_id_symbol(fn: str) -> dict:
    id_symbol = {}
    with open(fn, "r") as f:
        for line in f.readlines():
            line = line.strip().split("\t")
            if line[0] not in id_symbol:
                id_symbol[line[0]] = line[1]
            #else:
                #print(line)
    return id_symbol


def get_symbol(gid, id_symbol):
    try:
        return id_symbol[gid]
    except KeyError:
        return gid


def get_pcg_clusters(data:pd.DataFrame, labels:np.ndarray, pcgs:list) -> list:
    clusters = []
    for cid in list(set(sorted(labels.tolist()))):
        if len(get_pcg_in_cluster(cid, data, labels, pcgs)) > 0:
            clusters.append(cid)
    return sorted(clusters)


def get_pcg_clusters_smes(data:pd.DataFrame, labels:np.ndarray, pcgs:list) -> list:
    clusters = []
    for cid in list(set(sorted(labels.tolist()))):
        if len(get_pcg_in_smes(cid, data, labels, pcgs)) > 0:
            clusters.append(cid)
    return sorted(clusters)


def get_cluster_data(gene_mat_df:pd.DataFrame, cid:int)->pd.DataFrame:
    return gene_mat_df.iloc[np.where(labels==cid)]


def draw_cluster_heatmap(gene_mat_df:pd.DataFrame, cid:int, labels:np.ndarray, pcgs:list, path):
    #d = get_cluster_data(gene_mat_df, cid)
    genes = get_pcg_in_cluster(cid, gene_mat_df, labels, pcgs)
    vmin = np.min(gene_mat_df)
    vmax = np.max(gene_mat_df)
    if len(genes) >1:
        fig, axs = plt.subplots(nrows=len(genes), ncols=1,figsize=(6, 12))
        for i in range(len(genes)):
            sns.heatmap(gene_mat_df.loc[[genes[i]]], ax=axs[i], vmin=vmin, vmax=vmax, cmap="twilight_shifted")
    else:
        sns.heatmap(gene_mat_df.loc[genes], vmin=vmin, vmax=vmax, cmap="twilight_shifted")
    plt.savefig(os.path.join(path,f"cluster{cid}_pcg.png"))
    plt.close()


def sort_data(order:list, anno:pd.DataFrame, feature:str):
    """Subset and sort annotation data"""
    anno = anno[anno[feature].isin(order)]
    anno[feature] = anno[feature].astype('category')
    anno[feature].cat.set_categories(order, inplace=True)
    anno.sort_values(feature, inplace=True)
    return anno


def plot_time_series(df, symbol, genes, save, order):
    """Plot change of cluster as the time change for each gene"""
    for g in genes:
        fn = ''
        data = df.loc[df['gene']==g]
        data = sort_data(order, data, 'sample')
        clusters = data['cid']
        if g in list(symbol.keys()):
            g = symbol[g]
            if '/' in g:
                g = g.replace('/','_')
        if len(set(clusters)) == 0:
            continue
        elif len(set(clusters)) == 1:
            fn = os.path.join(save, f'{g}_no_cluster_change.png')
        elif len(set(clusters)) > 1:
            fn = os.path.join(save, f'{g}_cluster_changed.png')
        print(f'save to image file {fn}')
        plt.margins(x=0, y=0)
        plt.plot(data['sample'], clusters, '-o')
        plt.title(g)
        plt.xticks(data['sample'], rotation=90)
        plt.yticks(clusters)
        plt.xlabel('time')
        plt.ylabel('cluster id')
        plt.tight_layout()
        plt.savefig(fn)
        plt.close()


def save_gene_cluster(gene_label, fn, save):
    with open(os.path.join(save,fn), 'w') as f:
        for key, value in gene_label.items():
            f.write('%s\t%s\n' % (key, value))


def save_cluster_gene(cluster_genes, fn, save):
    with open(os.path.join(save,fn), 'w') as f:
        for key, value in cluster_genes.items():
            f.write('%s,%s\n' % (key, ",".join(map(str, value))))


# 2023-09-26
def update_recall_labels(labels, recall, total_genes):
    '''update cluster labels for all the previous -1 cluster genes'''
    for gene in recall['gene']:
        gindex = total_genes.index(gene)
        labels[gindex] = recall[recall.gene==gene]['cluster']
    return labels


def update_sorted_labels(corder:list, new_labels:np.ndarray)->np.ndarray:
    """update cluster labels by peak value order"""
    sorted_labels = new_labels.copy()
    for i, l in enumerate(corder):
        sorted_labels[sorted_labels==l] = i+999
    sorted_labels = sorted_labels - 999
    print(f'max cluster number is {sorted_labels.max()}')
    print(f'min cluster number is {sorted_labels.min()}')
    return sorted_labels


def sort_bin(df: pd.DataFrame):
    """
    sort data by binID
    :param df:
    :return:
    """
    og_index = list(df.index)
    new_index = sorted(og_index, key=lambda x: int(float(x.strip('bin'))))
    df = df.reindex(new_index)
    return df


def draw_umap(gene_mat, sorted_labels, non_cluster, save, name, order: list):
    sub_order = order.copy()
    sub_order.remove(non_cluster)
    X = gene_mat[sorted_labels != non_cluster]
    mapper = umap.UMAP().fit(X)
    umap.plot.points(mapper, labels=sorted_labels[sorted_labels != non_cluster])
    plt.legend(sub_order, ncol=2, bbox_to_anchor=(1, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(save, f"{name}_umap.pdf"), format='pdf')
    plt.close()
# 2023-09-26 end


def main():
    example_text = """example:
      --- flat clustering
      python pcg_pattern.py -d WT.csv -i WT
      python pcg_pattern.py -d WT.csv -i WT -fc 15
      --- common db clustering
      python pcg_pattern.py -d WT.csv -i WT -f False -w 0 
      python pcg_pattern.py -d WT.csv -i WT -f False -min 6 -c 0.8 -cor spearman
    """

    parser = argparse.ArgumentParser(description='Detecting PCG pattern along body axis', epilog=example_text, formatter_class=argparse.RawDescriptionHelpFormatter)
    # get data
    parser.add_argument('-d', '--data', metavar='expressionData', type=str, help='path to expression data. columns are genes, rows are features')
    parser.add_argument('-i', '--sample', metavar='sampleID', type=str, help='sample id. e.g. WT')
    # data pre-process
    parser.add_argument('-w', '--window', metavar='windowSize', type=int, default=5, help='[int] window size to filter genes. default is 5. no filtering when window size is set to 0. window size should not be bigger than bin number')
    parser.add_argument('-s', '--scale', nargs="*", default=[3,-3], help='[max, min]. default is [3,-3]')
    parser.add_argument('--pca', type=str, default='False', choices=['False', 'True'], help='if use pca to reduct dimentions of data')
    parser.add_argument('--smooth', type=str, default='True', choices=['False', 'True'], help='if gaussian smooth data')
    parser.add_argument('--sigma', type=int, default=3, help='sigma of the gaussian filter')
    # clustering
    parser.add_argument('-f', '--flat', metavar='flatClustering', type=str, default='False', choices=['False', 'True'], help='[bool] flat hdbscan clustering. default is True')
    parser.add_argument('-fc', '--cluster', metavar='flatClusterNumber', type=int, default=30, help='[int] flat hdbscan clusters number. default is 30')
    parser.add_argument('-min', '--min_cluster', metavar='minClusterSize', default=0, type=int, help='min cluster size of HDBSCAN. default is 0. when all four param are set to 0, HDBSCAN will run on default setting.')
    parser.add_argument('-ms', '--min_sample', metavar='minSampleSize', default=0, type=int, help='minimum sample size of HDBSCAN. default is 0. when all four param are set to 0, HDBSCAN will run on default setting.')
    parser.add_argument('-e', '--epsilon', metavar='epsilon', default=0, type=float, help='epsilon of HDBSCAN. default is 0. when all four param are set to 0, HDBSCAN will run on default setting.')
    # check pcg stats
    parser.add_argument('-p', '--pcg', metavar='pcgList', type=str, default='/dellfsqd2/ST_OCEAN/USER/liyao1/scripts/pcg_ap_pattern/input_data/gene_class.txt', help='path to pcg list')
    parser.add_argument('--total', metavar='totalPCGList', type=str, default='/dellfsqd2/ST_OCEAN/USER/liyao1/scripts/pcg_ap_pattern/input_data/unique_all_pcg_smes.txt', help='path to reddien pcg list')
    parser.add_argument('-y', '--symbol', metavar='geneSymbol', type=str, default='/dellfsqd2/ST_OCEAN/USER/liyao1/scripts/pcg_ap_pattern/input_data/pcg_symbol_smes.txt', help='path to symbol-gene list')
    parser.add_argument('-o', '--order', metavar='clusterOrderList', nargs="*", type=int, help='order of clusters to show in heatmap')
    # re-grouping genes
    parser.add_argument('-c', '--cutoff', metavar='cutoff', type=float, default=0.9, help='[float] cutoff of a strong pearson correlation coefficient value')
    parser.add_argument('-cor', '--correaltion', metavar='correaltionMethod', type=str, default='pearson', choices=['pearson', 'spearman'], help='method to calculate correlation coefficient b/w genes. options are pearson and spearman')
    # parse arguments 
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    
     
    print('--------------------------------------------------------')
    print(f'{args.sample} data: \nif use flat clustering is {args.flat}\nif smooth data is {args.smooth}\nif use PCA is {args.pca}\nCorrelation Coefficient method is {args.correaltion}')
    print('--------------------------------------------------------')
    #-------------------#
    # set up work space #
    #-------------------#
    print('set up work space')
    save = ''
    if args.flat == 'True':
        save = os.path.join(os.getcwd(), "output", f'{args.sample}_flat_{args.cluster}')
    elif args.flat == 'False':
        if args.min_cluster == 0 and args.epsilon == 0 and args.min_sample == 0: #and args.max_cluster == 0:
            save = os.path.join(os.getcwd(), "output", f'{args.sample}_default')
        else:
            save = os.path.join(os.getcwd(), "output", f'{args.sample}_min{args.min_cluster}_epsilon{args.epsilon}_minsample{args.min_sample}')
    if not os.path.exists(save):
        os.makedirs(save)
        print("create output directory")
    else:
        print("output directory exists")
    print('--------------------------------------------------------')


    #-----------#
    # load data #
    #-----------#
    print('Loading data...')
    #data = pd.read_csv(args.data)
    data = read_expr(args.data)
    if data.isnull().sum().sum() > 0:
        raise ValueError('warning: data contains NA value')
        #print('warning: data contains NA value')
        #return False
    
    
    #-----------------#
    #     PCGs        #
    #-----------------#
    # load PCGs data
    #pcgs = load_pcgs(args.pcg)
    #pcgs = read_list(args.pcg)
    pcg_class = pd.read_csv(args.pcg, delimiter='\t', header=None)
    pcg_class.columns = ['gene', 'class', 'symbol']
    pcgs = list(pcg_class['gene'])
    #id_symbol = load_id_symbol(args.symbol)
    total = read_list(args.total)


    #-----------------------#
    # filter genes and data #
    #-----------------------#
    print('--------------------------------------------------------')
    if args.window <= 0:  
        print('Filtering data: remove genes contain all zero values')
        # do not perform filtering
        # however, still need to remove all zero genes
        gene_remains = []
        for g in data.columns:
            if (data[g]==0).all():
                pass
            else:
                gene_remains.append(g)
    else:
        print(f"Filtering data: window size is {args.window}")
        gene_remains = filter_genes(data, args.window)
    print(f'{len(gene_remains)} genes remains')
    sub_data = data[gene_remains]
    genes = gene_remains
    if args.smooth == 'True':
        print('Smoothing data...')
        sub_data = smooth(sub_data, args.sigma)


    #-------------------------#
    # scaling expression data #
    #-------------------------#
    print("Normalizing data...")
    new_norm = norm_scale(sub_data, args.scale[0], args.scale[1])
    # smooth data
    #new_norm = smooth(new_norm, 3)
    new_norm.to_csv(os.path.join(save,"scaled_cell_norm.csv"), header=True, index=False)


    #------------#
    # clustering #
    #------------#
    gene_mat_df = new_norm.T
    gene_mat = gene_mat_df.to_numpy()
    if args.pca == 'False':
        print('no dimention reduction')
        cluster_gene_mat = gene_mat.copy()
    else:
        print('Dimention reduction via PCA...')
        # gene_mat: columns as bins, rows are genes
        cluster_gene_mat = pca_2d(new_norm)
    print('--------------------------------------------------------')
    print(f"Clustering using HDBSCAN...")
    if args.flat == 'True':
         clusterer = flat.HDBSCAN_flat(cluster_gene_mat, args.cluster)   #, preiction_data=True)
         #flat.approximate_predict_flat(clusterer, points_to_predict, 30)
    elif args.flat == 'False':
        clusterer = doHDBSCAN(cluster_gene_mat, int(args.min_cluster), args.epsilon, args.min_sample)
    print('clustering done')
    # get clusters properties
    labels = clusterer.labels_
    clusters = list(set(sorted(labels.tolist())))
    cluster_num = clusterer.labels_.max() + 1
    print(f'max cluster number is {clusterer.labels_.max()}')
    print(f'min cluster number is {clusterer.labels_.min()}')
    print(f'{len(labels[labels==-1])}/{len(labels)} genes were not assigned to any cluster')
    print('--------------------------------------------------------')
    gene_label = cal_gene_cluster(gene_remains, labels)
    cluster_genes = cal_cluster_genes(gene_mat_df, labels)
    draw_heatmap(gene_mat, labels, f'1.{args.sample}', save)
    np.savetxt(os.path.join(save,"og_labels.txt"), labels, fmt='%i', delimiter=',')
    with open(os.path.join(save,f"og_gene_cluster.txt"), 'w') as f:
        for key, value in gene_label.items():
            f.write('%s\t%s\n' % (key, value))
    with open(os.path.join(save,f"og_cluster_genes.csv"), 'w') as f:
        for key, value in cluster_genes.items():
            f.write('%s,%s\n' % (key, ",".join(map(str, value)))) 
    with open(os.path.join(save,'no_cluster_genes.txt'), 'w') as f:
        f.writelines('\n'.join(cluster_genes[-1]))
    meta = pd.DataFrame({'cluster':list(gene_label.values()), 'gene':list(gene_label.keys())})
    gene_mat_2d = doUMAP(gene_mat, labels)
    if 'gene_mat_2d' in locals():
        meta['umapx'] = gene_mat_2d[:,0]
        meta['umapy'] = gene_mat_2d[:,1]
    #meta['time'] = meta['gene'].str.split('_',expand=True)[1]
    #meta['gene'] = meta['gene'].str.split('_',expand=True)[0]
    meta.to_csv(os.path.join(save,'og_meta_data.csv'),index=False) 


    #----------------------# 
    # draw cluster in umap #
    #----------------------#
    gene_mat_2d = doUMAP(gene_mat, labels)
    draw_cluster(gene_mat_2d, clusterer, args.sample, save)
    X = gene_mat[labels!=-1]
    mapper = umap.UMAP().fit(X)
    umap.plot.points(mapper, labels=labels[labels!=-1])
    plt.savefig(os.path.join(save,f"{args.sample}_umap.png"))
    plt.close()


    #----------------------#
    # check ref pcg status #
    #----------------------#
    pcg_clusters = get_pcg_clusters(gene_mat_df, labels, pcgs)
    print(f'24 PCG belong to clusters {pcg_clusters}')
    for p in pcgs:
        print(p, get_cluster(p, gene_label))
    print('--------------------------------------------------------')
    
    
    #--------------------------------------------------------------#
    # re-clustering genes based on pearson correlation coefficient #
    #--------------------------------------------------------------#
    start_ = time.time()
    print('Calculating correlation coefficiency b/w unclustered genes and clustered genes')
    if len(labels[labels==-1]) > 1:
        new_labels = re_clustering(gene_mat_df, labels, cluster_genes, args.correaltion, args.cutoff)
        print(f'reassigned {len(np.where(labels == -1)[0]) - len(np.where(new_labels == -1)[0])} genes to a cluster')
        end_ = time.time()
        print(f"Time used to execute: {end_ - start_} seconds")
    else:
        new_labels = labels.copy()


    #---------------#
    # Visualization #
    #---------------#
    # Draw heatmap of all clusters
    draw_heatmap(gene_mat, new_labels, f'2.{args.sample}_recall_noorder', save)
    # sort clusters by the location of its max value
    print('--------------------------------------------------------')
    print('Sorting clusters')
    corder, non_cluster = sort_cluster(cluster_genes, new_norm)
    # update labels
    sorted_labels = sort_labels(corder, new_labels, save)
    
    
    #-------------------------------------------# 
    # draw cluster in umap by new sorted labels #
    #-------------------------------------------#
    gene_mat_2d = doUMAP(gene_mat, sorted_labels)
    draw_cluster(gene_mat_2d, clusterer, f'{args.sample}_sorted', save)
    X = gene_mat[sorted_labels!=non_cluster]
    mapper = umap.UMAP().fit(X)
    umap.plot.points(mapper, labels=sorted_labels[sorted_labels!=non_cluster])
    plt.savefig(os.path.join(save,f"{args.sample}_sorted_umap.png"))
    plt.close()


    #---------------------------------------------# 
    # update data based on new clustering result  #
    #---------------------------------------------#
    gene_label = cal_gene_cluster(gene_remains, sorted_labels)
    cluster_genes = cal_cluster_genes(gene_mat_df, sorted_labels)
    with open(os.path.join(save,f"gene_cluster.txt"), 'w') as f:
        for key, value in gene_label.items():
            f.write('%s\t%s\n' % (key, value))
    with open(os.path.join(save,f"cluster_genes.csv"), 'w') as f:
        for key, value in cluster_genes.items():
            f.write('%s,%s\n' % (key, ",".join(map(str, value)))) 


    #---------------#
    # Visualization #
    #---------------#
    #sorted_labels = np.loadtxt(os.path.join(save,"labels.txt"))
    new_order = list(range(int(sorted_labels.max())+1))
    draw_heatmap(gene_mat, sorted_labels, f'3.{args.sample}', save, new_order)


    
    #-------------get Meta data----------------    
    meta = pd.DataFrame({'cluster':list(gene_label.values()), 'gene':list(gene_label.keys())})
    if 'gene_mat_2d' in locals():
        meta['umapx'] = gene_mat_2d[:,0]
        meta['umapy'] = gene_mat_2d[:,1]
    #meta['time'] = meta['gene'].str.split('_',expand=True)[1]
    #meta['gene'] = meta['gene'].str.split('_',expand=True)[0]
    meta.to_csv(os.path.join(save,'meta_data.csv'),index=False) 
    #-----------------#
    #     PCGs        #
    #-----------------#
    print('--------------------------------------------------------')
    print('new cluster of pcgs')
    for p in pcgs:
        print(p, get_cluster(p, gene_label))
    pcg_meta = meta[meta['gene'].isin(pcgs)]
    pcg_meta = pcg_meta.merge(pcg_class, on='gene')
    pcg_meta.to_csv(os.path.join(save,'24_pcg_meta_data.csv'),index=False)
    print(f'{len(set(pcg_meta[pcg_meta["cluster"]==non_cluster]["gene"]))} pcgs are in -1 cluster')
    total_meta = meta[meta['gene'].isin(total)]
    total_meta.to_csv(os.path.join(save,'reddien_meta_data.csv'),index=False)
    print(f'{len(set(total_meta[total_meta["cluster"]==non_cluster]["gene"]))} reddien pcgs are in -1 cluster')
    print('--------------------------------------------------------')


if __name__ == '__main__':
    start = time.time()
    main()
    end = time.time()
    print(f"Time used to execute: {end - start}")


