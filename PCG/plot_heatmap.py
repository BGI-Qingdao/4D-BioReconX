import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import umap
import umap.plot
import seaborn as sns
import os

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
plt.rcParams["axes.grid"] = False


def plot_heatmap(gene_mat, labels:np.ndarray, name:str, orders:list, fntsz=6):
    c = 'RdYlBu_r'
    vmin = np.min(gene_mat)
    vmax = np.max(gene_mat)
    cluster_num = labels.max()+1
    # plotting
    cm = 1/2.54
    fig, axs = plt.subplots(nrows=len(orders), ncols=1,figsize=(4.1*cm, 5.8*cm))
    for i in range(len(orders)):
        sns.heatmap(gene_mat[np.where(labels==orders[i])], ax=axs[i], vmin=vmin, vmax=vmax, cmap=c, cbar=False, rasterized=True, edgecolor="none", linewidths=0.0)
        axs[i].set_xticks([])
        gene_num = len(labels[labels==orders[i]])
        #if i % 3 == 0:
        axs[i].set_yticks([gene_num//2])
        axs[i].set_yticklabels([gene_num],rotation=0,fontsize=fntsz)#fontsize = 'small')
        #else:
        #    axs[i].set_yticks([])
    #sns.heatmap(gene_mat[np.where(labels==orders[-1])], ax=axs[-1], vmin=vmin, vmax=vmax, cmap=c, cbar=False, rasterized=True)
    
    # show every out of 3
    #if (len(axs)-1) % 3 == 0:
    axs[-1].set_yticks([len(labels[labels==orders[-1]])//2])
    axs[-1].set_yticklabels([len(labels[labels==orders[-1]])],rotation=0,fontsize=fntsz)#fontsize = 'small')
    #else:
    #    axs[-1].set_yticks([])
    axs[-1].set_xticklabels(axs[-1].get_xticklabels(), rotation=90, fontsize=6) 
    # show every y tick labels
    #plt.yticks(np.arange(max(labels)+1,2))
    plt.savefig(f'{name}_final.pdf', bbox_inches='tight', format='pdf')


def draw_heatmap_gc(gc, expr, orders: list=None):
    vmin = np.min(gene_mat)
    vmax = np.max(gene_mat)
    clusters = set(gc[1])
    #clusters = sorted(list(clusters))
    # plotting
    fig, axs = plt.subplots(nrows=len(clusters), ncols=1,figsize=(6, 12))
    
    for i in range(len(clusters)-1):
        sns.heatmap(expr[list(gc[gc[1]==clusters[i]][0])].T.to_numpy(), ax=axs[i+1], vmin=vmin, vmax=vmax, cmap='twilight_shifted', cbar=False)
        axs[i+1].set_xticks([])
        axs[i+1].set_yticks([])
    sns.heatmap(expr[list(gc[gc[1]==clusters[-1]][0])].T.to_numpy(), ax=axs[-1], vmin=vmin, vmax=vmax, cmap='twilight_shifted', cbar=False)
    axs[-1].set_yticks([])
    plt.close()


'''
python plot_heatmap.py scaled_cell_norm.csv labels_sorted.txt ap 6
'''


if __name__ == '__main__':
    colors=["#00193a","#002b53","#023f73","#034780","#7a0213","#a10220","#bf0a26","#cd0c2b","#131313","#262626"]
    #twilight_shifted'
    c = 'RdYlBu_r'
    
    
    expr = pd.read_csv(sys.argv[1])
    gene_mat = expr.T.to_numpy()
    labels = np.loadtxt(sys.argv[2])
    mode = sys.argv[3]
    font_size = sys.argv[4]
    if mode == 'ap':
        ap_order = [1,0,3,4,5,6,7,8,9,10,11,12,13,14,16,15]
        #ap_order = [1,0,3,4,5,6,7,8,9,10,11]
        #ap_order = [1,0,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
        plot_heatmap(gene_mat, labels, 'ap_color', ap_order, font_size)
    elif mode == 'ml':
        #ml_order = [22, 21, 24, 20, 25, 23, 17, 30, 28, 29, 16, 19, 26, 27, 0, 1, 18, 5, 3, 8, 4, 7, 9, 11, 10, 12, 6, 14, 13, 15, 2]
        #ml_order = [5,4,25,27,28,24,15,1,3,23,10,8,9,2,17,30,26,11,12,14,13,18,19,20,29,21,6,16,0,22,7]
        ml_order = list(range(int(labels.max())+1))
        plot_heatmap(gene_mat, labels, 'ml', ml_order, font_size)
    elif mode == 'dv':
        order = list(range(int(labels.max())+1))
        plot_heatmap(gene_mat, labels, 'dv', order, font_size)

