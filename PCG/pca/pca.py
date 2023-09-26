import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from sklearn.decomposition import PCA
from util import read_expr, read_list, read_symbols, get_symbol
import argparse
import os
import matplotlib as mpl
mpl.rcParams['axes.spines.left'] = True
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.bottom'] = True
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'


def rank_genes(gene_mat:np.ndarray):
    """
    Rank candidate genes by egienvalues of their covariance
    """
    print('calculating covariance matrix...')
    cov = np.cov(gene_mat)
    print('calculating eigenvalues and eigenvectors...')
    evals, evecs = np.linalg.eig(cov)
    idx = evals.argsort()[::-1]
    evals = evals[idx]
    evecs = evecs[idx]
    return evals, evecs, idx


def draw_pcg_scatter(project_data, name=''):
    fig = plt.figure()
    plt.scatter(project_data[:, 0], project_data[:, 1], marker='o')
    plt.savefig(f'{name}pca.png')
    plt.close()


def draw_pc_plot(pdata, genes, pc1:int=1, pc2:int=2, sample='WT'):
    for i, label in enumerate(genes):
        plt.annotate(label, (pdata[i,pc1-1], pdata[i,pc2-1]))
    plt.scatter(pdata[:,pc1-1], pdata[:,pc2-1])
    plt.xlabel(f'PC{pc1}')
    plt.ylabel(f'PC{pc2}')
    plt.tight_layout()
    plt.savefig(f'{sample}_pc{pc1}_pc{pc2}.pdf', format='pdf')
    plt.close()


def draw_pc_plot_symbol(pdata, genes, symbols, pc1:int=1, pc2:int=2, sample='WT'):
    for i, label in enumerate(genes):
        # 2022-07-13
        if label in symbols.keys():
            plt.annotate(symbols[label], (pdata[i,pc1-1], pdata[i,pc2-1]))
    plt.scatter(pdata[:,pc1-1], pdata[:,pc2-1])
    plt.xlabel(f'PC{pc1}')
    plt.ylabel(f'PC{pc2}')
    plt.tight_layout()
    plt.savefig(f'{sample}_pc{pc1}_pc{pc2}.pdf', format='pdf')
    plt.close()


def plot_pcg_class(pcg_class, expr, pca, pc1_top=3, pc2_top=3):
    classes = list(set(pcg_class[1]))
    for c in classes:
        pcgs = list(pcg_class[pcg_class[1]==c][0])
        # intersection
        pcgs = [i for i in pcgs if i in expr.columns]
        try:
            pcg_expr = expr[pcgs].to_numpy()
        except KeyError:
            #print(pcgs)
            continue
        pp = pca.transform(pcg_expr.T)
        plt.scatter(pp[:,0], pp[:,1])
    plt.legend(classes, bbox_to_anchor=(1, 1), loc='upper left')
    plt.xlabel('PC1 (64.7%)')
    plt.ylabel('PC2 (24.7%)')
    plt.xlim([-15,15])
    plt.tight_layout()
    plt.savefig('pcg_classes.pdf', format='pdf')
    plt.close()


def plot_pcg_class_symbol(pcg_class, expr, pca, symbols):
    classes = list(set(pcg_class[1]))
    for c in classes:
        pcgs = list(pcg_class[pcg_class[1]==c][0])
        pcgs = [i for i in pcgs if i in expr.columns]
        try:
            pcg_expr = expr[pcgs].to_numpy()
        except KeyError:
            print(pcgs)
            continue
        pp = pca.transform(pcg_expr.T)
        plt.scatter(pp[:,0], pp[:,1])
        for i, label in enumerate(pcgs):
            # 2022-07-13
            if label in symbols.keys():
                plt.annotate(symbols[label], (pp[i,0], pp[i,1]))
    plt.legend(classes, bbox_to_anchor=(1, 1), loc='upper left')
    plt.xlabel('PC1 (64.7%)')
    plt.ylabel('PC2 (24.7%)')
    plt.xlim([-15,15])
    plt.tight_layout()
    plt.savefig('pcg_classes_symbol.pdf', format='pdf')
    plt.close()
   


def plot_class_timeseries(bkground, pcgs_5c, order, genes, expr, pca):
    """
    pcgs_5c
    order:
    genes: (list) total genes (of all samples)
    expr: total data (of all samples)
    """
    # get data of each class
    classes = list(set(pcgs_5c[1]))
    for c in classes:
        # plot background data
        bkcolors = plt.cm.rainbow(np.linspace(0, 1, 4))
        bk1 = bkground[np.where((bkground[:,0]<0)&(bkground[:,1]<0))]
        plt.scatter(bk1[:,0], bk1[:,1], color=bkcolors[0])
        bk2 = bkground[np.where((bkground[:,0]>0)&(bkground[:,1]<0))]
        plt.scatter(bk2[:,0], bk2[:,1], color=bkcolors[1])
        bk3 = bkground[np.where((bkground[:,0]<0)&(bkground[:,1]>0))]
        plt.scatter(bk3[:,0], bk3[:,1], color=bkcolors[2])
        bk4 = bkground[np.where((bkground[:,0]>0)&(bkground[:,1]>0))]
        plt.scatter(bk4[:,0], bk4[:,1], color=bkcolors[3])
        #plt.legend(['head','tail','pharynx-head','pharynx-tail'], bbox_to_anchor=(1.05, 1.05))
        cpcgs = list(pcgs_5c[pcgs_5c[1]==c][0])
        for pcg in cpcgs:
            sub_gene = [i for i in genes if str(pcg) in i]
            # sort in time order
            sub_gene = [p for o in order for p in sub_gene if o in p]
            if len(sub_gene) != 9:
                continue
            # get drawing data
            try:
                p_sub_data_mat = pca.transform(expr[sub_gene].T.to_numpy())
            except ValueError:
                print('plot pcg class series', pcg)
                continue
            # prepare quiver data
            x = p_sub_data_mat[:-1,0] # no need to input last point cuz the function will calculate its value by dx, dy
            y = p_sub_data_mat[:-1,1]
            Fx = np.diff(p_sub_data_mat[:,0])
            Fy = np.diff(p_sub_data_mat[:,1])
            colors = plt.cm.rainbow(np.linspace(0, 1, 9))
            plt.quiver(x, y, Fx, Fy, angles='xy', scale_units='xy', color=colors, scale=1, width=0.005) 

        # total plot attr
        plt.title(c)
        plt.tight_layout()
        plt.savefig(f'{c}_timeseries.pdf', format='pdf')
        plt.close()
    print('plot_class_timeseries DONE')


def plot_time_series(PCGS, order, genes, expr, pca, symbols):
    for pcg in PCGS:
        sub_gene = [i for i in genes if str(pcg) in i]
        # sort in time order
        sub_gene = [p for o in order for p in sub_gene if o in p]
        print(len(sub_gene))
        if len(sub_gene) == 0:
            continue
        
        #if len(sub_gene) != 9:
        #    continue
        # prepare drawing data
        sub_data = expr[sub_gene]
        sub_data_mat = expr[sub_gene].to_numpy()
        p_sub_data_mat = pca.transform(sub_data_mat.T)
        
        x = p_sub_data_mat[:-1,0] # no need to input last point cuz the function will calculate its value by dx, dy
        y = p_sub_data_mat[:-1,1]
        Fx = np.diff(p_sub_data_mat[:,0])
        Fy = np.diff(p_sub_data_mat[:,1])
        colors = plt.cm.rainbow(np.linspace(0, 1, 9))
        plt.scatter(p_sub_data_mat[:,0], p_sub_data_mat[:,1])
        for i, g in enumerate(sub_gene):
            plt.annotate(g.split('_')[1], (p_sub_data_mat[i,0], p_sub_data_mat[i,1]))
        
        
        
        plt.quiver(x, y, Fx, Fy, angles='xy', scale_units='xy', color=colors, scale=1, width=0.005) 
        
        
        #plt.xlim([-20,20])
        #plt.ylim([-20,20])
        plt.title(get_symbol(pcg, symbols))
        plt.xlabel('PC1 (64.7%)')
        plt.ylabel('PC2 (24.7%)')
        plt.tight_layout()
        plt.savefig(f'{get_symbol(pcg, symbols)}_timeseries.pdf', format='pdf')
        plt.close()
    print('plot_time_series DONE')


def plot_time_series2(PCGS, order, genes, expr, pca, symbols):
    for pcg in PCGS:
        sub_gene = [i for i in genes if str(pcg) in i]
        # sort in time order
        sub_gene = [p for o in order for p in sub_gene if o in p]
        if len(sub_gene) ==0 :
            continue
        # prepare drawing data
        sub_data = expr[sub_gene]
        sub_data_mat = expr[sub_gene].to_numpy()
        p_sub_data_mat = pca.transform(sub_data_mat.T)
        #x = p_sub_data_mat[:-1,0] # no need to input last point cuz the function will calculate its value by dx, dy
        #y = p_sub_data_mat[:-1,1]
        #Fx = np.diff(p_sub_data_mat[:,0])
        #Fy = np.diff(p_sub_data_mat[:,1])
        #colors = plt.cm.rainbow(np.linspace(0, 1, 9))
        plt.scatter(p_sub_data_mat[:,0], p_sub_data_mat[:,1], c='black', s=100)
        # smooth line


        # arrow
        #plt.quiver(x, y, Fx, Fy, angles='xy', scale_units='xy', color=colors, scale=1, width=0.005) 


        plt.title(get_symbol(pcg, symbols))
        plt.xlabel('PC1 (64.7%)')
        plt.ylabel('PC2 (24.7%)')
        plt.tight_layout()
        plt.savefig(f'{get_symbol(pcg, symbols)}_timeseries2.pdf', format='pdf')
        plt.close()


def cal(pcg, order, genes, expr, pca, symbols):
    sub_gene = [i for i in genes if str(pcg) in i]
    # sort in time order
    sub_gene = [p for o in order for p in sub_gene if o in p]
    p_sub_data_mat = pca.transform(expr[sub_gene].T.to_numpy())
    
    dists = {}
    for i in range(1, p_sub_data_mat.shape[0]-1):
        dist = np.linalg.norm(p_sub_data_mat[i+1,0:2]-p_sub_data_mat[i,0:2])
        dists[dist] = order[i+1]
    print(pcg, get_symbol(pcg, symbols), max(dists.keys()), dists[max(dists.keys())])


def plot_pca(vector, data, pc1=1, pc2=2, sample='90'):
    pdata = np.dot(data.to_numpy().T, vector.T)
    plt.scatter(pdata[:, pc1-1], pdata[:, pc2-1], marker='o')
    for i, label in enumerate(genes):
        plt.annotate(label, (pdata[i,pc1-1], pdata[i,pc2-1]))
    plt.xlabel(f'PC{pc1}')
    plt.ylabel(f'PC{pc2}')
    plt.tight_layout()
    plt.savefig(f'{sample}_pc{pc1}_pc{pc2}.pdf', format='pdf')
    plt.close()


def plot_pca_symbol(vector, data, anno, pc1=1, pc2=2, sample='90'):
    pdata = np.dot(data.to_numpy().T, vector[0:3, :].T) # correct way
    plt.scatter(pdata[:, pc1-1], pdata[:, pc2-1], marker='o')
    for i, label in enumerate(genes):
        plt.annotate(anno[label], (pdata[i,pc1-1], pdata[i,pc2-1]))
    plt.xlabel(f'PC{pc1}')
    plt.ylabel(f'PC{pc2}')
    plt.tight_layout()
    plt.savefig(f'{sample}_pc{pc1}_pc{pc2}.pdf', format='pdf')
    plt.close()


def main():
    example_text = '''example:
      python pca.py -i total.csv -t WT.csv -p pcg.txt -c pcg_class.txt -n test_WT
      python pca.py -i total.csv -t WT.csv -g filtered_gene.txt -p pcg.txt -c pcg_class.txt -n test_WT -s symbols.txt'''
    parser = argparse.ArgumentParser(description = 'PCA', epilog=example_text, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--name", '-n', type=str, help='sample name')
    parser.add_argument("--symbol", '-s', type=str, help='')
    # train pca model
    parser.add_argument("--trainingData", '-t', type=str, help='training dataset, cols are genes, rows are bins')
    # plot only pcg of WT sample
    parser.add_argument("--pcg", '-p', type=str, help='known pcg list file')
    parser.add_argument("--genes", '-g', type=str, help='target gene list file')
    parser.add_argument("--classes", '-c', type=str, help='pcg class file')
    # time series, total data
    parser.add_argument("--data", '-i', type=str, help='data expression file to project in PC space, cols are genes, rows are bins')
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    if not os.path.exists(args.name):
        os.makedirs(args.name)

    #----------#
    # Datasets #
    #----------#
    symbols = read_symbols(args.symbol, '\t')
    print(f'training data is {args.trainingData}')
    #print(f'predicting data is : {args.data}')
    expr = read_expr(args.data)
    gene_mat = expr.to_numpy()
    # interested genes
    genes = read_list(args.genes)
    print(genes)
    # known PCGs
    ref_pcg = read_list(args.pcg)
    total_ref_pcg = [j for i in ref_pcg for j in expr.columns if i in j]

    # WT data
    WT_data = read_expr(args.trainingData)
    ref_pcg = [i for i in ref_pcg if i in WT_data.columns]
    genes = [i for i in genes if i in WT_data.columns]
    print(genes)
    total_genes = ref_pcg + genes

    # real training dataset
    #expr_train = WT_data[genes]
    expr_train = WT_data[ref_pcg]
    gene_mat_train = expr_train.to_numpy()


    #---------------#
    #  sklearn PCA  #
    #---------------#
    print('performing pca')
    pca = PCA()
    p_data = pca.fit_transform(gene_mat_train.T)
    eigen_value = pca.explained_variance_
    eigen_vector = pca.components_
    ratio = pca.explained_variance_ratio_
    print(f'first two are {ratio[0:2]}')
    print(f'first two sum is {sum(ratio[0:2])}')
    np.savetxt(f'{args.name}/eigen_value_WT.txt', eigen_value)
    np.savetxt(f'{args.name}/eigen_vector_WT.txt', eigen_vector)


    #--------------------------#
    # Project data to PC space #
    #    PCG     #
    #------------#
    # new genes
    new_expr = WT_data[genes].to_numpy()
    pdata = pca.transform(new_expr.T)
    #draw_pc_plot(pdata, genes, 1,2, f'{args.name}/candi')
    pclass = pd.read_csv(args.classes, header=None, delimiter='\t')
    pcgs = list(set(pclass[0]).intersection(set(expr_train.columns)))
    print('making new pcgs pc plots...')
    #plot_pcg_class(pclass, WT_data, pca)
    #plot_pcg_class_symbol(pclass, WT_data, pca, symbols)


    #-------------------------------------------------------------#
    # II. Plot data of all samples ......
    #-------------------------------------------------------------#
    print('making time series plots...')
    order =  ["WT", "0hpa1", "12hpa2", "36hpa2", "3dpa2", "5dpa1", "7dpa2", "10dpa1", "14dpa1"]
    total_genes = list(expr.columns)
    #plot_class_timeseries(pdata, pclass, order, total_genes, expr, pca)
    #plot_time_series(genes, order, total_genes, expr, pca, symbols)
    plot_time_series(genes, order, total_genes, expr, pca, symbols)
    plot_time_series2(genes, order, total_genes, expr, pca, symbols)


if __name__ == '__main__':
    main()



