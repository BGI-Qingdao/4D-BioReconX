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
#plt.rcParams['font.size'] = '18'


def plot_pcg_class(pcg_class, expr, pca, symbols,fntsz=10):
    classes = list(set(pcg_class[1]))
    cm = 1/2.54
    #plt.figure(figsize=(6.4*cm,4.1*cm))
    for c in classes:
        pcgs = list(pcg_class[pcg_class[1]==c][0])
        # intersection
        pcgs = [i for i in pcgs if i in expr.columns]
        try:
            pcg_expr = expr[pcgs].to_numpy()
        except KeyError:
            continue
        pp = pca.transform(pcg_expr.T)
        plt.scatter(pp[:,0], pp[:,1])
        for i, label in enumerate(pcgs):
            if label in symbols.keys():
                plt.annotate(symbols[label], (pp[i,0], pp[i,1]), fontsize=fntsz)
    plt.legend(classes, bbox_to_anchor=(1, 1), loc='upper left', fontsize=fntsz)
    plt.xlabel('PC1 (64.7%)',fontsize=fntsz)
    plt.ylabel('PC2 (24.7%)',fontsize=fntsz)
    plt.xlim([-14.5,16])
    plt.xticks(fontsize=fntsz)
    plt.yticks(fontsize=fntsz)
    plt.tight_layout()
    plt.savefig('pcg_classes.pdf', format='pdf')
    plt.close()


def plot_time_series2(PCGS, order, genes, expr, pca, symbols,fntsz=10):
    for pcg in PCGS:
        sub_gene = [i for i in genes if str(pcg) in i]
        # sort in time order
        sub_gene = [p for o in order for p in sub_gene if o in p]
        if len(sub_gene) != 9:
            continue
        # prepare drawing data
        sub_data = expr[sub_gene]
        sub_data_mat = expr[sub_gene].to_numpy()
        p_sub_data_mat = pca.transform(sub_data_mat.T)
        #x = p_sub_data_mat[:-1,0] # no need to input last point cuz the function will calculate its value by dx, dy
        #y = p_sub_data_mat[:-1,1]
        #Fx = np.diff(p_sub_data_mat[:,0])
        #Fy = np.diff(p_sub_data_mat[:,1])
        plt.scatter(p_sub_data_mat[:,0], p_sub_data_mat[:,1], c='black', s=100)
        for i, g in enumerate(sub_gene):
            plt.annotate(g.split('_')[1], (p_sub_data_mat[i,0], p_sub_data_mat[i,1]),fontsize=fntsz)
        # smooth line
        #new_x = gaussian_filter1d(p_sub_data_mat[:,0],sigma=1)
        #new_y = gaussian_filter1d(p_sub_data_mat[:,1],sigma=1)
        #plt.plot(new_x, new_y,cmap='viridis')
        plt.title(get_symbol(pcg, symbols))
        plt.xlabel('PC1 (64.7%)',fontsize=fntsz)
        plt.ylabel('PC2 (24.7%)',fontsize=fntsz)
        plt.xticks(fontsize=fntsz)
        plt.yticks(fontsize=fntsz)
        plt.tight_layout()
        plt.savefig(f'{get_symbol(pcg, symbols)}_timeseries__.pdf', format='pdf')
        plt.close()



def main():
    example_text = '''example:
      python plot_pca.py -i total.csv -t WT.csv -p pcg.txt -c pcg_class.txt -n test_WT
      python plot_pca.py -i total.csv -t WT.csv -g filtered_gene.txt -p pcg.txt -c pcg_class.txt -n test_WT -s symbols.txt'''
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

    #if not os.path.exists(args.name):
    #    os.makedirs(args.name)

    #----------#
    # Datasets #
    #----------#
    symbols = read_symbols(args.symbol, '\t')
    print(f'training data is {args.trainingData}')
    expr = read_expr(args.data)
    gene_mat = expr.to_numpy()
    # interested genes
    genes = read_list(args.genes)
    # known PCGs
    ref_pcg = read_list(args.pcg)
    total_ref_pcg = [j for i in ref_pcg for j in expr.columns if i in j]

     # WT data
    WT_data = read_expr(args.trainingData)
    ref_pcg = [i for i in ref_pcg if i in WT_data.columns]
    genes = [i for i in genes if i in WT_data.columns]
    total_genes = ref_pcg + genes

     # real training dataset
    expr_train = WT_data[genes]
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

     #--------------------------#
     # Project data to PC space #
     #    PCG     #
     #------------#
     # new genes
    new_expr = WT_data[genes].to_numpy()
    pdata = pca.transform(new_expr.T)
    pclass = pd.read_csv(args.classes, header=None, delimiter='\t')
    pcgs = list(set(pclass[0]).intersection(set(expr_train.columns)))
    plot_pcg_class(pclass, WT_data, pca, symbols,14)


    #-------------------------------------------------------------#
    # II. Plot data of all samples ......
    #-------------------------------------------------------------#
    print('making time series plots...')
    order =  ["WT", "0hpa1", "12hpa2", "36hpa2", "3dpa2", "5dpa1", "7dpa2", "10dpa1", "14dpa1"]
    total_genes = list(expr.columns)
    plot_time_series2(ref_pcg, order, total_genes, expr, pca, symbols,14)


if __name__ == '__main__':
    main()

