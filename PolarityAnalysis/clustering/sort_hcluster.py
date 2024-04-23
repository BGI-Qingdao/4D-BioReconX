import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
import sys


def get_newick(node, parent_dist, leaf_names, newick='') -> str:
    """
    Convert sciply.cluster.hierarchy.to_tree()-output to Newick format.

    :param node: output of sciply.cluster.hierarchy.to_tree()
    :param parent_dist: output of sciply.cluster.hierarchy.to_tree().dist
    :param leaf_names: list of leaf names
    :param newick: leave empty, this variable is used in recursion.
    :returns: tree in Newick format
    """
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parent_dist - node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
        newick = get_newick(node.get_right(), node.dist, leaf_names, newick=",%s" % (newick))
        newick = "(%s" % (newick)
        return newick

#tree = hierarchy.to_tree(Z, False)
#get_newick(tree, tree.dist, leaf_names)

'''
python sort_hcluster.py pcgs.11
'''


if __name__ == '__main__':
    data = pd.read_csv('pcg.csv',sep='\t',header=0)
    pcgs_info = pd.read_csv(sys.argv[1], delimiter='\t', header=None)
    pcgs_info.columns = ['symbol', 'id', 'classes']
    time_list2 = ['WT', '0hpa','12hpa','36hpa','3dpa','5dpa','7dpa','10dpa','14dpa']

    time_list = ['WT','0hpa1','12hpa2','36hpa2','3dpa2','5dpa1','7dpa2','10dpa1','14dpa1']
    for c in list(set(pcgs_info.classes)):
        A_data = []
        P_data = []
        index = []
        class_pcgs_info = pcgs_info[pcgs_info.classes==c]
        for name in class_pcgs_info.symbol:
            index.append(name)
            A_list = []
            P_list = []
            for indv in time_list:
                ddd = data[(data['gene']==name)&(data['AP']=='AP_0')&(data['sample']==indv)]
                A_list.append(ddd['scaled_density'].to_numpy()[0])
                ddd = data[(data['gene']==name)&(data['AP']=='AP_9')&(data['sample']==indv)]
                P_list.append(ddd['scaled_density'].to_numpy()[0])
            A_data.append(A_list)
            P_data.append(P_list)
        A_array = np.array(A_data)
        linkage_data = linkage(A_array, metric='correlation') #optimal_ordering=True
        #dendrogram(linkage_data)
        #plt.savefig(f'{c}_Atree.pdf', format='pdf')
        A_array = pd.DataFrame(A_array, columns=time_list2, index=index)
        #A_array.to_csv(f'{c}_A_array.csv')
        sns.clustermap(A_array, vmin=-2, vmax=3, cmap="coolwarm", yticklabels=1, row_linkage=linkage_data, col_cluster=False)#, metric='correlation')
        plt.xticks(rotation=45)
        plt.title(c)
        plt.savefig(f'../{c}_Acluster.pdf', format='pdf')


        P_array = np.array(P_data)
        linkage_data = linkage(P_array, metric='correlation') #optimal_ordering=True
        #dendrogram(linkage_data)
        #plt.savefig(f'{c}_Ptree.pdf', format='pdf')
        P_array = pd.DataFrame(P_array, columns=time_list2, index=index)
        #P_array.to_csv(f'{c}_P_array.csv')
        sns.clustermap(P_array, vmin=-2, vmax=3, cmap="coolwarm", yticklabels=1, row_linkage=linkage_data, col_cluster=False)#, metric='correlation')
        plt.xticks(rotation=45)
        plt.title(c)
        plt.savefig(f'../{c}_Pcluster.pdf', format='pdf')


