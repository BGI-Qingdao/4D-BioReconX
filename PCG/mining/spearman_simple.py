"""
calculate spearman's cor
"""
import sys
import numpy as np
import pandas as pd
import scipy.stats
from util import get_max_of_dict, get_min_of_dict, read_list, read_expr
import matplotlib.pyplot as plt
from itertools import product
from filter_gene import *


if __name__ == '__main__':
    wt_expr = pd.read_csv(sys.argv[1])
    pcgs = read_list(sys.argv[2])

    target_genes = list(wt_expr.columns)
    print(f'target gene number: {len(target_genes)}')
    target_genes = remove_marker(target_genes)
    target_genes = remove_pcg(target_genes, pcgs)

    df_data = {'gene':[], 'pcg':[], 'scc':[]}
    for gene in target_genes:
        for pcg in pcgs:
            if pcg in wt_expr.columns and pcg != gene:
                df_data['gene'].append(gene)
                df_data['pcg'].append(pcg)
                scc = scipy.stats.spearmanr(wt_expr[gene], wt_expr[pcg]).correlation
                df_data['scc'].append(scc)


    results = pd.DataFrame(df_data)
    results = results.sort_values(by='scc',ascending=False)
    results.to_csv('spearman.csv', index=False)
    
