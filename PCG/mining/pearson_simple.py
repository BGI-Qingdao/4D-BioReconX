import pandas as pd
import numpy as np
import sys
from util import read_list, read_cluster_gene
from filter_gene import *


if __name__ == '__main__':
    wt_expr = pd.read_csv(sys.argv[1])
    pcgs = read_list(sys.argv[2])

    genes = list(wt_expr.columns)
    print(len(genes))
    genes = remove_marker(genes)
    genes = remove_pcg(genes,pcgs)
    df_data = {'gene':[], 'pcg':[], 'pcc':[]}
    for gene in genes:    
        for pcg in pcgs:
            if gene != pcg and pcg in wt_expr.columns:
                pcc = np.corrcoef(wt_expr[gene], wt_expr[pcg])[0,1]
                df_data['gene'].append(gene)
                df_data['pcg'].append(pcg)
                df_data['pcc'].append(pcc)

    results = pd.DataFrame(df_data)
    results = results.sort_values(by='pcc',ascending=False)
    results.to_csv('pearson.csv', index=False)

