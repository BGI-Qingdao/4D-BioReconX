# 2023-03-01


import pandas as pd
import numpy as np
from scipy.stats import entropy
import sys
from util import get_max_of_dict, get_min_of_dict, read_list, read_expr
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from filter_gene import *


wt_expr = pd.read_csv(sys.argv[1])
pcgs = read_list(sys.argv[2])

wt_expr = wt_expr + 3 + 1e-20
target_genes = list(wt_expr.columns)
print(len(target_genes))
target_genes = remove_pcg(target_genes,pcgs)
target_genes = remove_marker(target_genes)
print(len(target_genes))


df_data = {'gene':[], 'pcg':[], 'kl':[]}
for gene in target_genes:
    for pcg in pcgs:
        if pcg in wt_expr.columns and pcg != gene:
            qk = gaussian_filter1d(wt_expr[gene],sigma=3)
            pk = gaussian_filter1d(wt_expr[pcg],sigma=3)
            kl = entropy(qk=qk,pk=pk)
            df_data['gene'].append(gene)
            df_data['pcg'].append(pcg)
            df_data['kl'].append(kl)

results = pd.DataFrame(df_data)
results.to_csv('kl.csv', index=False)

