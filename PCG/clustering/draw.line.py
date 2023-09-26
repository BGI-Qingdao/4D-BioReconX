#!/usr/bin/env python3
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'



def scale(a):
    return  (a - np.mean(a)) / np.std(a)

data={}
pcgs_info = np.loadtxt(sys.argv[1], delimiter='\t',dtype=str)

for indv in ['WT','0hpa1','12hpa2','36hpa2','3dpa2','5dpa1','7dpa2','10dpa1','14dpa1']:
    density_file = f'{indv}_bin10/{indv}_10.density.txt'
    density_indv = pd.read_csv(density_file,sep='\t',header=0)
    data[indv] = density_indv

pd_list = []
for smes, name, c in pcgs_info:
    print(name)
    notskip=False
    for indv in ['WT','0hpa1','12hpa2','36hpa2','3dpa2','5dpa1','7dpa2','10dpa1','14dpa1']:
        if smes in data[indv].columns:
            print(f'{smes} is in {indv} file, not skip')
            notskip=True
            break
        else:
            print(f'{smes} is not in {indv} file, skip')
    if notskip == False:
        continue
    gene_pd_list = []
    for indv in ['WT','0hpa1','12hpa2','36hpa2','3dpa2','5dpa1','7dpa2','10dpa1','14dpa1']:
        print(f'{indv}',flush=True)
        curr_gene_pd = pd.DataFrame()
        curr_gene_pd['AP'] = data[indv]['AP'].to_numpy()
        if smes not in data[indv].columns:
            curr_gene_pd['density'] = np.zeros(len(curr_gene_pd))
        else :
            curr_gene_pd['density'] = data[indv][smes].to_numpy()
        curr_gene_pd['gene'] = [name]*len(curr_gene_pd)
        curr_gene_pd['sample'] = [indv]*len(curr_gene_pd)
        print(curr_gene_pd)
        gene_pd_list.append(curr_gene_pd)
    gene_pd = pd.concat(gene_pd_list,ignore_index=True)
    gene_pd['scaled_density'] = scale(gene_pd['density'])
    pd_list.append(gene_pd)

all_pd = pd.concat(pd_list,ignore_index=True)


cm = 1/2.54
plt.figure(figsize=(3.3*cm,26*cm))
sns.relplot(data=all_pd,x='AP',y='scaled_density', col='gene', hue='sample',col_wrap=5, kind="line")
plt.ylim(-1,4)
plt.tight_layout()
plt.savefig('5pcg.pdf', format='pdf')

