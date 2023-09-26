import sys
import numpy as np
import pandas as pd


def remove_marker(potential_genes):
    mfn = 'anno.v2.markers.topn.xls'
    marker = pd.read_csv(mfn, delimiter='\t')[['group', 'names']]
    marker.columns = ['group', 'gene']
    for m in marker.gene:
        if m in potential_genes:
            potential_genes.remove(m)
    return potential_genes


def remove_pcg(l,pcgs):
    if len(l) < len(pcgs):
        for i in l:
            if i in pcgs:
                l.remove(i)
    elif len(pcgs) < len(l):
        for i in pcgs:
            if i in l:
                l.remove(i)
    return l


def get_top_n(df_data:pd.DataFrame, topn:int, m:str, fn:str, ascending:bool=False):
    final = set()
    for ref in set(df_data['pcg']):
        data = df_data[df_data['pcg'] == ref]
        s = data.sort_values(by=[m],ascending=ascending)
        genes = list(s.head(topn).gene)
        for i in genes:
            final.add(i)
    ll = list(final)

    with open(fn,'w') as f:
        f.writelines('\n'.join(list(ll)))



if __name__ == '__main__':
    results = pd.read_csv(sys.argv[1])
    print(results.columns)
    top = int(sys.argv[2])
    method = sys.argv[3]
    get_top_n(results, top, method, f'{method}_top{top}.txt', False)
    #if sys.argv[4] == 1:
    #    get_top_n(results, top, method, f'{method}_top{top}.txt', True)
    #elif sys.argv[4] == 0:
    #    get_top_n(results, top, method, f'{method}_top{top}.txt', False)




