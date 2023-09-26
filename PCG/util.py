# 2022-05-27
# handle files
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def read_expr_(fn):
    expr = pd.read_csv(fn)
    expr = expr.dropna(how='all',axis=1)
    for c in expr.columns:
        if (expr[c]==0).all():
            expr = expr.drop(c,axis=1)
    return expr


def read_expr(fn):
    expr = pd.read_csv(fn)
    expr = expr.dropna(how='all',axis=1)
    for c in expr.columns:
        if (expr[c]==0).all():
            expr = expr.drop(c,axis=1)
    if expr.isna().sum().sum() > 0:
            print('this data contains NaN value!')
    return expr


def read_list(fn:str)->list:
    """extract list from a txt file where each line is one element of the list"""
    with open(fn, 'r') as f:
        list = f.read().splitlines()
    list = [i.strip() for i in list]
    return list


def read_symbols(fn:str,sep:str)->dir:
    symbol = {}
    with open(fn, "r") as f:
        for line in f.readlines():
            line = line.strip().split(sep)
            if len(line) != 2:
                continue
            else:
                symbol[line[0]] = line[1]
    return symbol


def get_symbol(gid:str, symbols:dict)->str:
    try:
        return symbols[gid]
    except KeyError:
        return gid


def read_gene_cluster(fn)->pd.DataFrame:
    """Load gene_cluster.txt file"""
    df = pd.read_csv(fn, delimiter='\t', header=None)
    time = df[0].str.split("_", expand=True)
    df[2] = time[0]
    df[3] = time[1]
    df.columns = ['id', 'cluster', 'gene', 'time']
    return df


def read_cluster_gene(fn:str)->dir:
    """Load previous saved cluster-gene.csv file"""
    cg = {}
    with open(fn, 'r') as f:
        for line in f.readlines():
            line = line.strip().split(',')
            cg[int(line[0])] = line[1:]
    return cg 


def plot_gene(genes:list, expr:pd.DataFrame, path='.'):
    for i, g in enumerate(genes):
        try:
            y = expr[g]
            plt.plot(y)
            plt.title(g)
            plt.savefig(f'{path}/{g}.png')
            plt.close()
        except IndexError:
            print(f'{g} not in wt data')


def get_max_of_dict(d):
    """max value of a flat dictionary"""
    maxv = max(d.values())
    maxk = [key for key, value in d.items() if value == maxv]
    return maxk, maxv


def get_min_of_dict(d):
    minv = min(d.values())
    mink = [key for key, value in d.items() if value == minv]
    return mink, minv

   
