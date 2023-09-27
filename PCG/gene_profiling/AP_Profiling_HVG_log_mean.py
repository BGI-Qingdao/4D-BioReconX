import scanpy as sc
import numpy as np
import pandas as pd
import sys
import argparse
import scipy
import anndata as ad

def ParseArgs():
    parser = argparse.ArgumentParser(
                 prog='AP_Profiling',
                 description='Create AP Profiling for genes.',
                 epilog='Best wishes'
                 )
    parser.add_argument('-i', '--input', type=str, required=True, help='the input adata, Notice: we expect raw count matrix!!!')
    parser.add_argument('-o', '--output', type=str, required=True, help = 'the output prefix')
    parser.add_argument('-m', '--min_bin', default=50, type=int, required=False, help='the minimum allowed cell number of a bin, default(50)')
    parser.add_argument('-n', '--num_HVG', default=5000, type=int, required=False, help='the number of HVG bin, default(5000)')
    parser.add_argument('-b', '--bin_num', default=100, type=int, required=False, help='the total bin number, default(100)')
    parser.add_argument('-s', '--spatial', default='spatial', type=str, required=False, help='the coordinate key in obsm default(spatial)')
    parser.add_argument('-a', '--axis', default=0, type=int, required=False, help='idx of spatial default(0)')
    args = parser.parse_args()
    print(f'input = {args.input}')
    print(f'output = {args.output}')
    print(f'min #cell in bin = {args.min_bin}')
    print(f'#bin = {args.bin_num}')
    return args

def GetValidConf(xPos, bin_num, min_bin):
    xmin, xmax = np.min(xPos), np.max(xPos)+1
    tmp_xmin, tmp_xmax = xmin, xmax
    while(True):
        binsize = (float(tmp_xmax)-float(tmp_xmin))/float(bin_num)
        ncell_first_bin = np.sum(((xPos>=tmp_xmin) & (xPos < tmp_xmin+binsize)))
        if ncell_first_bin < min_bin:
            tmp_xmin = tmp_xmin + 1
            continue
        ncell_last_bin = np.sum(((xPos>=tmp_xmax-binsize) & (xPos < tmp_xmax)))
        if ncell_last_bin < min_bin:
            tmp_xmax = tmp_xmax - 1
            continue
        break
    xmin, xmax = tmp_xmin, tmp_xmax
    binsize = (float(xmax)-float(xmin))/float(bin_num)
    print(f'final xmin = {xmin}')
    print(f'final xmax = {xmax}')
    print(f'final binsize = {binsize}')
    return xmin, xmax, binsize

def CreateBinLabels(xPos,xmin, xmax, binsize, bin_num):
    ret = np.zeros(xPos.shape)
    ret[:] = -1 # -1 refer to outliers
    for binid in range(bin_num):
        now_start = xmin + binid*binsize
        now_end = xmin + (binid+1)*binsize
        ret[((xPos>=now_start)&(xPos<now_end))]=binid
    return ret

def GetHVG(adata,nHVG):
    tmp_adata = adata.copy()
    sc.pp.normalize_total(tmp_adata)
    sc.pp.log1p(tmp_adata)
    ret = sc.pp.highly_variable_genes(tmp_adata,n_top_genes=nHVG,inplace=False)
    return tmp_adata.var.index[ret['highly_variable']].to_list()

def CreatePseudoBulk(adata,cell_labels,hvg_list,bin_num):
    adata.obs['bulkid'] = cell_labels
    adata = adata[adata.obs['bulkid']!=-1].copy()
    # norm total
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    # filter hvg
    adata = adata[:,hvg_list].copy()
    ret = pd.DataFrame()
    cell_in_bin = {}
    #for binid in adata.obs['bulkid'].unique():
    for binid in range(bin_num):
        bin_adata = adata[adata.obs['bulkid']==binid].copy()
        X = bin_adata.X
        if scipy.sparse.issparse(X):
            X = X.todense()
        mean = np.mean(X,axis=0)
        mean = np.array(mean)[0]
        ret[f'bin{binid}'] = mean
        cell_in_bin[f'bin{binid}'] = X.shape[0]
    ret_adata = ad.AnnData(X=ret.to_numpy().T,dtype=float,
                var=adata.var)
    ret_adata.obs_names = ret.columns
    return ret_adata, cell_in_bin

def main():
    ################################ Parse args
    args = ParseArgs()
    ################################ Load data
    adata = sc.read_h5ad(args.input)
    ################################ Get position data
    xPos = adata.obsm[args.spatial][:,args.axis]
    ################################ Get bin system
    xmin, xmax, binsize = GetValidConf(xPos,args.bin_num, args.min_bin)
    ################################ Create bin lables
    cell_labels = CreateBinLabels(xPos,xmin, xmax, binsize, args.bin_num)
    ################################ Create pseudo bulk adata
    hvg_genes = GetHVG(adata,args.num_HVG)
    rdata, cell_in_bin = CreatePseudoBulk(adata,cell_labels,hvg_genes,args.bin_num)
    rdata.uns['xmin'] = xmin
    rdata.uns['xmax'] = xmax
    rdata.uns['binsize'] = binsize
    rdata.uns['ncell_in_bin'] = cell_in_bin
    ################################ Save result
    rdata.write(f'{args.output}.h5ad')

if __name__ == '__main__':
    main()
