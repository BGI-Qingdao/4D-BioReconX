import sys
import numpy as np
import pandas as pd
import anndata as ad

for x in ['0hpa1','12hpa2','36hpa2','3dpa2','5dpa1','7dpa2','10dpa1','14dpa1']:
    for y in ['0hpa1','12hpa2','36hpa2','3dpa2','5dpa1','7dpa2','10dpa1','14dpa1']:
        if x == y :
            continue
        indata1 = f'../s00.indata/{x}.NsC.sct.scale.data.h5ad'
        indata2 = f'../s00.indata/{y}.NsC.sct.scale.data.h5ad'
        prefix = f'{x}_{y}'
                
        data1 = ad.read_h5ad(indata1)
        data2 = ad.read_h5ad(indata2)
        
        data1.obs['repcell'] =  data1.obs.apply(lambda row : row['hood'].split(',')[0],axis=1)
        data2.obs['repcell'] =  data2.obs.apply(lambda row : row['hood'].split(',')[0],axis=1)
        
        
        common_gene = np.intersect1d(data1.var.index.to_numpy(),data2.var.index.to_numpy())
        print(len(common_gene))
        data1 = data1[:, common_gene]
        print(len(data1.var))
        data2 = data2[:, common_gene]
        print(len(data2.var))
        
        cell1 = data1.obs['repcell'].to_numpy()
        cell2 = data2.obs['repcell'].to_numpy()
        
        
        R =  np.corrcoef(data1.X , data2.X)
        R = R[len(cell1):,:len(cell1)] # y as cell2 and x as cell1
        targets = np.argmax(R,axis=0)
        
        ret = pd.DataFrame()
        ret['from'] = cell1
        ret['to'] = cell2[targets]
        
        ret.to_csv(f'{prefix}.map.txt',sep='\t',header=True,index=False)
