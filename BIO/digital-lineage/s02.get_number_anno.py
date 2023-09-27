import sys
import numpy as np
import pandas as pd

label = pd.read_csv('NeibourCell.label.xyz.txt',sep='\t',header=0)
label = label.set_index('cell')
sample = ['0hpa1','12hpa2','36hpa2','3dpa2','5dpa1','7dpa2','10dpa1','14dpa1']

anno_list = [  'parenchymal',
               'epidermal',
               'gut',
               'Nb2',
               'parenchymal',
               'gut',
               'cathepsin',
               'parenchymal',
               'gut',
               'parenchymal',
               'cathepsin',
               'parenchymal',
               'neural',
               'parenchymal',
               'neural',
               'epidermal',
               'cathepsin',
               'muscle',
               'epidermal',
               'pharynx',
               'epidermal',
               'cathepsin',
               'parenchymal',
               'parenchymal',
               'epidermal',
               'neural',
               'parenchymal',
               'muscle',
               'neural',
               'pharynx',
               'epidermal',
               'epidermal',
               'parenchymal',
               'parenchymal',
               'protonephridia',
               'neural',]

label['anno'] = label.apply(lambda row: anno_list[row['label']], axis=1)
maplist = []


def update_map(mapdata,source,target):
    mapdata['from_cell'] = mapdata.apply(lambda row :f'{row["from"].split(".")[1]}_{row["from"].split(".")[0]}',axis=1)
    mapdata['to_cell'] = mapdata.apply(lambda row :f'{row["to"].split(".")[1]}_{row["to"].split(".")[0]}',axis=1)
    mapdata['from_cluster'] = mapdata.apply(lambda row : label.loc[ row['from_cell'] ]['anno'],axis=1)
    mapdata['to_cluster'] = mapdata.apply(lambda row : label.loc[ row['to_cell'] ]['anno'],axis=1)
    # Set source_to as from
    mapdata['From'] = mapdata.apply(lambda row :f'{source}_{row["to_cluster"]}',axis=1)
    mapdata['To'] = mapdata.apply(lambda row :f'{target}_{row["from_cluster"]}',axis=1)
    return mapdata
    

def find_map_min(row,rev_check):
    ck = rev_check[rev_check['To'] == row['From']]
    ck = ck[ck['From'] == row['To']]
    if ck.empty:
        return 0
    elif ck['number'].to_list()[0] < row['number']:
        return ck['number'].to_list()[0]
    else :
        return row['number']

for i in [ 1,2,3,4,5,6,7]:
    target = sample[i]
    source = sample[i-1]
    prefix= f'{target}_{source}'
    mapdata = pd.read_csv(f'../s01.get_map/{prefix}.map.txt',sep='\t',header=0)
    mapdata = update_map(mapdata,source,target)
    valid_data = mapdata[['From','To']].copy()
    valid_data = valid_data.groupby(['From','To']).size().reset_index(name='number')

    target = sample[i-1]
    source = sample[i]
    prefix= f'{target}_{source}'
    revdata = pd.read_csv(f'../s01.get_map/{prefix}.map.txt',sep='\t',header=0)
    revdata = update_map(revdata,source,target)
    rev_check = revdata[['From','To']].copy()
    rev_check = rev_check.groupby(['From','To']).size().reset_index(name='number')
    valid_data['number'] =  valid_data.apply( lambda row: find_map_min(row,rev_check),axis=1)

    maplist.append(valid_data)


final = pd.concat(maplist,ignore_index=True)
final.to_csv('Lines_anno.csv',sep='\t',header=True,index=False)
