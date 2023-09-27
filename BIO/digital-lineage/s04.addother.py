import pandas as pd
import numpy as np

strong_data = pd.read_csv('strong.csv',sep='\t',header=0)
spc_numer = pd.read_csv('SPC.num.txt',sep='\t',header=0)

def prev_key(key):
    indv_list = [ 'WT', '0hpa1', '12hpa2', '36hpa2', '3dpa2', '5dpa1', '7dpa2','10dpa1','14dpa1' ]
    prev_map = { '0hpa1' : 0, '12hpa2':1, '36hpa2':2, '3dpa2':3, '5dpa1':4, '7dpa2':5,'10dpa1':6,'14dpa1':7 }
    indv = key.split('_')[0]
    prev0 = indv_list[prev_map[indv]]
    prev = f'{prev0}_other'   
    return prev

def checkfull(key,count):
    left = 0
    right = 0
    temp = strong_data[ strong_data['To'] == key ]
    if not temp.empty:
        left = np.sum(temp['number'])
    draw = left
    if draw < count :
        return prev_key(key), key, count-draw;
    else :
        return '', key ,count-draw

froms = []
tos = []
counts = []
for index, row in spc_numer.iterrows():
    ret = checkfull(row['key'],row['SPC_num'])
    if ret[0] != '':
        froms.append(ret[0])
        tos.append(ret[1])
        counts.append(ret[2])

rets = pd.DataFrame()
rets['From'] = froms
rets['To'] = tos
rets['number'] = counts

rets.to_csv('other.txt',sep='\t',header=True,index=False) 

strong_data = strong_data.drop(columns = ['remain'])
draw_data =  pd.concat([ rets, strong_data],ignore_index=True)

draw_data.to_csv('draw.txt',sep='\t',header=True,index=False)
