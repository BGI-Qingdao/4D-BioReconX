import sys
import numpy as np
import pandas as pd
import anndata as ad

def usage():
    print("""
Usage: python3 s02.get_number_anno.py  <-i in.lst>
                                      [-o output folder, default ./]

example of in.lst

```
time0
time1
time2
time3
```
we will load time0.h5ad, time1.h5ad, time2.h5ad, time3.h5ad.
and data generate by s01 in output folder!

""")
def update_map(mapdata,source,target,s_adata,t_adata,anno_key):
    mapdata['from_cluster'] = mapdata.apply(lambda row : s_adata.obs.loc[ row['from'] ][anno_key],axis=1)
    mapdata['to_cluster'] = mapdata.apply(lambda row : t_adata.obs.loc[ row['to'] ][anno_key],axis=1)
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

def main(argv):
    ########################
    # no args equal to -h
    if len(argv) == 0 :
        usage()
        sys.exit(0)
    inf = ''
    prefixd = './'
    # parse args
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["help"])

    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif opt in ("-i" ):
            inf= arg
        elif opt in ("-o" ):
            prefixd = arg
    sample = np.loadtxt(inf,dtype=str).to_list()
    maplist = []
    for i in range(1,len(inlist)):
        target = sample[i]
        source = sample[i-1]
        prefix= f'{target}_{source}'
        mapdata = pd.read_csv(f'{prefixd}/{prefix}.map.txt',sep='\t',header=0)
        s_adata = ad.read_h5ad(f'{source}.h5ad')
        t_adata = ad.read_h5ad(f'{target}.h5ad')
        mapdata = update_map(mapdata,source,target,s_adata,t_adata,anno_key)
        valid_data = mapdata[['From','To']].copy()
        valid_data = valid_data.groupby(['From','To']).size().reset_index(name='number')

        target = sample[i-1]
        source = sample[i]
        prefix= f'{target}_{source}'
        revdata = pd.read_csv(f'{prefixd}/{prefix}.map.txt',sep='\t',header=0)
        revdata = update_map(revdata,source,target,t_adata,s_adata,anno_key)
        rev_check = revdata[['From','To']].copy()
        rev_check = rev_check.groupby(['From','To']).size().reset_index(name='number')
        valid_data['number'] =  valid_data.apply( lambda row: find_map_min(row,rev_check),axis=1)
        maplist.append(valid_data)

    final = pd.concat(maplist,ignore_index=True)
    final.to_csv(f'{prefixd}/Lines_anno.csv',sep='\t',header=True,index=False)

if __name__ == '__main__':
    main(sys.argv[1:])

