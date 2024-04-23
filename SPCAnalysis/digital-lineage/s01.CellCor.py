import sys
import getopt
import numpy as np
import pandas as pd
import anndata as ad

def usage():
    print("""
Usage   : python3 s01.CellCor.py  <-i in.lst>
                                  [-o output folder, default ./]

example of in.lst

```
time0
time1
time2
time3
```
and we will load time0.h5ad, time1.h5ad, time2.h5ad, time3.h5ad.

""")

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
    inlist = np.loadtxt(inf,dtype=str)
    for x in inlist:
        for y in inlist:
            if x == y :
                continue
            indata1 = f'{x}.h5ad'
            indata2 = f'{y}.h5ad'
            prefix = f'{x}_{y}'
            data1 = ad.read_h5ad(indata1)
            data2 = ad.read_h5ad(indata2)
            common_gene = np.intersect1d(data1.var.index.to_numpy(),data2.var.index.to_numpy())
            print(len(common_gene))
            data1 = data1[:, common_gene]
            print(len(data1.var))
            data2 = data2[:, common_gene]
            print(len(data2.var))
            cell1 = data1.index.to_numpy()
            cell2 = data2.index.to_numpy()

            R =  np.corrcoef(data1.X , data2.X)
            R = R[len(cell1):,:len(cell1)] # y as cell2 and x as cell1
            targets = np.argmax(R,axis=0)

            ret = pd.DataFrame()
            ret['from'] = cell1
            ret['to'] = cell2[targets]
            ret.to_csv(f'{prefixd}/{prefix}.map.txt',sep='\t',header=True,index=False)

if __name__ == '__main__':
    main(sys.argv[1:])
