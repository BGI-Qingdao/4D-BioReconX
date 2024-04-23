import sys
import getopt
import numpy as np
import pandas as pd

def usage():
    print("""
Usage: python3 s02.filter.py [-f filter threshold, default 0.05]
                                      [-o output folder, default ./]

we will load data generate by s02 in output folder!
""")

def main(argv):
    ########################
    # no args equal to -h
    if len(argv) == 0 :
        usage()
        sys.exit(0)
    prefixd = './'
    factor = 0.05
    # parse args
    try:
        opts, args = getopt.getopt(argv,"h:o:f:",["help"])

    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif opt in ("-f" ):
            factor = float(arg)
        elif opt in ("-o" ):
            prefixd = arg
    data = pd.read_csv(f'{prefixd}/Lines_anno.csv',sep='\t',header=0)
    from_counts = {}
    to_counts = {}
    from_names = np.unique(data['From'])
    to_names = np.unique(data['To'])
    for x in from_names:
        from_counts[x] = np.sum(data[data['From'] == x]['number'])

    for y in to_names:
        to_counts[y] = np.sum(data[data['To'] == y]['number'])

    def is_noise(row):
        if float(row['number'])/float(from_counts[row['From']]) < factor or float(row['number'])/float(to_counts[row['To']])  < factor:
            return False
        else :
            return True

    data['remain'] = data.apply(lambda row: is_noise(row), axis=1)
    data = data[data['remain'] == True].copy()
    data.to_csv(f'{prefixd}/strong.csv',header=True,sep='\t',index=None)

if __name__ == '__main__':
    main(sys.argv[1:])


