import numpy as np
import pandas as pd

data = pd.read_csv('Lines_anno.csv',sep='\t',header=0)

from_counts = {}
to_counts = {}

from_names = np.unique(data['From'])
to_names = np.unique(data['To'])

for x in from_names:
    from_counts[x] = np.sum(data[data['From'] == x]['number'])

for y in to_names:
    to_counts[y] = np.sum(data[data['To'] == y]['number'])


def is_noise(row):
    if float(row['number'])/float(from_counts[row['From']]) < 0.05 or float(row['number'])/float(to_counts[row['To']])  < 0.05:
        return False
    else :
        return True

data['remain'] = data.apply(lambda row: is_noise(row), axis=1)

data = data[data['remain'] == True].copy()

data.to_csv('strong.csv',header=True,sep='\t',index=None)
