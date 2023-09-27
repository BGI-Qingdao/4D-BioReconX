import numpy as np
import pandas as pd
import scipy.ndimage as nd
from skimage import io as skio

import sys

SP=sys.argv[1]

binsize = 20.0
infos = pd.read_csv(f'{SP}.csv',sep=' ',header=None)
infos.columns = ['filename','zvalue']
infos['zvalue'] = infos['zvalue'].astype(int)
infos['zvalue'] = infos['zvalue'] - 1 # start from 0
slices = {}

annos = pd.read_csv(f'../01.prepare_anno/{SP}.cellbin.SPC.txt',sep=',',header=0)

for i , row in infos.iterrows():
    cellmask = np.loadtxt(row['filename'],dtype=int)
    y, x = np.nonzero(cellmask)
    tmp_draw = pd.DataFrame()
    tmp_draw['x'] = x
    tmp_draw['y'] = y
    tmp_draw['cell'] = cellmask[y,x]
    cellmask[y,x] = 0
    tmp_draw= tmp_draw.set_index('cell')
    slice_anno = annos[ annos['slice_id'] == int(row['zvalue']+1) ].copy()
    slice_anno = slice_anno.set_index('cell_id')
    tmp_draw['anno'] = slice_anno['anno_id']
    tmp_draw = tmp_draw[tmp_draw['anno']!='NA'].copy()
    cellmask[tmp_draw['y'],tmp_draw['x']] = 100
    pharynx = tmp_draw[tmp_draw['anno']=='c21']
    gut = tmp_draw[tmp_draw['anno']=='c1']
    neural = tmp_draw[tmp_draw['anno']=='c33']
    cellmask[pharynx['y'],pharynx['x']] = 150
    cellmask[gut['y'],gut['x']] = 200
    cellmask[neural['y'],neural['x']] = 250

    h,w = cellmask.shape
    affine = np.matrix(np.array([[1.0/binsize,0,0],[0,1.0/binsize,0],[0,0,1]]))
    binimage = nd.affine_transform(cellmask.T,affine.I,output_shape=(int(w/binsize),int(h/binsize)),order=0)
    slices[row['zvalue']] = binimage.T
    H = int(h/binsize)
    W = int(w/binsize)

zmax = infos['zvalue'].max()
image_buff = np.zeros((zmax+1,H,W),dtype='uint8')
for x in slices:
    image_buff[x,:,:] = slices[x]

skio.imsave(f'{SP}.tif',image_buff)

