#!/usr/bin/env python3

import sys
import getopt
import numpy as np
import pandas as pd
import scipy.ndimage as nd
import skimage.morphology as sm
from skimage import io,filters
from pathlib import Path
import cv2

import clahe_enhance

def get_mask(ssdna_file,radius,prefix):
    ssdna = io.imread(ssdna_file)

    if len(ssdna.shape) == 3:  # RGB tiff to 8 bit gray tiff
        new_data = np.zeros((ssdna.shape[0], ssdna.shape[1]), dtype=int)
        new_data = new_data + ssdna[:, :, 0]
        new_data = new_data + ssdna[:, :, 1]
        new_data = new_data + ssdna[:, :, 2]
        ssdna = new_data.astype('uint8')

    # ssdna =sfr.enhance_contrast(ssdna, sm.disk(9))
    # ssdna =sfr.enhance_contrast(ssdna, sm.disk(3))

    #convert to mask
    thre = filters.threshold_otsu(ssdna)
    ssdna[ssdna >= thre] = 255
    ssdna[ssdna < thre] = 0

    #expand pixel
    ssdna = sm.dilation(ssdna, sm.disk(radius))
    # ssdna_dilation=sm.dilation(ssdna_dilation,sm.square(5))
    # edges = edges.astype('uint8')
    # print(edges)
    # edges = edges.astype('uint8')'''

    io.imsave(f'{prefix}.mask.tif', ssdna)
    return ssdna


from skimage.filters import rank
from skimage.morphology import disk, ball
from skimage.restoration import denoise_tv_chambolle

def remove_outliers(image, radius=2, threshold=50):
    #footprint_function = disk if image.ndim == 2 else ball
    #footprint = footprint_function(radius=radius)
    kernel = np.ones((radius, radius), np.uint8)
    median_filtered = rank.median(image, kernel)
    outliers = (
            (image > median_filtered + threshold)
            | (image < median_filtered - threshold)
            )
    output = np.where(outliers, median_filtered, image)
    return output

def cvthreshold(ssDNA, prefix):
    image = cv2.imread(ssDNA, 0)

    #threshold = np.quantile(image, .75) + 10
    #threshold = np.mean(image)
    #ret, th = cv2.threshold(image, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    #th = cv2.adaptiveThreshold(image, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, 11, 2)
    #print(f'threshold cutoff {threshold}')
    ret, th = cv2.threshold(image, 20, 255, cv2.THRESH_BINARY)

    #kernel = np.ones((8, 8), np.uint8)
    #th = cv2.erode(th, kernel, iterations = 1)
    #cv2.imwrite(f'{prefix}.erode.tif', th)

    #th = cv2.morphologyEx(th, cv2.MORPH_CLOSE, kernel)
    #cv2.imwrite(f'{prefix}.close.tif', th)
    #th = cv2.medianBlur(th, 31)
    #cv2.imwrite(f'{prefix}.mb.tif', th)
    #th = remove_outliers(th, 31)
    #cv2.imwrite(f'{prefix}.ro.tif', th)

    cv2.imwrite(f'{prefix}.mask.tif', th)
    return th

ssDNA = sys.argv[1]
prefix = Path(ssDNA).name.replace('.tif', '')
print(prefix)

#get_mask(sys.argv[1], 10, prefix)
cvthreshold(ssDNA, prefix)

