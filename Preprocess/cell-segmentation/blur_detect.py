#!/usr/bin/env python3

import os
import sys

import numpy as np
import cv2
from skimage import color, data, restoration
from scipy.signal import convolve2d as conv2

def color_masking(image):
    # background was detected as blur region
    # I need to segmentate tissue region firstly
    # here I used color masking for segmentation on green channel
    b, g, r = cv2.split(image)

    # detect based on green channel
    light_green = 10
    dark_green = 255
    mask = cv2.inRange(g, light_green, dark_green)
    kernel = np.ones((10, 10), np.uint8)
    mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, kernel)
    #result = cv2.bitwise_and(image, image, mask=mask)
    return mask

def blur_detect(image, chunk_size=3, method='laplacian', top_svd=30):
    blur_image = np.zeros(shape=image.shape, dtype=np.uint8)
    for (x, y), value in np.ndenumerate(mask):
        if value == 0:
            continue
        chunk = image[x:x+chunk_size, y:y+chunk_size]
        # small value indicate blur region
        if method == 'laplacian':
            blur_value = cv2.Laplacian(chunk, cv2.CV_64F).var()
        elif method == 'svd':
            u, sigma, vt = np.linalg.svd(img)
            blur_value = sum(sigma[:top_svd]) / sum(sigma)
        blur_image[x, y] = blur_value
        #print(blur_value)
    return blur_image

if len(sys.argv[1:]) != 2:
    sys.stderr.write('usage: {} image prefix'.format(__file__))
    sys.exit()

image = cv2.imread(sys.argv[1])
prefix = sys.argv[2]

mask = color_masking(image)
#cv2.imwrite('masked.png', masked_image)

image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

blur_image = blur_detect(image, chunk_size=3, method='laplacian')
#mask_inv = 1 - (mask // 255)
#blur_image[mask_inv] = 0
np.savetxt('{}_blurValue.txt'.format(prefix), blur_image, fmt='%d')

blur_rgb_image = cv2.applyColorMap(blur_image, cv2.COLORMAP_JET)
cv2.imwrite('{}_blurHeatmap.png'.format(prefix), blur_rgb_image)

black = np.zeros(shape=image.shape, dtype=np.uint8)
blur_mask = np.where(blur_image < 30, mask, black)
#kernel = np.ones((30, 30), np.uint8)
#blur_mask = cv2.morphologyEx(blur_mask, cv2.MORPH_OPEN, kernel)
cv2.imwrite('{}_blurMask.png'.format(prefix), blur_mask)

sys.exit()

stains = np.where(blur_image < 30, mask, image)
cv2.imwrite('stains.png', stains)

#inpaint = cv2.inpaint(image, blur_mask, 1, cv2.INPAINT_NS) # or cv2.INPAINT_TELEA
inpaint = cv2.inpaint(image, blur_mask, 1, cv2.INPAINT_TELEA) # or cv2.INPAINT_TELEA
cv2.imwrite('inpaint.png', inpaint)

psf = np.ones((5, 5)) / 25
image = conv2(image, psf, 'same')

image += 0.1 * image.std() * np.random.standard_normal(image.shape)

#deconvolved = restoration.wiener(image, psf, 1, clip=False)
#deconvolved, _ = restoration.unsupervised_wiener(image, psf)
deconvolved = restoration.richardson_lucy(image, psf, iterations=10)
cv2.imwrite('deconv.png', deconvolved)



