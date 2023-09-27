"""usage: {} --cell <cell_mask> --tissue <tissue_mask> --blur <blur_mask> --image <orig_img.tif>
"""

import os
import sys
import getopt
import copy

from collections import Counter

import numpy as np
import pandas as pd
import cv2

def getfootprint(struc, a, b=None):
    from skimage.morphology import (
            square, 
            rectangle, 
            diamond, 
            disk, 
            octagon, 
            star)

    struc_lib = {
            'square': square, 
            'rectangle': rectangle, 
            'diamond': diamond, 
            'disk': disk, 
            'octagon': octagon, 
            'star': star
            }

    morph = struc_lib[struc]

    if struc in ['rectangle', 'octagon']:
        if b is None:
            sys.stderr.write('two args required\n')
            sys.exit()
        return morph(a, b)
    else:
        if b is not None:
            sys.stderr.write('only one arg required\n')
            sys.exit()
        return morph(a)

class Stoarr:
    def __init__(self, matrix):
        if isinstance(matrix, str):
            if matrix.endswith('.txt'):
                matrix = np.loadtxt(matrix)
            elif matrix.endswith(('.tif', '.png')):
                matrix = cv2.imread(matrix, cv2.IMREAD_UNCHANGED)
        self.matrix = matrix.astype(int)
        
    def to_triplet(self, name='mask'):
        import scipy.sparse

        mtx= scipy.sparse.csc_matrix(self.matrix)
        mtx = mtx.tocoo()
        tmp = []
        for x, y, mask in zip(mtx.row, mtx.col, mtx.data):
            tmp.append([x, y, int(mask)])
        triplet = pd.DataFrame(tmp, columns=['x', 'y', name])
        return triplet

    def binning(self, bin_size):

        sys.stdout.write('binning ... ')
        sys.stdout.flush()

        triplet = self.to_triplet()
    
        triplet['xbin'] = (triplet.x / bin_size).astype(int) * bin_size
        triplet['ybin'] = (triplet.y / bin_size).astype(int) * bin_size
        triplet['bin'] = triplet.xbin.astype(str) + '_' + triplet.ybin.astype(str)

        index = [(-i, x) for i, x in enumerate(triplet['bin'].unique())]
        index = pd.DataFrame(index, columns=['N', 'bin'])

        triplet = triplet.merge(index, how='left', on='bin')
    
        matrix = np.zeros(shape=self.matrix.shape, dtype=int)
        matrix[triplet['x'], triplet['y']] = triplet['N']

        sys.stdout.write('done\n')
        return Stoarr(matrix)

    def to_binary(self):
        obj = copy.deepcopy(self)
        mask = np.isin(obj.matrix, [0], invert=True)
        obj.matrix[mask] = 1
        return obj

    def subtract(self, other):

        sys.stdout.write('subtracting ... ')
        sys.stdout.flush()

        obj = copy.deepcopy(self)

        obj = obj.to_binary()
        other = other.to_binary()

        obj.matrix = obj.matrix - other.matrix

        sys.stdout.write('done\n')
        return obj
    
    def intersection(self, other, label_area_cutoff=0.3):
        """intersection of label mask and binary mask
        * mask: binary matrix
        * label_area_cutoff: labels with greater area will be dropped
        """
        
        sys.stdout.write('intersection ... ')
        sys.stdout.flush()

        obj = copy.deepcopy(self)

        if isinstance(other, Stoarr):
            other = other.to_binary()

        values = np.unique(obj.matrix)
        if len(values) == 2:
            mask = cv2.bitwise_and(obj.matrix, other.matrix)
            mask = np.invert(mask.astype(bool))
        else:
            binary = self.to_binary()
            
            mask = cv2.bitwise_and(binary.matrix, other.matrix)
            mask = np.invert(mask.astype(bool))

            orig_counter = Counter(obj.matrix.flatten())

            filter_part = obj.matrix[mask]
            filter_counter = Counter(filter_part.flatten())

            filter_labels = []
            for label, pixels in filter_counter.items():
                if label == 0:
                    continue
                ratio = pixels / orig_counter[label]
                if ratio < label_area_cutoff:
                    continue
                filter_labels.append(label)

            filter_labels = list(set(filter_labels))
            mask = np.isin(obj.matrix, filter_labels)

        obj.matrix[mask] = 0

        sys.stdout.write('{} labels removed\n'.format(len(filter_labels)))
        return obj

    def relabel(self, label_map=None):
        if label_map is None:
            unique_labels, labels = np.unique(self.matrix, return_inverse=True)
            matrix = labels.reshape(self.matrix.shape)

            #obj = Mask(matrix)
            #obj.unique_labels = unique_labels
            #obj.labels = labels
            return Stoarr(matrix)
        else:
            triplet = self.to_triplet()

            triplet = triplet.merge(label_map, how='left', 
                    left_on='mask', right_index=True)
    
            matrix = np.zeros(shape=self.matrix.shape, dtype=int)
            matrix[triplet['x'], triplet['y']] = triplet['mask_y']
            return Stoarr(matrix)
    
    def retrieve(self):
        if not self.unique_labels and not self.labels:
            return

        matrix = self.unique_labels[self.labels]
        matrix = matrix.reshape(self.shape)
        obj = Stoarr(matrix)

        return obj

    def minimum_filter(self, footprint='octagon', ksize=(4, 4), iterations=2):

        sys.stdout.write('minimum filter ... ')
        sys.stdout.flush()

        obj = copy.deepcopy(self)
        obj.matrix = obj.matrix.astype(np.uint8)
        #obj.matrix = cv2.applyColorMap(
        #        obj.matrix,
        #        cv2.COLORMAP_JET
        #        )
        
        try:
            n, m = ksize
        except:
            n = ksize
            m = None
        footprint = getfootprint(footprint, n, m)
        obj.matrix = cv2.erode(
                obj.matrix, 
                kernel=footprint, 
                iterations=iterations
                )
        #cv2.imwrite('blur.png', obj.matrix)

        sys.stdout.write('done\n')
        return obj
    
    def filter_by_matrix(self, on=None, min_value=None, max_value=None, 
            draw=False, prefix=None):
        """label mask method
        * on: filter by minimum value of the input matrix
        """
        
        sys.stdout.write('filter by matrix ... ')
        sys.stdout.flush()

        obj = copy.deepcopy(self)

        triplet = obj.to_triplet()
        ref = on.to_triplet()
        triplet = triplet.merge(ref, how='left', on=('x', 'y'))
        triplet = triplet.fillna(0)

        medians = triplet.groupby('mask_x')['mask_y'].median()
        medians = medians.to_frame()

        if draw:
            fig = self.relabel(medians)
            cv2.imwrite(f'{prefix}.median.png', fig.matrix)

        if min_value:
            filter_labels = medians[medians['mask_y'] < min_value].index.values
        if max_value:
            filter_labels = medians[medians['mask_y'] > max_value].index.values
        
        mask = np.isin(obj.matrix, filter_labels)
        obj.matrix[mask] = 0

        sys.stdout.write('{} labels removed\n'.format(len(filter_labels)))
        return obj

    def filter_by_diameter(self, min_size=1, max_size=None):
        """label mask method
        * min_size: max circo radius
        """

        sys.stdout.write('filter by diameter ... ')
        sys.stdout.flush()

        from skimage.measure import regionprops

        obj = copy.deepcopy(self)
        #obj.matrix = obj.matrix.astype(np.uint8)
        
        filter_labels = []
        regions = regionprops(obj.matrix)
        for index, props in enumerate(regions):
            if props.minor_axis_length <= 8 and (props.minor_axis_length * 5 
                    <= props.major_axis_length):
                # abnormity cell with large aspect ratio
                filter_labels.append(props.label)
                continue
            if props.area > 1000 or props.area < 6:
                # extreme large cell caused by non-detected blur region
                # extreme small cell original segmentation fault
                filter_labels.append(props.label)
                continue
            if props.extent < 0.3:
                filter_labels.append(props.label)
                continue
            if props.minor_axis_length < min_size:
                # extreme thin cell
                filter_labels.append(props.label)
                continue
            if max_size and props.major_axis_length > max_size:
                # extreme fat cell
                filter_labels.append(props.label)
                continue

        mask = np.isin(obj.matrix, filter_labels)
        obj.matrix[mask] = 0

        sys.stdout.write('{} labels removed\n'.format(len(filter_labels)))
        return obj

    def merge(self, other, how='left'):
        
        sys.stdout.write('merge mix labels ... ')
        sys.stdout.flush()

        if how == 'left':
            obj = copy.deepcopy(self)
            mask1 = obj.to_binary()
            mask2 = copy.deepcopy(other)
        elif how == 'right':
            obj = copy.deepcopy(other)
            mask1 = obj.to_binary()
            mask2 = copy.deepcopy(self)
        else:
            pass

        intersection = cv2.bitwise_and(mask1.matrix, mask2.matrix)

        mask2.matrix[intersection] = 0

        obj.matrix += mask2.matrix

        sys.stdout.write('done\n')
        return obj
    
    def save(self, prefix='out'):
        
        np.savetxt(f'{prefix}.mask.txt', self.matrix, fmt='%d')

        return 

    def overlayoutlines(self, image=None, prefix=None):

        sys.stdout.write('draw outlines ... ')
        sys.stdout.flush()
        
        import skimage.io
        import skimage.segmentation

        if isinstance(image, str):
            image = skimage.io.imread(image)

        outlines = skimage.segmentation.mark_boundaries(
                image, 
                self.matrix, 
                color=(1, 0, 0),
                mode='inner',
                )
        b, g, r = cv2.split(outlines)

        sys.stdout.write('{} labels\n'.format(len(np.unique(self.matrix))))

        mask = np.isin(b, [1])
        image[mask] = 255
        
        if prefix:
            np.savetxt(f'{prefix}.outlines.txt', b, fmt='%d')
            cv2.imwrite(f'{prefix}.outlines.png', image)
        return b, image

def thres_mask(image, out_prefix=None):
    image = cv2.imread(image, 0)

    _, th = cv2.threshold(image, 20, 255, cv2.THRESH_BINARY)

    if out_prefix:
        cv2.imwrite(f'{prefix}.mask.tif', th)
    return th

def mixture_seg(cell_mask, tissue_mask, blur_mask, image=None, prefix='out',):

    cell_mask = Stoarr(cell_mask)
    tissue_mask = Stoarr(tissue_mask)
    blur_mask = Stoarr(blur_mask)

    blur_mask = blur_mask.minimum_filter(
            footprint='octagon',
            ksize=(7, 4)
    )

    orig_cell_mask = cell_mask.intersection(
            tissue_mask, 
            label_area_cutoff=0.3
    )

    cell_mask = orig_cell_mask.filter_by_matrix(
            on=blur_mask, 
            max_value=90, 
            draw=True,
            prefix=prefix
    )

    cell_mask = cell_mask.filter_by_diameter(
            min_size=3,
            max_size=None,
    )

    tissue_mask = orig_cell_mask.subtract(cell_mask)

    bin_mask = tissue_mask.binning(
            bin_size=20
    )
    
    mix_mask = cell_mask.merge(
            bin_mask, 
            how='left'
    )

    mix_mask.save(prefix=prefix)

    outlines, image = mix_mask.overlayoutlines(
            image=image, 
            prefix=prefix
    )

    return outlines, image


class Options:
    sargs = 'hc:t:b:i:p:'
    largs = ['help', 'cell=', 'tissue=', 'blur=', 'image=', 'prefix=']
    def __init__(self, cmdline):
        self.cmdline = cmdline
        try:
            self.options, self.not_options = getopt.gnu_getopt(self.cmdline, 
                    self.sargs, self.largs)
        except getopt.GetoptError as err:
            sys.stderr.write('*** Error: {}\n'.format(err))
            sys.exit()
        if (not self.options) or self.not_options:
            sys.stderr.write('*** Error: {}\n'.format(err))
            sys.stdout.write(self.usage)
            sys.exit()
        for opt, arg in self.options:
            if opt in ['-h', '--help']:
                sys.stdout.write(self.usage)
                sys.exit()
            elif opt in ['-c', '--cell']:
                self.cell_mask = arg
            elif opt in ['-t', '--tissue']:
                self.tissue_mask = arg
            elif opt in ['-b', '--blur']:
                self.blur_mask = arg
            elif opt in ['-i', '--image']:
                self.image = arg
            elif opt in ['-p', '--prefix']:
                self.prefix = arg
    @property
    def usage(self):
        return __doc__.format(__file__)

if __name__ == '__main__':

    opts = Options(sys.argv[1:])

    mixture_seg(opts.cell_mask, 
                opts.tissue_mask, 
                opts.blur_mask, 
                image=opts.image, 
                prefix=opts.prefix
                )
