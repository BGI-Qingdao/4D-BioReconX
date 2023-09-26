#!/usr/bin/env python3
"""usage: {} --image test.tif --gem test.gem --outdir ./ --prefix out [--threads 10] [--skip]

options:
    -i, --image     <STR>       tiff image filename
    -g, --gem       <STR>       gem filename
    -o, --outdir    <STR>       directory saving the output, default ./
    -p, --prefix    <STR>       output prefix, default 'out'
    -t, --threads   <INT>       cpu numbers, default 10
    --skip          <BOOL>      skip the filtering of neibor
    --source        <STR>       image source, one of 'leica' and 'motic', default 'leica'
                                PLS ask manager for help if you are not sure
    --range         <STR>       range of cell size in um, default '8,10'
    --cpi           <STR>       custom cellprofile pipeline file
    --name          <STR>       image name in cpi file
    --object        <STR>       object name which will saved as final detected object
"""
import os
os.environ['OPENCV_IO_MAX_IMAGE_PIXELS'] = pow(2, 40).__str__()
import sys
import getopt
from pathlib import Path
import pandas as pd
import numpy as np
from itertools import repeat
#from multiprocessing import set_start_method
#set_start_method('spawn')
#import multiprocessing
#from multiprocessing import get_context
import cv2
import skimage
import scipy.sparse

import bioformats.formatreader
import cellprofiler_core.pipeline
import cellprofiler_core.preferences
import cellprofiler_core.utilities.zmq
import cellprofiler_core.utilities.java
from cellprofiler_core.setting.subscriber import LabelSubscriber
from cellprofiler_core.setting.range import IntegerRange

def cp_main(image, name='DNA', cpi=None, skip=False, psize=None, saved_object='IdentifySecondaryObjects'):

    # load pipeline from cpi file
    print('load pipeline from {}'.format(cpi))
    pipeline = cellprofiler_core.pipeline.Pipeline()
    pipeline.load(cpi)

    # get modules list
    modules = pipeline.modules()

    # setup image_set
    image_set = cellprofiler_core.image.ImageSet(0, {'name':name}, name)
    #x = skimage.io.imread(image, as_gray=True)
    x = cv2.imread(str(image), 0)
    x[x > 230] = 230
    image_x = cellprofiler_core.image.Image(x, path_name=image.parent, file_name=image.name)
    image_set.add(name, image_x)
    #skimage.io.imsave('ssDNA.png', image_set.get_image('ssDNA').pixel_data,)

    # init empty object_set
    object_set = cellprofiler_core.object.ObjectSet()
    #for name in ['Nuclei', 'FilterNuclei', 'Cell', 'Cytoplasm']:
    #    objects = cellprofiler_core.object.Objects()
    #    object_set.add_objects(objects, name)

    measurements = cellprofiler_core.measurement.Measurements()

    workspace = cellprofiler_core.workspace.Workspace(
            pipeline, 
            modules,
            image_set,
            object_set,
            measurements,
            [image_set]
            )

    for module in modules:
        if skip:
            if module.module_name in ['MeasureObjectNeighbors', 
                    'FilterObjects', 'MeasureObjectSizeShape']:
                continue
            elif module.module_name == 'IdentifySecondaryObjects':
                module.x_name.value = 'Nuclei'
                module.size_range = IntegerRange('size range for cell', value=tuple(psize))
            elif module.module_name == 'IdentifyTertiaryObjects':
                module.primary_objects_name = LabelSubscriber(
                        'Select the smaller identified objects', 
                        'Nuclei')
        print('\t' + module.module_name, flush=True)
        module.run(workspace)

    objects = workspace.object_set.get_objects(saved_object)
    try:
        celloutlines = workspace.image_set.get_image('CellOutlines')
    except:
        sys.stderr.write('cell outlines not get\n')
        celloutlines = None

    return objects, celloutlines

class Options:
    sargs = 'hi:n:o:p:t:m:s:'
    largs = ['help', 'image=', 'name=', 'outdir=', 'prefix=', 'suffix=', 'skip', 
            'threads=', 'cpi=', 'object=', 'mask=', 'source=']
    cpi = str(Path(__file__).parent.resolve() / 'default.cppipe')
    
    # how many pixels for per um
    pixel_size = {
        'motic': 2,
        'leica': 2.0819
    }

    def __init__(self, cmdline):
        self.cmdline = cmdline
        self.image = None
        self.image_name = 'DNA'
        self.outdir = Path('./').resolve()
        self.prefix = 'out'
        self.suffix = 'txt'
        self.skip = False
        self.threads = 10
        self.source = 'leica'
        self.cell_size = [8, 10]
        self.saved_object = 'IdentifySecondaryObjects'
        self.mask = None
        try:
            self.options, self.not_options = getopt.gnu_getopt(self.cmdline, 
                    self.sargs, self.largs)
        except getopt.GetoptError as err:
            sys.stderr.write('***Error: {}\n'.format(err))
            sys.exit()
        if not self.options:
            sys.stderr.write(self.usage)
            sys.exit()
        for opt, arg in self.options:
            if opt in ['-h', '--help']:
                sys.stderr.write(self.usage)
                sys.exit()
            elif opt in ['-i', '--image']:
                self.image = Path(arg).resolve()
            elif opt in ['-m', '--mask']:
                self.mask = arg
            elif opt in ['-n', '--name']:
                self.image_name = arg
            elif opt in ['-o', '--outdir']:
                self.outdir = Path(arg).resolve()
            elif opt in ['-p', '--prefix']:
                self.prefix = arg
            elif opt in ['-s', '--suffix']:
                self.suffix = arg
            elif opt in ['--skip']:
                self.skip = True
            elif opt in ['--threads']:
                self.threads = int(arg)
            elif opt in ['--source']:
                assert arg in self.pixel_size
                self.source = arg
            elif opt in ['--range']:
                self.cell_size = [int(x) for x in arg.split(',')]
            elif opt in ['--cpi']:
                self.cpi = arg
            elif opt in ['--object']:
                self.saved_object = arg
        
        pixels_per_um = self.pixel_size[self.source]
        self.pixel_range = [x * pixels_per_um for x in self.cell_size]
    @property
    def usage(self):
        return __doc__.format(__file__)

def main():
    opts = Options(sys.argv[1:])
    if not opts.outdir.exists():
        opts.outdir.mkdir(parents=True, exist_ok=True)

    if opts.image and not opts.mask:
        try:
            cellprofiler_core.preferences.set_headless()
            cellprofiler_core.preferences.set_temporary_directory(opts.outdir)
            cellprofiler_core.preferences.set_default_output_directory(opts.outdir)

            cellprofiler_core.utilities.java.start_java()

            print('Starting cellprofiler identify ...', flush=True)
            objects, celloutlines = cp_main(
                    opts.image, name=opts.image_name, 
                    cpi=opts.cpi, skip=opts.skip, psize=opts.pixel_range, 
                    saved_object=opts.saved_object
                    )
            print('Cell objects and outlines generated', flush=True)
        except Exception as err:
            sys.stderr.write('***Error: {}\n'.format(err))
        finally:
            cellprofiler_core.utilities.zmq.join_to_the_boundary()
            bioformats.formatreader.clear_image_reader_cache()
            cellprofiler_core.utilities.java.stop_java()

        print('Saving labled cells ...', flush=True)
        mask_file = opts.outdir / '{}_mask.txt'.format(opts.prefix)
        np.savetxt(mask_file, objects.segmented, fmt='%d')
        mask = objects.segmented
    elif opts.image and opts.mask:
        mask = np.loadtxt(opts.mask)
        celloutlines = None

    celloutlines_file = str(opts.outdir / f'{opts.prefix}_CellOutlines.{opts.suffix}')
    if celloutlines is not None:
        print('Saving cell outlines ...', flush=True)
        b, g, r = cv2.split(celloutlines.pixel_data)
        np.savetxt(celloutlines_file, b, fmt='%d')
    else:
        print('Detecting cell outlines ...', flush=True)
        import skimage.segmentation
        image = cv2.imread(str(opts.image))
        outlines = skimage.segmentation.mark_boundaries(
                image,
                mask,
                color=(1, 0, 0),
                mode='inner',
                )
        b, g, r = cv2.split(outlines)
        if opts.suffix == 'txt':
            np.savetxt(celloutlines_file, b, fmt='%d')
        elif opts.suffix == 'npy':
            np.save(celloutlines_file, b)
    print('Saving masked image ...', flush=True)
    b = np.isin(b, [1])
    image = cv2.imread(str(opts.image))
    image[b] = (255, 0, 0)
    celloutlines_png = str(opts.outdir / '{}_CellOutlines.png'.format(opts.prefix))
    cv2.imwrite(celloutlines_png, image)
    
    return 

if __name__ == '__main__':
    main()

