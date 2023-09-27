
import os
os.environ['OPENCV_IO_MAX_IMAGE_PIXELS'] = pow(2, 40).__str__()
import sys
from pathlib import Path
import numpy as np
import cv2

import bioformats.formatreader
import cellprofiler_core.pipeline
import cellprofiler_core.preferences
import cellprofiler_core.utilities.zmq
import cellprofiler_core.utilities.java
from cellprofiler_core.setting.subscriber import LabelSubscriber
from cellprofiler_core.setting.range import IntegerRange

def _pycellprofilter(image, name='DNA', cpi=None, saved_object='IdentifySecondaryObjects'):

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

    # init workspace
    object_set = cellprofiler_core.object.ObjectSet()

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
        sys.stdout.write(f'... {module.module_name}\n')
        module.run(workspace)

    objects = workspace.object_set.get_objects(saved_object)
    try:
        celloutlines = workspace.image_set.get_image('CellOutlines')
    except:
        sys.stderr.write('cell outlines not get\n')
        celloutlines = None

    return objects, celloutlines

def pycellprofiler(image, mask, save_prefix='test', img_boundary=True,
        cpi='./default.cppipe', 
        saved_object='IdentifySecondaryObjects',
        outdir='./outdir', tmpdir='./tmpdir', ):

    outdir, tmpdir = Path(outdir), Path(tmpdir)
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)

    try:
        cellprofiler_core.preferences.set_headless()
        cellprofiler_core.preferences.set_temporary_directory(outdir)
        cellprofiler_core.preferences.set_default_output_directory(outdir)

        cellprofiler_core.utilities.java.start_java()

        sys.stdout.write('Starting cellprofiler identify ...\n')
        objects, celloutlines = _pycellprofilter(
                image, 
                name=image_name, 
                cpi=cpi, 
                saved_object=saved_object
                )
        sys.stdout.write('Cell objects and outlines generated\n')
    except Exception as err:
        sys.stderr.write('***Error: {}\n'.format(err))
    finally:
        cellprofiler_core.utilities.zmq.join_to_the_boundary()
        bioformats.formatreader.clear_image_reader_cache()
        cellprofiler_core.utilities.java.stop_java()

    sys.stdout.write('Saving labled cells ...\n')

    mask = objects.segmented
    b, g, r = cv2.split(celloutlines.pixel_data)
    if save_prefix:
        mask_file = str(outdir / f'{prefix}_mask.txt')
        np.savetxt(mask_file, mask, fmt='%d')

        boundary_file = str(outdir / f'{prefix}_boundary.txt')
        np.savetxt(boundary_file, b, fmt='%d')
    
    if img_boundary:
        image = img_outliner(image, boundary=b, 
                    save=f'{prefix}.celloutlines.png'
                    )
    
    return mask, b

def boundary_detect(mask, image, save_prefix='cell'):
    import skimage.segmentation

    image = cv2.imread(str(image))
    outlines = skimage.segmentation.mark_boundaries(
            image,
            mask,
            color=(1, 0, 0),
            mode='inner',
            )
    
    b, g, r = cv2.split(outlines)

    if save:
        np.savetxt(f'{prefix}.boundary.txt', b, fmt='%d')
    
    image = img_outliner(image, boundary=b, 
                         save=f'{prefix}.celloutlines.png'
                         )
    
    return b

def img_outliner(image, boundary, save='celloutlines.png'):
    if isinstance(image, str):
        image = cv2.imread(image)
    
    mask = np.isin(boundary, [1])
    image[mask] = (255, 0, 0)
    
    if save:
        cv2.imwrite(save, image)
    return image
