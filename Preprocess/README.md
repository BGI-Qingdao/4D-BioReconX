# Data preprocess and 3D atlas construction

## Basic data preprocess

see ```data-preprocess``` folder

## The WACCA pipeline

* step01 - cell_segmentation: foreach slice, generate cell segmentation via ssDNA image.
* step02 - mirror: foreach slice, registrate ssDNA to the GEM coordinate system.
* step03 - seam: align multiple slices to one unified 3D coordinate system.
* step04 - meshgen: generate surface model by the stacked, annotated images. 

