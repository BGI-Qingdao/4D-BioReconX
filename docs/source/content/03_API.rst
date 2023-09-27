API
===

.. toctree::
   :maxdepth: 1

   031_API_GEM3D_toolkit

The `GEM3D_toolkit` module includes input and output parts. 

The input part supports reading files in multiple formats to generate StereoExpData object, \
such as GEM, GEF, H5ad (from Scanpy and Seurat) and so on. Stereo-seq sequencing data were preprocessed using \
SAW to generate spatial gene expression matrices in GEM format. Image files, which are usually in TIFF format, \
generated during the stereo-seq library construction process should be ready for further process together with the GEM files.


.. toctree::
   :maxdepth: 1

   032_API_cell-segmentation

  A tool for single-cell segmentation for ST data. We started by segmenting the nuclei by ssDNA staining image \
  and later transferred the cell boundaries to the aligned gene expression image. \
  For the more complex tasks of nuclei segmentation due to crowded cell environments or out-of-focus regions, \
  we used the Laplacian algorithm to calculate the blur value for each pixel based on sliding kernels within 3 neighbor pixels surrounded. \
  Pixels with blur value less than 30 were labeled as regions that failed to be identified and retained bin15 \
  (10 µm in width and height) partition instead of boundary determinization. \
  With the precise 2D registration of the gene expression image to the pairwise ssDNA-staining image, \
  the boundaries detected for each cell or bin were mapped onto the gene expression spatial color-code image. \
The DNB level gene expression matrices were further aggregated into putative cells ready for downstream analysis.


.. toctree::
   :maxdepth: 1
