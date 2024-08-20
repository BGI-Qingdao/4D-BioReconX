.. _`data-preprocess`:
======================
Data preprocessing
======================
This section contains description of preprocessing of Stereo-seq data.

Stereo-seq Analysis Workflow
----------------------------
Stereo-seq sequencing data were preprocessed using `SAW <https://github.com/STOmics/SAW>`_ to generation spatial gene expression matrices in `**GEM** format <https://stereopy.readthedocs.io/en/latest/Tutorials/IO.html#GEM>`_.

Image files, which are usually in **TIFF** format, generated during the stereo-seq library construction process should be ready for further process together with the GEM files.

.. note::
	Data of each section were packed with a sinlge **Tissue Section** ID in STOmicsDB database, including staining image *(.tif)*, GEM file *(.gem)* as well as relevant annotation information *(.txt)*.
		* `[STTS0000461 - 515] <https://db.cngb.org/stomics/project/STT0000028>`_

GEM file and TIFF file for each chip were manually cropped according to the ROI region to extract the expression and image data of each individual in the tissue section. 

Image registration
------------------
Each pair of GEM file and TIFF file was then went through image registration following :ref:`/tutorials/mirror.rst`. 

#. input:
   * gem file
   * image file
#. output:
   * affine transformation matrix

Cell boundary detection
-----------------------
Each ssDNA image was went through object detection and cell segmentation following :ref:`/tutorials/cell-segmentation.rst`. 

#. input:
   * image file
#. output:
   * cell mask matrix

DNB aggregation
---------------
With cell mask and affine transformation matrix, the original gem file would be aggregated into putative cells following `cell aggregation <https://spacipy.readthedocs.io/en/latest/intro/gem_process.html#>`_.

#. input:
   * affine transformation matrix
   * cell mask matrix
   * gem file
#. output:
   * cell-gene expression matrix in h5d format

Cell clustering
---------------
The generated Cell-by-gene matrix of each individual was processed following the `Seurat integration <https://satijalab.org/seurat/archive/v4.3/integration_introduction>`_ workflow, from which the clustering results were manually annotated to obtain lineage information of each cluster.

#. input:
   * cell-gene matrix
#. output:
   * cell-lineage annotation

3D alignment
------------
Serial annotation images were aligned by the similarities of morphology and annotation color code pair by pair folling :ref:`/tutorials/seam.rst`.

#. input:
   * cell-lineage annotation
   * images
#. output:
   * aligned 3D coordinates

3D mesh building
----------------
With the aligned 3D coordinates and cell annotation, tissue level triangular meshes were constructed folling :ref:`/tutorials/meshgen.rst`.

#. input:
   * 3D coordinates
#. output:
   * triangular meshes

Data availability
---------------------------------
All data generated in this study were deposited at CNGB Nucleotide Sequence Archive (accession code: STT0000028).

Processed data (h5ad format) can be interactively explored and doanload from our `PRISTA4D <https://db.cngb.org/stomics/prista4d>`_ database. 
