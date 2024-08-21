.. _`downstream-analysis`:

.. sidebar:: quick index

	.. rubric:: SPC clustering
	| `SPC <https://github.com/lskfs/SPC>`_: Spatial Proximity based Clustering ☀ 
	| :ref:`/tutorials/spc.rst`

	.. rubric:: Spatial gradient gene profiling
	| ::ref:`/tutorials/gene_profiling.rst` 

	.. rubric:: Spatial pattern clustering
	| ::ref:`/tutorials/clustering.rst` 

	.. rubric:: Potential PCGs detection
	| ::ref:`/tutorials/mining.rst` 

	.. rubric:: Principal component analysis of PCGs
	| ::ref:`/tutorials/pca.rst` 

	.. rubric:: Cell-cell connection analysis
	| ::ref:`/tutorials/CCC.rst` 

	.. rubric:: Cell colocation estimation
	| ::ref:`/tutorials/CellCellColocation3D.rst` 

	.. rubric:: Micro-environment estimation
	| ::ref:`/tutorials/calculate-niches-on-coords.ipynb` 

======================
Downstream Analysis
======================
This section contains description of downstream analysis including SPC clustering, PCG analysis, cell-cell connection analysis, cell colocation estimation and meta-niche analysis.

SPC clustering
--------------
With the cell-level clustering results, non-continous segmentation on the 3D spatial coordinates was performed to aggratate cells with functional and spatial proximity following :ref:`/tutorials/SPC.rst`.

Spatial gradient gene profiling
-------------------------------
Each 3D individual was digitally split into body fragments with equal length along A/P, M/L, and D/V axes, from which the gene expression of each fragment was calculated by averaging sctransform-based data following :ref:`/tutorials/gene_profiling.rst`.

.. note:: 
    The digitally split was performed after the body straightening following :ref:`/tutorials/body_straightening.rst`.

Spatial pattern clustering
--------------------------
Hierarchical density-based clustering algorithm was applied to divide genes into groups sharing similar expression and positional profiles along the body axis. The unclustered genes were then mapped to the specific group with the largest possibility given by linear regression. Detailed process was following :ref:`/tutorials/clustering.rst`. 

Potential PCGs detection
------------------------
Spatial gradient genes with high confidence were selected as potential PCGs following :ref:`/tutorials/mining.rst`. 

Principal component analysis of PCGs
------------------------------------
Principal component analysis was applied to the binned expressions of PCGs along three body axes following :ref:`/tutorials/pca.rst`.

Cell-cell connection analysis
-----------------------------
Spatial cell-cell interaction was estimated based on ligand-receptor pairs from WNT, BMP, COLLAGEN and NOTCH pathways in each 3D planarian individual in the context of each SPC cell to its neighbors within five SPC cell layers. Spatial connection network of SPC cells based on their relationships of cell-cell interactions following :ref:`/tutorials/CCC.rst`.

Cell colocation estimation
--------------------------
Global spatial distribution similarity was estimated by KL in 3D coordinates following :ref:`/tutorials/CellCellColocation3D.rst`.

Micro-environment estimation
----------------------------
The composition of neighbor cells for neoblast in the 3D space was used to evaluate the microenvironment of neoblast in homeostasis and regeneration following :ref:`/tutorials/calculate-niches-on-coords.ipynb`.

