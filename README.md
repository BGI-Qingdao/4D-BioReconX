[![stars](https://img.shields.io/github/stars/BGI-Qingdao/4D-BioReconX?logo=GitHub&color=yellow)](https://img.shields.io/github/stars/BGI-Qingdao/4D-BioReconX) 
[![docs](https://img.shields.io/static/v1?label=docs&message=4d-bioreconx&color=green)](https://4D-bioreconx.readthedocs.io/en/latest/index.html)

## **4D-BioReconX**: a bioinformatic framework for reconstructing 4D spatial transcriptomics atlas and spatiotemporal analyses+

[Installation](https://4d-bioreconx.readthedocs.io/en/latest/content/00_Installation.html) - 
[Quick start](https://4d-bioreconx.readthedocs.io/en/latest/Tutorials/index.html) - 
[Documentation](https://4d-bioreconx.readthedocs.io/en/latest/index.html)

<img src="https://github.com/BGI-Qingdao/4D-BioReconX/blob/main/docs/source/_static/4D-BioReconX_workflow_v1.0.0.jpg" width="100%" height="100%">


[comment]: <> (![4D-BioReconX]&#40;./docs/source/_static/4D-BioReconX.PNG&#41;)


**4D-BioReconX** is a scalable and versatile framework designed to construct a comprehensive four-dimensional (4D) spatial transcriptomics atlas of whole animals at single-cell resolution. This approach enables us to track the intricate spatiotemporal dependencies of morphogenetic gradients and regenerative patterning. It runs based on spatial transcriptomics data, such as Stereo-seq (spatial enhanced resolution omics sequencing) data. Notably, we are still working on the improvement of performance and calculation efficiency.

## Discussion 
Please use GitHub [issues](https://github.com/BGI-Qingdao/4D-BioReconX/issues) tracker for reports and discussions of:
 - Bug reports
 - Document and data issues
 - Feature requirements
 - ...

## Contribution 
**4D-BioReconX** is in a stage of rapid development so that we will carefully consider all aspects of your proposal. We hope future input will be given by both users and developers.


## File tree
**To enhance clarity and prevent any potential confusion, we have organized this GitHub repository into a structured file tree, complete with detailed annotation for each code segment, tailored to their respective purposes.

.
├── LICENSE # MIT License
├── OtherModules    # Independent analysis modules for specifc requirements. Please see the README in each folder.
│   ├── README.md   # General README for OtherModules.
│   ├── blastemadetection   #  Blastema region detection.
│   │   ├── BlastemaByWound_v2.py   # Python script for blastemadetection.
│   │   └── README.md   # README for blastemadetection.
│   └── bodystrengthening   # Body straightening.
│       ├── README.md   # README for bodystrengthening.
│       ├── adjust_APDV_example.ipynb   # Jupyter Notebook for bodystrengthening along AP-DV plane.
│       └── adjust_APML_example.ipynb   # Jupyter Notebook for bodystrengthening along AP-ML plane.
├── PolarityAnalysis    # To investigate spatial patterns of gene expression capture positional gradients, we conducted a qualitative analysis of positional control genes (PCGs) spanning the planarian body axis.
│   ├── README.md   # General README for PolarityAnalysis.
│   ├── newpcgmining    # Find genes that could be PCGs by calculating and evaluating correlation coefficient values between known PCGs (from published studies) and target genes.
│   │   ├── filter_gene.py  # Python script for filtering genes.
│   │   ├── kl_simple.py    # Python script for searching new PCGs by KL divergence.
│   │   ├── pearson_simple.py    # Python script for searching new PCGs by Pearson correlation coefficient.
│   │   └── spearman_simple.py   # Python script for searching new PCGs by Spearman rank correlation coefficient.
│   ├── patternclustering   # perform clustering via hdbscan and LogicRegression.
│   │   ├── cluster_change.py   # Python script for detecting cluster change with time.
│   │   ├── clustering.py   # Python script for hdbscan clustering.
│   │   ├── draw.heatmap.py # Python script for drawing heatmap.
│   │   ├── draw.line.py    # Python script for drawing line.
│   │   ├── gen_data.sh # Shell script for preparing data.
│   │   ├── pcg_pattern.py  # Python script for computing pattern.
│   │   ├── plot_heatmap.py # Python script for plotting heatmap.
│   │   ├── run_pcgs_oneindv.sh # Shell script for running task individually.
│   │   └── sort_hcluster.py    # Python script for sorting clusters.
│   ├── pca # PCG data in PCA space
│   │   ├── pca.py  # Python script for pca.
│   │   └── plot_pca.py # Python script for plotting trace in PCA space.
│   ├── pcgprofiling    # Gene spatial profiling.
│   │   ├── AP_Profiling_HVG_log_mean.py    # Python script for spatial profiling of HVGs along A/P axis.
│   │   └── README.md   # README for pcgprofiling.
│   └── util.py # Python script for utility.
├── README.md   # General README for 4D-BioReconX.
├── SPCAnalysis # On the basic assumption that a functional unit must be comprised of a group of cells rather than a single cell, and in the case of neoblast cells in planarian dispersedly distributed in space, we used a functional-spatial-proximity-based strategy to define hubs for cells in similar states or with similar functions.
│   ├── README.md   # README for SPCAnalysis.
│   ├── __init__.py # Python script for initialization.
│   ├── cell-cellcolocalization # Calculate cell-cell colocation array by KL in 3D coordinate and then generate its MST.
│   │   ├── CellColocation.py   # Python script for cell-cellcolocalization.
│   │   └── README.md   # README for cell-cellcolocalization.
│   ├── digitallineagetracing   # Lineage similarity calculation.
│   │   ├── README.md   # README for digitallineagetracing.
│   │   ├── s01.CellCor.py  # Python script for calulating the simularity array for each pair of cells between thoes two datasets and saving the closest cell in another dataset for every cells.
│   │   ├── s02.get_number_anno.py  # Python script for calculating the mapping number of cells for all clusters.
│   │   ├── s03.filter.py   # Python script for filtering weak connections.
│   │   └── s04.Sankey_plot.ipynb   # Jupyter Notebook for plotting the final sankey plot.
│   ├── meta-nicheanalysis # Calculating surrounding cell component.
│   │   ├── calculate-niches-on-coords.ipynb    # Jupyter Notebook for meta-nicheanalysis.
│   │   ├── demo.txt    # Demo data.
│   │   └── niche.py    # Python script for calculating niches.
│   ├── spc-basedccc    # SPC-based cell-cell communication analysis.
│   │   ├── CellChatDB.human.interaction.planarian.txt  # Modified ligand-receptor database from CellChatDB.
│   │   ├── NICHES_3D.R # R package for NICHES analysis.
│   │   └── single_cell_cross_talk.R    # R package for CellChat analysis.
│   └── spcclustering   # Spatial Proximity based Clustering.
│       ├── README.md   # README for spcclustering.
│       └── cell2spc.py # Python script for spcclustering.
├── WACCA   #
│   ├── 3dmeshreconstruction    # Mesh generating pipeline.
│   │   ├── README.md   # README for 3dmeshreconstruction.
│   │   ├── gen_stacked_tif.py  # Python script for generating stacked slices.
│   │   └── reset_obj.py    # Python script for resetting object.
│   ├── README.md   # General README for WACCA
│   ├── cellsegmentation    # Cell segmentation based on ssDNA.
│   │   ├── __pycache__ # pycache files.
│   │   │   ├── objseg.cpython-311.pyc
│   │   │   ├── objseg.cpython-312.pyc
│   │   │   └── objseg.cpython-38.pyc
│   │   ├── blur_detect.py  # Python script for blurred region detection.
│   │   ├── cell-segmente-on-ssDNA.ipynb    # Jupyter Notebook for cellsegmentation.
│   │   ├── cellsegment.py  # Python script for cellsegmentation.
│   │   ├── clahe_enhance.py    # Python script for enhancing image contrast by CLAHE.
│   │   ├── default.cppipe  # CellProfiler script.
│   │   ├── demo.tif    # Demo figure.
│   │   ├── objseg.py   # Python script for image processing such as coloring.
│   │   └── segarr.py   # Python script for image masking.
│   ├── datapreprocessing   # Data preprocessing.
│   │   ├── GEM3D_toolkit.py    # Python script for main entries.
│   │   ├── README.md   # README for datapreprocessing.
│   │   ├── affine_gem.py   # Python script for affine transformation for GEM files.
│   │   ├── affine_h5ad.py  # Python script for affine transformation for H5AD files.
│   │   ├── affine_ssdna.py # Python script for affine transformation for ssDNA image files.
│   │   ├── affine_txt.py   # Python script for affine transformation for txt files.
│   │   ├── apply_cells.py  # Python script for applying transformation to cells.
│   │   ├── apply_registration.py   # Python script for applying registration.
│   │   ├── chop_gem.py # Python script for chopping GEM files. 
│   │   ├── chop_image.py   # Python script for chopping images.
│   │   ├── chop_paste.py   # Python script for chopping and pasting.
│   │   ├── draw_heatmap.py # Python script for drawing heatmap.
│   │   ├── gem_to_h5ad.py  # Python script for converting GEM to H5AD.
│   │   ├── gem_xy.py   # Python script for GEM files.
│   │   ├── h5ad_dataframe.py   # Python script for H5AD files.
│   │   ├── image_blend.py  # Python script for blending images.
│   │   ├── mask_gem.py # Python script for masking GEM files.
│   │   ├── mask_h5ad.py    # Python script for masking H5AD files.
│   │   ├── merge_h5ad.py   # Python script for merging H5AD files.
│   │   ├── save_miscdf.py  # Python script for saving checking.
│   │   ├── slice_dataframe.py  # Python script for converting dataframe.
│   │   ├── split_gem.py    # Python script for spliting GEM files.
│   │   └── trakEM2_to_affine.py    # Python script for alignment using trakEM2.
│   ├── mirrorregistration  # the mirro pipeline aim to registrate ssDNA image to GEM coordinate system with minimum registration error.
│   │   ├── MIRROR.py   # Python script for mirrorregistration.
│   │   ├── README.md   # README for mirrorregistration.
│   │   ├── gem_to_gemc.py  # Python script for converting GEM to GEMC.
│   │   ├── gemc_to_h5ad.py # Python script for converting GEMC to H5AD.
│   │   ├── prepare_registration_heatmap.py # Python script for preparing gene heatmap for registration.
│   │   ├── prepare_registration_ssdna.py   # Python script for preparing ssDNA images for registration.
│   │   ├── second_registration.py  # Python script for second-round registration.
│   │   └── slice_dataframe.py  # Python script for converting dataframe.
│   └── seamalignment   # seam pipeline aims to align and merge serial slices (each slice in an H5AD file) into one 3D atlas (one H5AD file with aligned 3D.
│       ├── README.md   # README for seamalignment.
│       ├── SEAM.py # Python script for seamalignment.
│       ├── apply_alignment.py  # Python script for applying alignment to all coordinates.
│       ├── get_xml_matrix.py   # Python script for obtaining transformation matrix.
│       └── prepare_alignment_image.py  # Python script for preparing images for 3D alignment.
├── docs    # Files for readthedocs.
│   ├── Makefile
│   ├── make.bat
│   ├── requirements.txt
│   └── source
│       ├── Tutorials
│       │   ├── BIO.rst
│       │   ├── CCC.rst
│       │   ├── CellCellColocation3D.rst
│       │   ├── Meta-Niche.rst
│       │   ├── Other.rst
│       │   ├── PCG.rst
│       │   ├── Preprocess.rst
│       │   ├── SPC.rst
│       │   ├── adjust_APDV_example.ipynb
│       │   ├── adjust_APML_example.ipynb
│       │   ├── assign_blastema_region.rst
│       │   ├── body_straightening.rst
│       │   ├── calculate-niches-on-coords.ipynb
│       │   ├── cell-segmentation.rst
│       │   ├── cell-segmente-on-ssDNA.ipynb
│       │   ├── clustering.rst
│       │   ├── data-preprocess.rst
│       │   ├── digital-lineage-s1-s3.rst
│       │   ├── digital-lineage-s4.rst
│       │   ├── digital-lineage.rst
│       │   ├── gene_profiling.rst
│       │   ├── index.rst
│       │   ├── meshgen.rst
│       │   ├── mining.rst
│       │   ├── mirror.rst
│       │   ├── pca.rst
│       │   ├── s04.Sankey_plot.ipynb
│       │   └── seam.rst
│       ├── _static
│       │   ├── 4D-BioReconX_workflow_v1.0.0.jpg
│       │   ├── assign_blastema_region_2d.png
│       │   ├── assign_blastema_region_aligned.png
│       │   ├── assign_blastema_region_final.png
│       │   ├── assign_blastema_region_raw.png
│       │   ├── assign_blastema_region_workflow.png
│       │   ├── fig.png
│       │   ├── gene_profiling_workflow.png
│       │   ├── meshgen_workflow.png
│       │   ├── mirror_ipo.png
│       │   ├── mirror_workflow.png
│       │   ├── sankey-01.png
│       │   ├── sankey-02.png
│       │   ├── seam_ipo.png
│       │   ├── seam_workflow.png
│       │   ├── spc_annotation.png
│       │   ├── spc_compare.jpg
│       │   ├── spc_leiden.png
│       │   └── spc_leiden.spc.png
│       ├── conf.py
│       ├── content
│       │   ├── 00_Installation.rst
│       │   ├── 01_Basic_Usage.rst
│       │   └── 03_References.rst
│       └── index.rst
├── planarian manuscript    # Detailed usage of this framework in the planarian manuscript.
│   └── mappingtoplanarian.jpg  # Mapping relations.
├── pyproject.toml  # Readthedocs building.
└── requirements.txt    # Installation requirements.



