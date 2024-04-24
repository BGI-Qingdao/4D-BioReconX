.
├── LICENSE
├── OtherModules
│   ├── README.md
│   ├── blastemadetection
│   │   ├── BlastemaByWound_v2.py
│   │   └── README.md
│   └── bodystrengthening
│       ├── README.md
│       ├── adjust_APDV_example.ipynb
│       └── adjust_APML_example.ipynb
├── PolarityAnalysis
│   ├── README.md
│   ├── newpcgmining
│   │   ├── filter_gene.py
│   │   ├── kl_simple.py
│   │   ├── pearson_simple.py
│   │   └── spearman_simple.py
│   ├── patternclustering
│   │   ├── cluster_change.py
│   │   ├── clustering.py
│   │   ├── draw.heatmap.py
│   │   ├── draw.line.py
│   │   ├── gen_data.sh
│   │   ├── pcg_pattern.py
│   │   ├── plot_heatmap.py
│   │   ├── run_pcgs_oneindv.sh
│   │   └── sort_hcluster.py
│   ├── pca
│   │   ├── pca.py
│   │   └── plot_pca.py
│   ├── pcgprofiling
│   │   ├── AP_Profiling_HVG_log_mean.py
│   │   └── README.md
│   └── util.py
├── README.md
├── SPCAnalysis
│   ├── README.md
│   ├── __init__.py
│   ├── cell-cellcolocalization
│   │   ├── CellColocation.py
│   │   └── README.md
│   ├── digitallineagetracing
│   │   ├── README.md
│   │   ├── s01.CellCor.py
│   │   ├── s02.get_number_anno.py
│   │   ├── s03.filter.py
│   │   └── s04.Sankey_plot.ipynb
│   ├── meta-nicheanalysis
│   │   ├── calculate-niches-on-coords.ipynb
│   │   ├── demo.txt
│   │   └── niche.py
│   ├── spc-basedccc
│   │   ├── CellChatDB.human.interaction.planarian.txt
│   │   ├── NICHES_3D.R
│   │   └── single_cell_cross_talk.R
│   └── spcclustering
│       ├── README.md
│       └── cell2spc.py
├── WACCA
│   ├── 3dmeshreconstruction
│   │   ├── README.md
│   │   ├── gen_stacked_tif.py
│   │   └── reset_obj.py
│   ├── README.md
│   ├── cellsegmentation
│   │   ├── __pycache__
│   │   │   ├── objseg.cpython-311.pyc
│   │   │   ├── objseg.cpython-312.pyc
│   │   │   └── objseg.cpython-38.pyc
│   │   ├── blur_detect.py
│   │   ├── cell-segmente-on-ssDNA.ipynb
│   │   ├── cellsegment.py
│   │   ├── clahe_enhance.py
│   │   ├── default.cppipe
│   │   ├── demo.tif
│   │   ├── objseg.py
│   │   └── segarr.py
│   ├── datapreprocessing
│   │   ├── GEM3D_toolkit.py
│   │   ├── README.md
│   │   ├── affine_gem.py
│   │   ├── affine_h5ad.py
│   │   ├── affine_ssdna.py
│   │   ├── affine_txt.py
│   │   ├── apply_cells.py
│   │   ├── apply_registration.py
│   │   ├── chop_gem.py
│   │   ├── chop_image.py
│   │   ├── chop_paste.py
│   │   ├── draw_heatmap.py
│   │   ├── gem_to_h5ad.py
│   │   ├── gem_xy.py
│   │   ├── h5ad_dataframe.py
│   │   ├── image_blend.py
│   │   ├── mask_gem.py
│   │   ├── mask_h5ad.py
│   │   ├── merge_h5ad.py
│   │   ├── save_miscdf.py
│   │   ├── slice_dataframe.py
│   │   ├── split_gem.py
│   │   └── trakEM2_to_affine.py
│   ├── mirrorregistration
│   │   ├── MIRROR.py
│   │   ├── README.md
│   │   ├── gem_to_gemc.py
│   │   ├── gemc_to_h5ad.py
│   │   ├── prepare_registration_heatmap.py
│   │   ├── prepare_registration_ssdna.py
│   │   ├── second_registration.py
│   │   └── slice_dataframe.py
│   └── seamalignment
│       ├── README.md
│       ├── SEAM.py
│       ├── apply_alignment.py
│       ├── get_xml_matrix.py
│       └── prepare_alignment_image.py
├── docs
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
├── planarian manuscript
│   └── mappingtoplanarian.jpg
├── pyproject.toml
├── requirements.txt
└── tree.md

28 directories, 154 files
