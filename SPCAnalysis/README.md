# SPC Analysis
On the basic assumption that a functional unit must be comprised of a group of cells rather than a single cell, and in the case of neoblast cells in planarian dispersedly distributed in space, we used a functional-spatial-proximity-based strategy to define hubs for cells in similar states or with similar functions.


This section mainly contains:

## SPC Clustering (spcclustering) 
For advanced usage and more demo example, please see
[here](https://github.com/BGI-Qingdao/SPC)

![mouse brain demo](https://github.com/BGI-Qingdao/SPC/blob/main/demo/compare.jpg)


## Meta-niche Analysis (meta-nicheanalysis)

See the detailed [tutorial](https://github.com/BGI-Qingdao/4D-BioReconX/blob/main/SPCAnalysis/Meta-Niche/calculate-niches-on-coords.ipynb)

## Cell-cell colocalization (cell-cellcolocalization)
### Brief

Calcaulate cell-cell colocation array by KL in 3D coordinate and then generate its MST.

### Usage

```
Usage: CellColocation.py -i <in.h5ad>
                         -o <prefix>
                         -b [binsize (for each slice), default 50 (refer to 50um if 1 unit=1um).]
                            Please set -b based on your coordinate system !!!
                            Notice: -b must be integer
                         -m [min cell number for a cell type, default 100.]
                         -f [sample fraction, default 0.8. The bootstrap fraction.]
                         -l [loop number, default 100. The number of iterations for bootstrap.]
                         -s [spatial key in obsm, default 'spatial']
                         -a [annotation key in obs, default 'annotation']
```

*  Ref

Wei, R., He, S., Bai, S. et al. Spatial charting of single-cell transcriptomes in tissues. Nat Biotechnol 40, 1190–1199 (2022). https://doi.org/10.1038/s41587-022-01233-1


## Digital Lineage Tracing (digitallineagetracing)

### workflow

* step01: foreach pair of timepoint, a) calulate the simularity array for each pair of cells between thoes two datasets and b) save the closest cell in another dataset for every cells.


```
Usage   : python3 s01.CellCor.py  <-i in.lst>
                                  [-o output folder, default ./]

example of in.lst
-------------  
time0
time1
time2
time3
-------------

and we will load time0.h5ad, time1.h5ad, time2.h5ad, time3.h5ad.

```

* step02: foreach pair of timepoint, calculate the mapping number of cells for all clusters.

```
Usage: python3 s02.get_number_anno.py  <-i in.lst>
                                      [-o output folder, default ./]

we will load time0.h5ad, time1.h5ad, time2.h5ad, time3.h5ad.
and data generate by s01 in output folder!
```

* step03: filter weak connections.

```
Usage: python3 s02.filter.py [-f filter threshold, default 0.05]
                                      [-o output folder, default ./]

we will load data generate by s02 in output folder!

```

* step04: now we plot the final sankey plot.

See our html for visualization codes and results!



## SPC-based Gene Regulatory Newtwork Inference
For advanced usage and more demo example, please see
[SpaGRN] (https://github.com/BGI-Qingdao/SpaGRN))

![heatmap of cell-specific regulons](https://github.com/BGI-Qingdao/SpaGRN/blob/main/resource/E14-16h_hotspot_clusters_heatmap_top5.png)

![2D spatial distribution map of a regulon](https://github.com/BGI-Qingdao/SpaGRN/blob/main/resource/Egr3.png)
 
![3D spatial distribution map of a regulon](https://github.com/BGI-Qingdao/SpaGRN/blob/main/resource/grh_L3.png)


## SPC-based Cell-cell Communications (spc-basedccc)

We modified [CellChat](https://github.com/sqjin/CellChat) and [NICHES](https://github.com/msraredon/NICHES) for the 3D database.