## Brief

Calcaulate cell-cell colocation array by KL in 3D coordinate and then generate its MST.


## Usage

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

## Ref

Wei, R., He, S., Bai, S. et al. Spatial charting of single-cell transcriptomes in tissues. Nat Biotechnol 40, 1190–1199 (2022). https://doi.org/10.1038/s41587-022-01233-1