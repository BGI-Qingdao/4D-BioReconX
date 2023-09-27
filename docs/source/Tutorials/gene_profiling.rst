.. _`gene_profiling`:
========================================
PCG Gene Profiling Pipeline
========================================

.. note:: 
The basic idea here is converting 3D data into 1D pseudoBulk data.To achieve this, \
we choose a body axis ( an idx of spatial coordinate ) and group data into equal-width bins.\
Cells of each bin then merge into one meta cell.


Usage
========================================
usage: AP_Profiling [-h] -i INPUT -o OUTPUT [-m MIN_BIN] [-n NUM_HVG] [-b BIN_NUM]
  

Create AP Profiling for genes.

optional arguments
========================================
  
  -h, --help            show this help message and exit
  
  -i INPUT, --input INPUT
  
                the input adata, Notice: we expect raw count matrix!!!

  -o OUTPUT, --output OUTPUT

                the output prefix

  -m MIN_BIN, --min_bin MIN_BIN

                the minimum allowed cell number of a bin, default(50)

  -n NUM_HVG, --num_HVG NUM_HVG

                the number of HVG bin, default(5000)

  -b BIN_NUM, --bin_num BIN_NUM

                the total bin number, default(100)
