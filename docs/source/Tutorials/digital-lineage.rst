.. _`digital-lineage`:
========================================
Digital Lineage Tracing Pipeline
========================================

workflow
========================================

step01: for each pair of time points, a) calculate the similarity array for each pair of cells between those two datasets and b) save the closest cell in another dataset for every cells.

step02: for each pair of time points, calculate the mapping number of cells for all clusters.

step03: filter weak connections.

step04: for plot Sankey, we create a pseudo cluster to accept the extra cells from the previous timepoint.

step05: now we plot the final Sankey plot.
