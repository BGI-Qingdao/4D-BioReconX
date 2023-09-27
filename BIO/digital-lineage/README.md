# the Digital lineage tracing pipeline

## workflow

* step01: foreach pair of timepoint, a) calulate the simularity array for each pair of cells between thoes two datasets and b) save the closest cell in another dataset for every cells.
* step02: foreach pair of timepoint, calculate the mapping number of cells for all clusters.
* step03: filter weak connections.
* step04: for plot sankey, we create pseudo cluster to accept the extra cells from previous timepoint.
* step05: now we plot the final sankey plot.
