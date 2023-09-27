# the Digital lineage tracing pipeline

## workflow

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
