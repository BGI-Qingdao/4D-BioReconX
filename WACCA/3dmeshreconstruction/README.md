# Mesh generating pipeline

## workflow
![image](https://github.com/BGI-Qingdao/4D-BioReconX/assets/8720584/0e550dc1-1faf-4dd5-b2f0-e1a3af632b74)

## usage

### gen stacked tif 

```
Usage   : python3 gen_stacked_tif.py < -i mask.lst>
                                     < -a anno.txt>
                                     < -o output prefix>
                                     [ -b binsize default 20]

```
* mask.lst example

the first column is z value ( must start from 0 ),the second column is cellbin mask file
```
0 s0.mask
1 s1.mask
...
n sn.mask
```

* anno.txt example

```
slice_id,cell_id,anno_id
1,11,1
2,12,2
```
the slice id start from 1

### reset obj

```
Usage   : python3 reset_obj.py  < -i in.obj>
                                 -o output prefix>
                                [ -b binsize default 20]
```

