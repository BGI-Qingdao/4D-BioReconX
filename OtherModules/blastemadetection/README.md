# Detect blastema by pigment difference

## workflow
![image](https://github.com/BGI-Qingdao/4D-BioReconX/assets/8720584/fb2c54a3-0eff-4d92-99bd-f76516a1898f)

## example images

* raw photo
  
![image](https://github.com/BGI-Qingdao/4D-BioReconX/assets/8720584/96c193ee-9044-4973-9516-853098897a4e)

* 2D mapping photo of 3D atlas

![image](https://github.com/BGI-Qingdao/4D-BioReconX/assets/8720584/2859f355-ed30-4ec4-917c-de96f00f2339)

* aligned photo and detect wound mask

![image](https://github.com/BGI-Qingdao/4D-BioReconX/assets/8720584/57d56423-66c2-459d-9cc8-ef2c2fc912f8)

* now assign region to each cells

```
./BlastemaByWound_v2.py

Usage   : python3 BlastemaByWound.py < -p prefix>
                                     [ -o output prefix, default output]
                                     [ -e exponential number, default 2]
                                     [--only_wound yes/no, default no]
                                     [--ll left wound left extern distance, default 20]
                                     [--lr left wound right extern distance, default 20]
                                     [--rl right wound left extern distance, default 20]
                                     [--rr right wound right extern distance, default 20]

Notice  : the unit of distance is 3 micron, so the default 10 refer to 60 microns.

Example :
          example 01: python3 BlastemaByWound.py -p 12hpa1
          example 02: python3 BlastemaByWound.py -p WT -o test_WT
          example 02: python3 BlastemaByWound.py -p 5dpa1 -o test_5dpa1 -e 3
          example 03: python3 BlastemaByWound.py -p 3dpa1 -o test_3dpa1_lr15 --lr 15

Output label :
          1 -- [red]     left blastema
          2 -- [green]   left margin of left wound
          3 -- [magenta] right margin of left wound
          4 -- [yellow]  body
          5 -- [white]   left margin of right wound
          6 -- [cyan]    right margin of right wound
          7 -- [orange]  right blastema

Output label in only_wound mode:
          3 -- [magenta] left wound region, similar to 1+2+3 in blastema mode
          4 -- [yellow]  body
          5 -- [white]   right wound region, similar to 5+6+7 in blastema mode
```

* final region image

![image](https://github.com/BGI-Qingdao/4D-BioReconX/assets/8720584/ae275b44-47d7-4ec3-bb53-64e880606a75)

