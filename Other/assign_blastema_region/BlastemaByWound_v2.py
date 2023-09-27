#!/usr/bin/env python3

import sys
import getopt
import numpy as np
import pandas as pd
from skimage import io as skio

########################################################
# main logic function
#
def get_blastema(sample,prefix,ll,lr,rl,rr,exponential_number,wound_only):
    base_file = f'/dellfsqd2/ST_OCEAN/USER/guolidong/wochong_3d/blastema/wound_info/{sample}/{sample}_photo_bin3space_masked.tif'
    left_edge_file = f'/dellfsqd2/ST_OCEAN/USER/guolidong/wochong_3d/blastema/wound_info/{sample}/wound_left.tif'
    right_edge_file = f'/dellfsqd2/ST_OCEAN/USER/guolidong/wochong_3d/blastema/wound_info/{sample}/wound_right.tif'
    pos_file = f'/dellfsqd2/ST_OCEAN/USER/guolidong/wochong_3d/blastema/raw_pos/PosDB.{sample}.txt'

    base = skio.imread(base_file)
    y, x ,z = np.nonzero(base)
    ymin = np.min(y)
    ymax = np.max(y)
    xmax = np.max(x)
    H,W,D = base.shape
    
    yvals = np.array(range(ymin, ymax+1),dtype=int)
    def mask_one_edge( fname , yvals , ymin,exponential_number,xmax):
        ########################################################
        # get edge points
        #
        masked_edge_data = skio.imread(fname)
        y, x = np.nonzero(masked_edge_data)
        fit_ymin = np.min(y)
        fit_ymax = np.max(y)
        ########################################################
        # fit edge line
        #
        # let x->y y->x to fit vertical line
        flist = np.polyfit(y, x, exponential_number)
        ########################################################
        # plot final line
        #
        #fit_yvals = np.array(range(fit_ymin, fit_ymax+1),dtype=int)
        #xvals = np.polyval(flist, fit_yvals)
        #xret = np.zeros(len(yvals),dtype=int)
        #for i in range(len(xvals)):
        #    xret[fit_ymin-ymin+i] = xvals[i]
        xret = np.polyval(flist, yvals).astype(int)
        xret [ xret < 0 ] = 0
        xret [ xret > xmax ] = xmax
        return flist, xret , y , x
    
    
    flist ,  x_l , sy_l, sx_l = mask_one_edge(left_edge_file, yvals,ymin, exponential_number ,W-1)
    #print(f'the left polyfit parameter is : {flist}',flush=True)
    linedf = pd.DataFrame()
    linedf['x'] = x_l
    linedf['y'] = yvals
    linedf['x'] = linedf['x']*3
    linedf['y'] = linedf['y']*3
    #linedf['y'] = yvals
    linedf.to_csv(f'{prefix}.left_line.csv',header=True,sep='\t',index=False)
    flist ,  x_r , sy_r, sx_r = mask_one_edge(right_edge_file, yvals,ymin ,exponential_number,W-1)
    #print(f'the right polyfit parameter is : {flist}',flush=True)
    linedf = pd.DataFrame()
    linedf['x'] = x_l
    linedf['y'] = yvals
    linedf['x'] = linedf['x']*3
    linedf['y'] = linedf['y']*3
    #linedf['y'] = yvals
    linedf.to_csv(f'{prefix}.right_line.csv',header=True,sep='\t',index=False)
    raw = base.copy()
    base[sy_l , sx_l, 1] = 255 
    base[sy_r , sx_r, 1] = 255 
    base[yvals, x_l, 0] = 255
    base[yvals, x_r, 0] = 255
    skio.imsave(f'{prefix}.edge.tif',base.astype('uint8'))
     
    labels = np.zeros((H,W),dtype='uint8')   
    label_array = [0,1,2,3,4,5,6,7]
    if wound_only:
        label_array = [0,3,3,3,4,5,5,5]
    for i in range(len(yvals)):
        y_index = yvals[i]-ymin
        if x_l[y_index] != 0 and  x_r[y_index] != 0 :
            labels[yvals[i],  0                 : x_l[y_index] - ll ] = label_array[1]
            labels[yvals[i],  x_l[y_index] - ll : x_l[y_index]      ] = label_array[2]
            labels[yvals[i],  x_l[y_index]      : x_l[y_index] + lr ] = label_array[3]
            labels[yvals[i],  x_l[y_index] + lr : x_r[y_index] - rl ] = label_array[4]
            labels[yvals[i],  x_r[y_index] - rl : x_r[y_index]      ] = label_array[5]
            labels[yvals[i],  x_r[y_index]      : x_r[y_index] + rr ] = label_array[6]
            labels[yvals[i],  x_r[y_index] + rr :                   ] = label_array[7]
        elif x_l[y_index] != 0 and  x_r[y_index] == 0 :                               
            labels[yvals[i],  0                 : x_l[y_index] - ll ] = label_array[1]
            labels[yvals[i],  x_l[y_index] - ll : x_l[y_index]      ] = label_array[2]
            labels[yvals[i],  x_l[y_index]      : x_l[y_index] + lr ] = label_array[3]
            labels[yvals[i],  x_l[y_index] + lr :                   ] = label_array[4]
        elif x_l[y_index] == 0 and  x_r[y_index] != 0 :                               
            labels[yvals[i],                    : x_r[y_index] - rl ] = label_array[4]
            labels[yvals[i],  x_r[y_index] - rl : x_r[y_index]      ] = label_array[5]
            labels[yvals[i],  x_r[y_index]      : x_r[y_index] + rr ] = label_array[6]
            labels[yvals[i],  x_r[y_index] + rr :                   ] = label_array[7]
        else:                                                                         
            labels[yvals[i],                    :                   ] = label_array[4]
    
    label_all_df = pd.DataFrame()
    #raw = base.copy()
    base[:,:,:]=0
    
    label = np.where(labels==1)
    base[label[0],label[1],0] = 255 #draw 1 as red

    label_df = pd.DataFrame()
    label_df['y'] = label[0]
    label_df['x'] = label[1]
    label_df['l'] = np.ones(len(label[0])) *1
    label_all_df = label_df  
   
    label = np.where(labels==2)
    base[label[0],label[1],1] = 255 #draw 2 as green

    label_df = pd.DataFrame()
    label_df['y'] = label[0]
    label_df['x'] = label[1]
    label_df['l'] = np.ones(len(label[0])) *2
    label_all_df = pd.concat([ label_all_df, label_df ], ignore_index=True)
    
    label = np.where(labels==3)
    base[label[0],label[1],0] = 255 #draw 3 as magenta 
    base[label[0],label[1],2] = 255 #draw 3 as magenta 
    label_df = pd.DataFrame()
    label_df['y'] = label[0]
    label_df['x'] = label[1]
    label_df['l'] = np.ones(len(label[0])) *3
    label_all_df = pd.concat([ label_all_df, label_df ], ignore_index=True)
    
    label = np.where(labels==4)
    base[label[0],label[1],0] = 255 #draw 4 as yellow
    base[label[0],label[1],1] = 255 #draw 4 as yelow
    label_df = pd.DataFrame()
    label_df['y'] = label[0]
    label_df['x'] = label[1]
    label_df['l'] = np.ones(len(label[0])) *4
    label_all_df = pd.concat([ label_all_df, label_df ], ignore_index=True)
    
    label = np.where(labels==5)
    base[label[0],label[1],:] = 255 #draw 5 as white
    raw[y,x,:] = base[y,x,:]
    label_df = pd.DataFrame()
    label_df['y'] = label[0]
    label_df['x'] = label[1]
    label_df['l'] = np.ones(len(label[0])) *5
    label_all_df = pd.concat([ label_all_df, label_df ], ignore_index=True)
    
    label = np.where(labels==6)
    base[label[0],label[1],1] = 255 #draw 6 as cyan 
    base[label[0],label[1],2] = 255 #draw 6 as cyan 
    label_df = pd.DataFrame()
    label_df['y'] = label[0]
    label_df['x'] = label[1]
    label_df['l'] = np.ones(len(label[0])) *6
    label_all_df = pd.concat([ label_all_df, label_df ], ignore_index=True)

    label = np.where(labels==7)
    base[label[0],label[1],0] = 255 #draw 7 as orange
    base[label[0],label[1],1] = 69  #draw 7 as orange
    label_df = pd.DataFrame()
    label_df['y'] = label[0]
    label_df['x'] = label[1]
    label_df['l'] = np.ones(len(label[0])) *7
    label_all_df = pd.concat([ label_all_df, label_df ], ignore_index=True)


    raw[y,x,:] = base[y,x,:]
    skio.imsave(f'{prefix}.class.tif',raw.astype('uint8'))

    pos_file_df = pd.read_csv(pos_file,sep=',',header=None)
    pos_file_df.columns = ['label','cell','x','y','z']
    pos_file_df['x'] = pos_file_df['x'] / 3
    pos_file_df['x'] = pos_file_df['x'].astype(int)   
    pos_file_df['y'] = pos_file_df['y'] / 3
    pos_file_df['y'] = pos_file_df['y'].astype(int)   
    
    output = pos_file_df.merge(label_all_df, on=['x' ,'y'])
    output['l'] = output['l'].astype(int)
    output.to_csv(f'{prefix}.class.csv',sep='\t',columns=['label','cell','l'],index=None) 

########################################################
# usage
#
def usage():
    print("""
Usage	: python3 BlastemaByWound.py < -s sample_name> 
                                     [ -o output prefix, default output]
                                     [ -e exponential number, default 2]
                                     [--only_wound yes/no, default no]
                                     [--ll left wound left extern distance, default 20]
                                     [--lr left wound right extern distance, default 20]
                                     [--rl right wound left extern distance, default 20]
                                     [--rr right wound right extern distance, default 20]

Notice	: the unit of distance is 3 micron, so the default 10 refer to 60 microns.

Example : 
	  example 01: python3 BlastemaByWound.py -s 12hpa1 
	  example 02: python3 BlastemaByWound.py -s WT -o test_WT
	  example 02: python3 BlastemaByWound.py -s 5dpa1 -o test_5dpa1 -e 3
          example 03: python3 BlastemaByWound.py -s 3dpa1 -o test_3dpa1_lr15 --lr 15

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
""",flush=True)

########################################################
# main
#
def main(argv):
    ########################
    # no args equal to -h
    if len(argv) == 0 :
        usage()
        sys.exit(0)

    ########################
    # default values
    sample = ''
    prefix = 'output'
    ll = lr = rl = rr = 20
    exponential_number = 2
    wound_only = False
    ########################
    # parse args
    try:
        opts, args = getopt.getopt(argv,"hs:o:",["help",
                                                 "ll=",
                                                 "lr=",
                                                 "rl=",
                                                 "rr=",
                                                 "only_wound=",
                                                 ])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif opt in ("-s" ):
            sample = arg   
        elif opt in ("-o" ):
            prefix = arg
        elif opt in ("--ll" ):
            ll = int(arg)
        elif opt in ("--lr" ):
            lr = int(arg)
        elif opt in ("--rl" ):
            rl = int(arg)
        elif opt in ("--rr" ):
            rr = int(arg)
        elif opt in ("--only_wound"):
            if arg == 'yes':
                wound_only = True

    ########################
    # sanity check
    if sample not in [ '0hpa1','0hpa2','10dpa1','10dpa2','12hpa1','12hpa2','14dpa1','14dpa2','36hpa1', '36hpa2','3dpa1', '3dpa2', '5dpa1','5dpa2', '7dpa1' ,'7dpa2','WT' ]:
        print(f'Error : invalid sample name : \" -s {sample}\"',flush=True)
        usage()
        sys.exit(1)
     
    ########################
    # do the job 
    get_blastema(sample,prefix,ll,lr,rl,rr,exponential_number,wound_only)          

if __name__ == "__main__":
    main(sys.argv[1:])

