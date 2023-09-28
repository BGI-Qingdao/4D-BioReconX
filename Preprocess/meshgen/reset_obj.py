import os,sys,json
import numpy as np
import pandas as pd
from io import StringIO
def usage():
    print("""
Usage   : python3 reset_obj.py  < -i in.obj>
                                 -o output prefix>
                                [ -b binsize default 20]
""")


def MesheStr(objfile):
    mesh_str = ''
    file1 = open(objfile, 'r')
    Lines = file1.readlines()
    file1.close()
    for line in Lines:
        if len(line)>2 and ( line[0] == 'v' or line[0] == 'f' ) and line[1] == ' ':
            mesh_str = mesh_str + line
    return mesh_str

def add_mesh_str(objfile,binsize,outfile):
    mesh_io = StringIO(MesheStr(objfile))
    cache = pd.read_csv(mesh_io, sep='\s+',header=None, compression='infer', comment='#')
    cache.columns = ['type','v1','v2','v3']
    vectors = cache[cache['type'] == 'v'].copy()
    vectors = vectors[['v1','v2','v3']].copy()
    vectors.columns = ['x','y','z']
    #print(vectors)
    vectors = vectors.astype(float)
    vectors['x'] = vectors['x'] * binsize
    vectors['y'] = vectors['y'] * binsize
    vectors['z'] = vectors['z'] + 1
    vectors['z'] = (vectors['z'] -40 )*20
    vectors['z'] = vectors['z'] * -1
    vectors = vectors.astype(int)
    vectors['type'] = 'v'
    faces = cache[cache['type'] == 'f'].copy()
    if faces.dtypes['v1'] == object:
        faces['i'] = faces.apply(lambda row: int(row['v1'].split('/')[0]), axis=1)
        faces['j'] = faces.apply(lambda row: int(row['v2'].split('/')[0]), axis=1)
        faces['k'] = faces.apply(lambda row: int(row['v3'].split('/')[0]), axis=1)
    else:
        faces['i'] = faces['v1']
        faces['j'] = faces['v2']
        faces['k'] = faces['v3']
    faces = faces[['i','j','k']].copy()
    faces = faces.astype(int)
    faces['type'] = 'f'
    faces.columns = ['x','y','z','type']
    vectors = vectors[['type','x','y','z']].copy()
    faces = faces[['type','x','y','z']].copy()
    obj_new = pd.concat([vectors,faces],ignore_index=True)
    obj_new.to_csv(f'{outfile}.obj',sep=' ',header=False,index=False)


def main(argv):
    ########################
    # no args equal to -h
    if len(argv) == 0 :
        usage()
        sys.exit(0)

    ########################
    # default values
    inobj = ''
    binsize = 20
    prefix = ''
    ########################
    # parse args
    try:
        opts, args = getopt.getopt(argv,"hi:o:b:",["help"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif opt in ("-i" ):
            inobj = arg
        elif opt in ("-o" ):
            prefix = arg
        elif opt in ("-b" ):
            binsize = int(arg)

    add_mesh_str(inobj,binsize,prefix)


if __name__ == "__main__":
    main(sys.argv[1:])


