#!/usr/bin/env python3

import sys
import getopt
import time

############################################################################
# Main gateway: 
#############################################################################
# usage
def main_usage():
    print("""
Usage : SEAM.py action [options]

Actions:

---------------------------------------------------------------------
prepare_alignment_image       generate 8bit spot-level binary/annatation image for 3D alignment.
apply_alignment               set 3D coordinate for GEM(C)/h5ad/ssDNA/cell.mask.
get_xml_matrix                get matrix from trakEM2 xml file.

---------------------------------------------------------------------
-h/--help               show this short usage
""")

# logic codes
if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] in ( "-h" , "--help" ):
        main_usage()
        exit(0)
    elif len(sys.argv) < 2 or not sys.argv[1] in ( "prepare_alignment_image",
                                                   "apply_alignment",
                                                   "get_xml_matrix",
                                                   ) :
        main_usage()
        exit(1)
    if sys.argv[1] == "prepare_alignment_image" :
        from prepare_alignment_image import prepare_alignment_image_main
        prepare_alignment_image_main(sys.argv[2:])
        exit(0)
    elif sys.argv[1] == "apply_alignment" :
        from apply_alignment import apply_alignment_main
        apply_alignment_main(sys.argv[2:])
        exit(0)
    elif sys.argv[1] == "get_xml_matrix" :
        from get_xml_matrix import get_xml_matrix_main
        get_xml_matrix_main(sys.argv[2:])
        exit(0)
    else:
        main_usage()
        exit(1)
