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
Usage : MIRROR.py action [options]

Actions:

---------------------------------------------------------------------
    prepare_registration_heatmap  generate 8bit spot-level heatmap with highlighted tracklines.
    prepare_registration_ssdna    generate 8bit close-spot-level ssDNA iamge with highlighted tracklines.
    second_registration           second round registration.
    gem_to_gemc                   convert GEM into GEMC based on cellbin result and registration results.
    gemc_to_h5ad                  convert GEMC into h5ad.
 
    -----------------------------------------------------------------
    -h/--help               show this short usage
""")

# logic codes
if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] in ( "-h" , "--help" ):
        main_usage()
        exit(0)
    elif len(sys.argv) < 2 or not sys.argv[1] in ( "second_registration",
                                                   "prepare_registration_heatmap",
                                                   "prepare_registration_ssdna",
                                                   "gem_to_gemc",
                                                   "gemc_to_h5ad",
                                                   ) :
        main_usage()
        exit(1)
    elif sys.argv[1] in ( "secondregistration" , "second_registration"):
        from second_registration import secondregistration_main
        secondregistration_main(sys.argv[2:])
        exit(0)
    elif sys.argv[1] in ("prepareregistrationheatmap","prepare_registration_heatmap"):
        from prepare_registration_heatmap import prepareregistrationheatmap_main 
        prepareregistrationheatmap_main(sys.argv[2:])
        exit(0)
    elif sys.argv[1] in ("prepareregistrationssdna" , "prepareregistrationdapi","prepare_registration_ssdna"):
        from prepare_registration_ssdna import prepareregistrationdapi_main 
        prepareregistrationdapi_main(sys.argv[2:])
        exit(0)
    elif sys.argv[1] in ("gem_to_cfm", "gem_to_gemc"):
        from gem_to_gemc import gem_to_cfm_main
        gem_to_cfm_main(sys.argv[2:])
        exit(0)
    elif  sys.argv[1] == "gemc_to_h5ad" :
        from gemc_to_h5ad import gemc_to_h5ad_main 
        gemc_to_h5ad_main(sys.argv[2:])
        exit(0)
    else:
        main_usage()
        exit(1)
