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
Usage : GEM_toolkit.py action [options]

Actions:

---------------------------------------------------------------------

 Format coverting tools:
    gem_to_h5ad                   convert GEM into h5ad by a certain binsize.
 
 Affine tools:
    affine_gem                    modify the 2D coordinate in GEM(C) by user-defined affine matrix.
    affine_h5ad                   modify the 2D coordinate in GEM(C) by user-defined affine matrix.
    affine_ssdna                  affine the ssdna image by user-defined affine matrix.
    affine_txt                    affine txt like cell.mask by user-defined affine matrix.
    apply_registration            use registration result(with/without ROI) to update ssdna/mask/gem ...
    apply_cells               	  add cells column to gem based on registered mask file.

 Region of interest(ROI) tools: 
    chop_image                    chop region of interests from whole image.
    chop_gem                      chop region of interests from GEM(C).

 Mask tools:
    mask_gem                      mask GEM(C) by mask image.
    mask_h5ad                     mask h5ad data by mask image.
  
 Visualization tools:
    draw_heatmap                  draw heatmap of expression counts in bin1 resolution with/without cellbin and with/without ssDNA.
    image_blend                   merge image(like heatmap/annotation image) with ssDNA and border image
 
 Other tools:
    chop_paste                    chop or paste ssDNA image. This tools is useful for ultra-large ssDNA image.
    trakEM2_to_affine             covert trakEM2_matrix to standart affine matrix.
    get_xml_matrix                get matrix from trakEM2 xml file.
    split_gem                     split gem by x or y coordinate.
    merge_h5ad                    merge files of h5ad.
    gem_xy                        get xmin ymin of gem

    -----------------------------------------------------------------
    -h/--help               show this short usage
""")

# logic codes
if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] in ( "-h" , "--help" ):
        main_usage()
        exit(0)
    elif len(sys.argv) < 2 or not sys.argv[1] in ( "second_registration",
                                                   "secondregistration",
                                                   "prepareregistrationheatmap",
                                                   "prepare_registration_heatmap",
                                                   "prepare_registration_ssdna",
                                                   "prepareregistrationssdna",
                                                   "prepareregistrationdapi",
                                                   "gem_to_cfm",
                                                   "gem_to_gemc",
                                                   "gemc_to_h5ad",
                                                   "gem_to_h5ad",
                                                   "gem_xy",
                                                   "chop_image",
                                                   "chop_gem",
                                                   "mask_gem",
                                                   "affine_gem",
                                                   "affine_ssdna",
                                                   "affine_txt",
                                                   "affine_h5ad",
                                                   "mask_h5ad",
                                                   "prepare_alignment_image",
                                                   "apply_alignment",
                                                   "apply_registration",
                                                   "apply_cells",
                                                   "trakEM2_to_affine",
                                                   "chop_paste",
                                                   "image_blend",
                                                   "draw_heatmap",
                                                   "get_xml_matrix",
                                                   "split_gem",
                                                   "merge_h5ad"
                                                   ) :
        main_usage()
        exit(1)
    if  sys.argv[1] == "chop_image" :
        from chop_image import chopimages_main
        chopimages_main(sys.argv[2:])
        exit(0)
    elif  sys.argv[1] == "chop_gem" :
        from chop_gem import chopgems_main
        chopgems_main(sys.argv[2:])
        exit(0)
    elif  sys.argv[1] == "mask_gem" :
        from mask_gem import mask_gem_main 
        mask_gem_main(sys.argv[2:])
        exit(0)
    elif  sys.argv[1] == "gem_to_h5ad" :
        from gem_to_h5ad import gem_to_h5ad_main 
        gem_to_h5ad_main(sys.argv[2:])
        exit(0)
    elif  sys.argv[1] == "gemc_to_h5ad" :
        from gemc_to_h5ad import gemc_to_h5ad_main 
        gemc_to_h5ad_main(sys.argv[2:])
        exit(0)
    elif sys.argv[1] == "affine_ssdna" :
        from affine_ssdna import affine_ssdna_main
        affine_ssdna_main(sys.argv[2:])
        exit(0)
    elif sys.argv[1] == "affine_gem" :
        from affine_gem import affine_gem_main
        affine_gem_main(sys.argv[2:])
        exit(0)
    elif sys.argv[1] == "affine_txt" :
        from affine_txt import affine_txt_main
        affine_txt_main(sys.argv[2:])
        exit(0)
    elif sys.argv[1] == "affine_h5ad" :
        from affine_h5ad import affine_h5ad_main
        affine_h5ad_main(sys.argv[2:])
        exit(0)
    elif sys.argv[1] == "mask_h5ad" :
        from mask_h5ad import mask_h5ad_main
        mask_h5ad_main(sys.argv[2:])
        exit(0)
    elif sys.argv[1] == "chop_paste" :
        from chop_paste import chop_paste_main
        chop_paste_main(sys.argv[2:])
        exit(0)
    elif sys.argv[1] == "trakEM2_to_affine" :
        from trakEM2_to_affine import trakEM2_to_affine_main
        trakEM2_to_affine_main(sys.argv[2:])
        exit(0)
    elif sys.argv[1] == "image_blend" :
        from image_blend import image_blend_main
        image_blend_main(sys.argv[2:])
        exit(0)
    elif sys.argv[1] == "draw_heatmap" :
        from draw_heatmap import heatmap_main
        heatmap_main(sys.argv[2:])
        exit(0) 
    elif sys.argv[1] == "split_gem" :
        from split_gem import split_gem_main
        split_gem_main(sys.argv[2:]) 
    elif sys.argv[1] == "merge_h5ad" :
        from merge_h5ad import merge_h5ad_main
        merge_h5ad_main(sys.argv[2:])
        exit(0)
    elif sys.argv[1] == "apply_registration" :
        from apply_registration import apply_registration_main
        apply_registration_main(sys.argv[2:])
        exit(0)
    elif sys.argv[1] == "gem_xy" :
        from gem_xy import gem_xy_main
        gem_xy_main(sys.argv[2:])
        exit(0)
    elif sys.argv[1] == "apply_cells" :
        from apply_cells import apply_cells_main
        apply_cells_main(sys.argv[2:])
        exit(0)
    else:
        main_usage()
        exit(1)
