
import numpy as np
import cv2

def clahe(image):
    #-----Reading the image-----------------------------------------------------
    if not isinstance(image, np.ndarray):
        image = cv2.imread(image, 1)

    #-----Converting image to LAB Color model----------------------------------- 
    lab = cv2.cvtColor(image, cv2.COLOR_BGR2LAB)

    #-----Splitting the LAB image to different channels-------------------------
    l, a, b = cv2.split(lab)

    #-----Applying CLAHE to L-channel-------------------------------------------
    clahe = cv2.createCLAHE(clipLimit=2, tileGridSize=(8,8))
    cl = clahe.apply(l)

    #-----Merge the CLAHE enhanced L-channel with the a and b channel-----------
    limg = cv2.merge((cl,a,b))

    #-----Converting image from LAB Color model to RGB model--------------------
    final = cv2.cvtColor(limg, cv2.COLOR_LAB2BGR)

    #_____END_____#
    #return cl
    return final

if __name__ == '__main__':

    import sys
    from pathlib import Path

    image = sys.argv[1]
    out = Path(image).name
    
    for i in range(5):
        image = clahe(image)

    cv2.imwrite(out, image)

