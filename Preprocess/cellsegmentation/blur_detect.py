
import numpy as np
import cv2

def blur_detect(image, channel='g', chunk_size=3, method='laplacian', top_svd=30, 
        outfile=None, show_in_rgb=None, show_in_grey=None):
    """
    Calculte blur values with stepwise slide chunks for RGB image

    Parameters
    ------------------------------
    image: np.ndarray, image
        image matrix with three channels
    channel: {'r', 'g', 'b'}, default g
        which channel to be used
    chunk_size: int
        pixel number for each chunk
    method: {'laplacian', 'svd'}, default laplacian
        which method to calculate blur value
    top_svd: int
        top N svd used for svd method
    outfile: str
        write the blur matrix into file
    show_in_rgb: str
        display the blur value in rgb image
    show_in_grey: str
        display the blur value in grey image
    """
    
    # background was detected as blur region
    # I need to segmentate tissue region firstly
    # here I used color masking for segmentation on green channel
    b, g, r = cv2.split(image)
    # detect based on green channel
    light = 10
    dark = 255
    if channel == 'r':
        channel = r
    elif channel == 'g':
        channel = g
    elif channel == 'b':
        channel = b
    mask = cv2.inRange(channel, light, dark)
    kernel = np.ones((10, 10), np.uint8)
    mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, kernel)

    blur_image = np.zeros(shape=image.shape, dtype=np.uint8)
    for (x, y), value in np.ndenumerate(mask):
        if value == 0:
            continue
        chunk = image[x:x+chunk_size, y:y+chunk_size]
        # small value indicate blur region
        if method == 'laplacian':
            blur_value = cv2.Laplacian(chunk, cv2.CV_64F).var()
        elif method == 'svd':
            u, sigma, vt = np.linalg.svd(img)
            blur_value = sum(sigma[:top_svd]) / sum(sigma)
        blur_image[x, y] = blur_value
    
    if outfile:
        np.savetxt(outfile, blur_image, fmt='%d')
    
    if show_in_rgb:
        blur_rgb_image = cv2.applyColorMap(blur_image, cv2.COLORMAP_JET)
        cv2.imwrite(show_in_rgb, blur_rgb_image)
    
    if show_in_grey:
        black = np.zeros(shape=image.shape, dtype=np.uint8)
        blur_mask = np.where(blur_image < 30, mask, black)
        cv2.imwrite(show_in_grey, blur_mask)

    return blur_image


