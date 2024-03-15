import numpy as np
import cv2
import netCDF4 as nc    
from skimage import exposure
import matplotlib.pyplot as plt

    

def compressImg(img, map = 'gray', sat = 1e-2, gamma = 1, pixsize = [1,1], bit = 8):
    """Adjusts the image based on given saturation and gamma levels, before compressing it to either 8-bit or 16-bit, and finally resizing based on given pixel size.
    If bit specified is neither 8 nor 16, defaults to 8-bit. Any dummy data in the image is converted to value of -1.

    Args:
        img: the complex image
        map: colourscale of the compressed image. Defaults to 'gray'.
        sat: saturation level of the compressed image. Defaults to 1e-3.
        gamma: gamma level of the compressed image. Defaults to 1.
        pixsize: ratio of row and column pixels. Defaults to [1,1].
        bit: bit compression level. Defaults to 8.

    Returns:
        img_resized: The final compressed and resized image.
    """
    #find the magnitude of each complex value
    if np.any(np.iscomplex(img)):
        img = abs(img)
    #find the min and max intensity values using saturation percentage
    # low, high = exposure.rescale_intensity(img, (sat, 1-sat)).min(), \
    # exposure.rescale_intensity(img, (sat, 1-sat)).max()
    row, col = np.where(img >= 9999)
    low, high = np.percentile(img, (sat)*100), np.percentile(img, (1-sat)*100)
    # print(low, high)
    #adjust saturation and gamma and output to [0,1] range
    img_adjusted = exposure.rescale_intensity(img, in_range = (low,high),out_range= (0,1))
    img_adjusted = exposure.adjust_gamma(img_adjusted, gamma)

    # fit to bit requirements
    if bit == 16:
        img_compressed = (img_adjusted*65535).astype(float)
    #if not 16bit, defaults to 8bit
    else:
        img_compressed = (img_adjusted*255).astype(float)
    img_compressed[row,col] = -1
    #scale by pixel ratio
    img_resized = cv2.resize(img_compressed, (img_compressed.shape[1]*pixsize[1], img_compressed.shape[0]*pixsize[0]), interpolation= cv2.INTER_NEAREST)
    # plt.imshow(img_resized, cmap=map)
    # plt.show()
    return img_resized



if __name__ == '__name__':
    #open and read the file
    filePath = 'SAR_CPLX_20190823071330_9.6G_VH_11_pres_2_fdc_246.sar.rgo.sig.nc'
    file = nc.Dataset(filePath, 'r')
    realimg = file.variables['SigmaImageSingleLookRealPart'][:]
    imagimg = file.variables['SigmaImageSingleLookImaginaryPart'][:]
    GBPGridInfo = file.variables['GBPGridInfo'][:]
    ModelTransformationTag = file.variables['ModelTransformationTag'][:]
    lat = file.variables['LatImage'][:]
    long = file.variables['LonImage'][:]
    realimg = realimg[:] 
    imagimg = imagimg[:]
    calimage = file.variables['CalImage'][:]
    fullimg = realimg + 1j*imagimg

    #compress and return 8bit image
    fullimg = compressImg(fullimg)
    plt.imshow(fullimg)
    plt.show()

