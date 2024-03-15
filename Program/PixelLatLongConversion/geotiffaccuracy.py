import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import rasterio
import seaborn as sns
import os
from utm import to_latlon

import sys
sys.path.append(os.path.join(sys.path[0],".."))
from PixelLatLongConversion.pixel_LatLon import getLatLonFromGeotiff
from SLCtogeotiff.SLCtogeotiff import convertToGeotiff_affine
from SLCtogeotiff.adjustSLC import compressImg
from processimage.validimage import getValidImage


def getAccuracyHeatmap(filePath, gbpgridinfo, outputPath):
    img = rasterio.open(filePath)
    img_gt = img.transform
    geotifflatmatrix = np.zeros((latimage.shape[0], latimage.shape[1]))
    geotifflonmatrix = np.zeros((latimage.shape[0], latimage.shape[1]))
    for y in range(0, fullimg.shape[0]):
        for x in range(0, fullimg.shape[1]):
            #convert geotransform data to long and lat
            lon = img_gt[0] * x + img_gt[1] * y + img_gt[2]
            lat = img_gt[3] * x + img_gt[4] * y + img_gt[5]
            lat, lon = to_latlon(lon, lat, int(gbpgridinfo[3]), 'S' if gbpgridinfo[4] == 1 else 'N')
            # print(lat, lon, x, y)
            geotifflatmatrix[y][x] = lat
            geotifflonmatrix[y][x] = lon

    print(geotifflatmatrix.shape, latimage.shape)
    latdiffmatrix = geotifflatmatrix - latimage
    londiffmatrix = geotifflonmatrix - lonimage

    fig, axs = plt.subplots(1,2)
    sns.heatmap(latdiffmatrix, cbar=True, ax = axs[0])
    axs[0].set_xlabel('Column')
    axs[0].set_xlabel('Row')
    axs[0].set_title('Latitude Difference')

    sns.heatmap(londiffmatrix, cbar=True, ax = axs[1])
    axs[1].set_xlabel('Column')
    axs[1].set_xlabel('Row')
    axs[1].set_title('Longitude Difference')

    plt.tight_layout()
    plt.s
    if os.path.exists(outputPath + ".png"):
        os.remove(outputPath + ".png")
    plt.savefig(outputPath + ".png")
    plt.show()

def test_coord_accuracy(geotiffPath, x, y, latimage, lonimage):
    lat, lon = getLatLonFromGeotiff(geotiffPath,x, y) 
    imglat = latimage[y][x]
    imglon = lonimage[y][x]
    print(f'for x and y {x}, {y}, the latlon in the geotiff is {lat}, {lon} and the actual latlon in the netCDF is {imglat}, {imglon}')
    difflat = abs(lat - imglat)
    difflon = abs(lon - imglon)
    if(difflat > 0.00003 or difflon > 0.00003):
        print(f' the difference in latlon is {difflat}, {difflon}')
    
if __name__ == '__main__':

    filePath = 'data/SAR_CPLX_20190823071330_9.6G_VH_11_pres_2_fdc_246.sar.rgo.sig.nc'
    file = nc.Dataset(filePath, 'r')
    modeltransformationtag = file.variables['ModelTransformationTag'][:]
    gbpgridinfo = file.variables['GBPGridInfo'][:]
    latimage = file.variables['LatImage'][:]
    lonimage = file.variables['LonImage'][:]
    realimg = file.variables['SigmaImageSingleLookRealPart'][:]
    imagimg = file.variables['SigmaImageSingleLookImaginaryPart'][:]
    fullimg = realimg + 1j*imagimg
    demimage = file.variables['DEMImage'][:]

    

    img_adjusted = compressImg(fullimg)
    convertToGeotiff_affine(img_adjusted, modeltransformationtag, gbpgridinfo, 'compare/newslc.tif')
    img = rasterio.open('compare/newslc.tif')
    print(img.shape, latimage.shape)
    #checking 4 corners
    test_coord_accuracy('compare/newslc.tif', 0, 0, latimage, lonimage)
    test_coord_accuracy('compare/newslc.tif', 0, latimage.shape[0]-1, latimage, lonimage)
    test_coord_accuracy('compare/newslc.tif', latimage.shape[1]-1, 0, latimage, lonimage)
    test_coord_accuracy('compare/newslc.tif', latimage.shape[1]-1, latimage.shape[0]-1, latimage, lonimage)

    test_coord_accuracy('compare/newslc.tif', 100, 200, latimage, lonimage)
    test_coord_accuracy('compare/newslc.tif', 1800, 2000, latimage, lonimage)
    test_coord_accuracy('compare/newslc.tif', 2000, 1859, latimage, lonimage)
    test_coord_accuracy('compare/newslc.tif', 1256, 50, latimage, lonimage)

    validfullimg, validmodeltransformationtag = getValidImage(fullimg, modeltransformationtag)
    [row, col] = np.where(fullimg.real != -9999)   
    minrow = row.min()   
    maxrow = row.max()
    mincol = col.min()
    maxcol = col.max()
    latimage = latimage[minrow:maxrow, mincol:maxcol]
    lonimage = lonimage[minrow:maxrow, mincol:maxcol]
    
    validimg_adjusted = compressImg(validfullimg)
    convertToGeotiff_affine(validimg_adjusted, validmodeltransformationtag, gbpgridinfo, 'compare/newslccropped.tif')
    test_coord_accuracy('compare/newslccropped.tif', 0, 0, latimage, lonimage)
    test_coord_accuracy('compare/newslccropped.tif', 0, latimage.shape[0]-1, latimage, lonimage)
    test_coord_accuracy('compare/newslccropped.tif', latimage.shape[1]-1, 0, latimage, lonimage)
    test_coord_accuracy('compare/newslccropped.tif', latimage.shape[1]-1, latimage.shape[0]-1, latimage, lonimage)

    test_coord_accuracy('compare/newslccropped.tif', 100, 200, latimage, lonimage)
    test_coord_accuracy('compare/newslccropped.tif', 1800, 2000, latimage, lonimage)
    test_coord_accuracy('compare/newslccropped.tif', 2000, 1859, latimage, lonimage)
    test_coord_accuracy('compare/newslccropped.tif', 1256, 50, latimage, lonimage)
    print('\n')

    file_Xband = nc.Dataset('data/SAR_MLOOK_20181115053707_9.6G_VV_11_pres_8_fdc_-181_-91_90_180_271.sar.sig.pow.db.nc')
    file_KUband = nc.Dataset('data/SAR_MLOOK_20181115053707_17.2G_VV_11_pres_8_fdc_-272_-181_-91_90_180_271.sar.sig.pow.db.nc')

    latimage_X = file_Xband.variables['LatImage'][:]
    latimage_KU = file_KUband.variables['LatImage'][:]
    lonimage_X = file_Xband.variables['LonImage'][:]
    lonimage_KU = file_KUband.variables['LonImage'][:]
    image_X = file_Xband.variables['SigmaImageAmplitude'][:]
    # print(image_X.shape, latimage_X.shape, lonimage_X.shape)
    x_X = int(latimage_X.shape[1]/2)
    y_X = int(latimage_X.shape[0]/2)
    x_KU = int(latimage_KU.shape[1]/2)
    y_KU = int(latimage_KU.shape[0]/2)

    
    test_coord_accuracy('compare/newfullku.tif', 0,0, latimage_KU, lonimage_KU)
    test_coord_accuracy('compare/newfullx.tif', 0,0, latimage_X, lonimage_X)
    test_coord_accuracy('compare/newfullku.tif', 1019, 649, latimage_KU, lonimage_KU)
    test_coord_accuracy('compare/newfullx.tif', 1019, 649, latimage_X, lonimage_X)
    test_coord_accuracy('compare/newfullku.tif', x_KU,y_KU, latimage_KU, lonimage_KU)
    test_coord_accuracy('compare/newfullx.tif', x_X,y_X, latimage_X, lonimage_X)
