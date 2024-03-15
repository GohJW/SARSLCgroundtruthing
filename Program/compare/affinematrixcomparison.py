import netCDF4 as nc
import numpy as np
import sys, os
sys.path.append(os.path.join(sys.path[0],".."))
from SLCtogeotiff.SLCtogeotiff import convertToGeotiff_affine
from SLCtogeotiff.adjustSLC import compressImg
from PixelLatLongConversion.LatLon_pixel import getPixelfromLatLonSLC, getPixelfromGeotiff
from PixelLatLongConversion.pixel_LatLon import getLatLonFromPixelSLC, getLatLonFromGeotiff

file_Xband = nc.Dataset('data/SAR_MLOOK_20181115053707_9.6G_VV_11_pres_8_fdc_-181_-91_90_180_271.sar.sig.pow.db.nc')
file_KUband = nc.Dataset('data/SAR_MLOOK_20181115053707_17.2G_VV_11_pres_8_fdc_-272_-181_-91_90_180_271.sar.sig.pow.db.nc')

latimage_X = file_Xband.variables['LatImage'][:]
latimage_KU = file_KUband.variables['LatImage'][:]
lonimage_X = file_Xband.variables['LonImage'][:]
lonimage_KU = file_KUband.variables['LonImage'][:]
# print(latimage_X.shape, latimage_KU.shape)
matrix_X = file_Xband.variables['ModelTransformationTag'][:]
matrix_KU = file_KUband.variables['ModelTransformationTag'][:]
print(matrix_X, '\n',matrix_KU)
gbp_X = file_Xband.variables['GBPGridInfo'][:]
gbp_KU = file_KUband.variables['GBPGridInfo'][:]
amplitude_X = file_Xband.variables['SigmaImageAmplitude'][:]
amplitude_KU = file_KUband.variables['SigmaImageAmplitude'][:]
# print(matrix_X, '\n', matrix_KU)
#dimensions of the KU band and X band are not aligned

file_HH = nc.Dataset('data/SAR_CPLX_20190823071330_9.6G_HH_12_pres_2_fdc_246.sar.rgo.sig.nc')
file_VH = nc.Dataset('data/SAR_CPLX_20190823071330_9.6G_VH_11_pres_2_fdc_246.sar.rgo.sig.nc')

latimage_HH = file_HH.variables['LatImage'][:]
latimage_VH = file_VH.variables['LatImage'][:]
# print(np.any(np.not_equal(latimage_VH, latimage_HH)))
lonimage_HH = file_HH.variables['LonImage'][:]
lonimage_VH = file_VH.variables['LonImage'][:]
# print(np.any(np.not_equal(lonimage_HH, lonimage_VH)))
matrix_HH = file_HH.variables['ModelTransformationTag'][:]
matrix_VH = file_VH.variables['ModelTransformationTag'][:]
# print(np.any(np.not_equal(matrix_HH, matrix_VH)))
gbpgridinfo_HH = file_HH.variables['GBPGridInfo'][:]
#same image with same band, diff polarisation has the exact same transformation matrix given its the same dimensions



#generate_geotiff from image
img_KU = compressImg(amplitude_KU)
convertToGeotiff_affine(amplitude_KU, matrix_KU, gbp_KU, 'images/img_KU')
img_X = compressImg(amplitude_X)
convertToGeotiff_affine(amplitude_X, matrix_X, gbp_X, 'images/img_X')

#lat and lon obtained from KU geotiff image
lat, lon = getLatLonFromGeotiff('images/img_KU.tif', 3335, 2788)
print(lat,lon)
print(getLatLonFromPixelSLC(3335, 2788, matrix_KU, gbp_KU))
print(latimage_KU[2788][3335], lonimage_KU[2788][3335])

print(gbp_KU[3], gbp_KU[4])

#using this latlon, obtain the x y from KU MLK
x, y = getPixelfromLatLonSLC(lat, lon, matrix_KU, gbp_KU)
print(x, y)
#get the latlon back in the KU SLC
new_lat_KU = latimage_KU[y][x]
new_lon_KU = lonimage_KU[y][x]
print(new_lat_KU, new_lon_KU)

#use this same latlon to obtain x y from X MLK
x, y = getPixelfromLatLonSLC(lat,lon, matrix_X, gbp_X)
print(x,y)
#get the latlon back in the X SLC
# new_lat_X = latimage_X[y][x]
# new_lon_X = lonimage_X[y][x]
# print(new_lat_X, new_lon_X)

realimg_HH = file_HH.variables['SigmaImageSingleLookRealPart'][:]
imagimg_HH = file_HH.variables['SigmaImageSingleLookImaginaryPart'][:]
img_HH = compressImg(realimg_HH + 1j*imagimg_HH)
convertToGeotiff_affine(img_HH, matrix_HH, gbpgridinfo_HH, 'images/img_HH')	

lat, lon = getLatLonFromGeotiff('images/img_HH.tif', 1019, 649)
print(lat,lon)

x,y = getPixelfromLatLonSLC(lat, lon, matrix_HH, gbpgridinfo_HH)
new_lat_HH = latimage_HH[y][x]
new_lon_HH = lonimage_HH[y][x]
print(new_lat_HH, new_lon_HH)


