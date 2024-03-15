import netCDF4 as nc
import numpy as np
import os, sys

sys.path.append(os.path.join(sys.path[0],".."))
from SLCtogeotiff.SLCtogeotiff import convertToGeotiff_affine
from SLCtogeotiff.adjustSLC import compressImg
from processimage.validimage import getValidImage
from chips.cutchip import getChips

# xband = nc.Dataset('data\SAR_MLOOK_20181115053707_9.6G_VV_11_pres_8_fdc_-181_-91_90_180_271.sar.sig.pow.db.nc')
# kuband = nc.Dataset('data\SAR_MLOOK_20181115053707_17.2G_VV_11_pres_8_fdc_-272_-181_-91_90_180_271.sar.sig.pow.db.nc')

# ximg = xband['SigmaImageAmplitude'][:]
# xmodel = xband['ModelTransformationTag'][:]
# xgbp = xband['GBPGridInfo'][:]
# xorbitlat = xband['OrbitLatitude'][:]
# xorbitlon = xband['OrbitLongitude'][:]
# xorbitheight = xband['OrbitHeight'][:]


# kuimg = kuband['SigmaImageAmplitude'][:]
# kumodel = kuband['ModelTransformationTag'][:]
# kugbp = kuband['GBPGridInfo'][:]
# kuorbitlat = kuband['OrbitLatitude'][:]
# kuorbitlon = kuband['OrbitLongitude'][:]
# kuorbitheight = kuband['OrbitHeight'][:]

# validximg, validxmodel = getValidImage(ximg, xmodel)
# validkuimg, validkumodel = getValidImage(kuimg, kumodel)

# outputpath = 'chips/multipleimagetest2/' 

# compressedximg = compressImg(validximg)
# compressedkuimg = compressImg(validkuimg)
# xgeotiffpath = os.path.join(outputpath, 'ximage.tif')
# kugeotiffpath = os.path.join(outputpath, 'kuimage.tif')
# convertToGeotiff_affine(compressedximg,validxmodel,xgbp, xgeotiffpath)
# convertToGeotiff_affine(compressedkuimg,validkumodel,kugbp, kugeotiffpath)

# geotiffpathlist = [xgeotiffpath, kugeotiffpath]
# slclist = [compressedximg, compressedkuimg]
# modellist = [validxmodel, validkumodel]
# gbplist = [xgbp, kugbp]
# orbitlatitudelist = [xorbitlat, kuorbitlat]
# orbitlongitudelist = [xorbitlon, kuorbitlon]
# orbitheightlist = [xorbitheight, kuorbitheight]

# getChips(37.4079280,126.7344027, 20, geotiffpathlist, slclist, modellist, gbplist, orbitlatitudelist, orbitlongitudelist, orbitheightlist, os.path.join(outputpath, 'chip1'))
# getChips(37.4089280,126.7346027, 20, geotiffpathlist, slclist, modellist, gbplist, orbitlatitudelist, orbitlongitudelist, orbitheightlist, os.path.join(outputpath, 'chip2'))


hh = nc.Dataset('data\SAR_CPLX_20190823071330_9.6G_HH_12_pres_2_fdc_246.sar.rgo.sig.nc')
vh = nc.Dataset('data\SAR_CPLX_20190823071330_9.6G_VH_11_pres_2_fdc_246.sar.rgo.sig.nc')
outputpath = 'chips/multipleimagetest/' 


hhreal = hh['SigmaImageSingleLookRealPart'][:]
hhimag = hh['SigmaImageSingleLookImaginaryPart'][:]
hhgbp = hh['GBPGridInfo'][:]
hhmodel = hh['ModelTransformationTag'][:]
hhorbitlat = hh['OrbitLatitude'][:]
hhorbitlon = hh['OrbitLongitude'][:]
hhorbitheight = hh['OrbitHeight'][:]
hhfullimg = hhreal + 1j*hhimag
validhhimg, validhhmodel = getValidImage(hhfullimg, hhmodel)
compressedhhimg = compressImg(validhhimg)
hhimagepath = os.path.join(outputpath, 'hhimage.tif')
# convertToGeotiff_affine(compressedhhimg, validhhmodel, hhgbp, hhimagepath)  


vhreal = vh['SigmaImageSingleLookRealPart'][:]
vhimag = vh['SigmaImageSingleLookImaginaryPart'][:]
vhgbp = vh['GBPGridInfo'][:]
vhmodel = vh['ModelTransformationTag'][:]
vhorbitlat = vh['OrbitLatitude'][:]
vhorbitlon = vh['OrbitLongitude'][:]
vhorbitheight = vh['OrbitHeight'][:]
vhfullimg = vhreal + 1j*vhimag
validvhimg, validvhmodel = getValidImage(vhfullimg, vhmodel)
compressedvhimg = compressImg(validvhimg)
vhimagepath = os.path.join(outputpath, 'vhimage.tif')
# convertToGeotiff_affine(compressedvhimg, validvhmodel, vhgbp, vhimagepath)

geotiffpathlist = [hhimagepath, vhimagepath]
slclist = [compressedhhimg, compressedvhimg]
modellist = [validhhmodel, validvhmodel]
gbplist = [hhgbp, vhgbp]
orbitlatitudelist = [hhorbitlat, vhorbitlat]
orbitlongitudelist = [hhorbitlon, vhorbitlon]
orbitheightlist = [hhorbitheight, vhorbitheight]

# getChips(51.8694067,4.4809171, 20, geotiffPathlist=geotiffpathlist, slcimagelist=slclist, modeltransformationtaglist=modellist, gbpgridinfolist=gbplist, outputPath=os.path.join(outputpath, 'chip1'), \
#     orbitlatitudelist=orbitlatitudelist, orbitlongitudelist=orbitlongitudelist, orbitheightlist=orbitheightlist)

getChips(51.8696047,4.4809621, 20, geotiffPathlist=geotiffpathlist, slcimagelist=slclist, modeltransformationtaglist=modellist, gbpgridinfolist=gbplist, outputPath=os.path.join(outputpath, 'chip2'), \
    orbitlatitudelist=orbitlatitudelist, orbitlongitudelist=orbitlongitudelist, orbitheightlist=orbitheightlist)

  
