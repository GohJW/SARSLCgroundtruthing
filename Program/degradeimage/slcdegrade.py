import netCDF4 as nc
import numpy as np
import os, sys
import matplotlib.pyplot as plt
from math import radians, sin

sys.path.append(os.path.join(sys.path[0],".."))
from iceyeutils.utils import dr_coherent_multiple_inwtoutwt_preserveintensity
from anglecalculations.anglecalculation import getIntegrationAngle, getGrazingAngle
from SLCtogeotiff.adjustSLC import compressImg
from SLCtogeotiff.SLCtogeotiff import convertToGeotiff_affine
from processimage.validimage import getValidBoundaries, getValidImage


def adjustPhase(SLCchip, latimage, lonimage, orblatimage, orblonimage, orbheightimage, wavelength):
    """Adjusts the phase of the SLCchip before degradation is done. Uses 2D perpendicular distance. Note all images need to have the same dimensions as the SLCchip.

    Args:
        SLCchip: the SLC image of the chip
        latimage: the latitude image of the chip
        lonimage: the longitude image of the chip
        orblatimage: the latitude image of the satellite with respect to the chip
        orblonimage: the longitude image of the satellite with respect to the chip
        orbheightimage: the height image of the satellite with respect to the chip
        wavelength: the wavelength of the transmitted signal

    Returns:
        img_comp2: the phase adjusted SLC chip
    """

    r = 6378137.0
    alt = 0 #Assumes height of pixel is 0
    orbalt = np.mean(orbheightimage)
    r1 = r + alt
    r2 = r + orbalt
    
    y1 = r1 * np.sin(lonimage) * np.sin(latimage)
    z1 = r1 * np.cos(lonimage)

    y2 = r2 * np.sin(orblonimage) * np.sin(orblatimage)
    z2 = r2 * np.cos(orblonimage)

    R_perp = np.linalg.norm([y1 - y2, z1 - z2], axis = 0)
    wl = np.array(wavelength)
    phase_comp2 = np.exp(-1j * 4 * np.pi / wl * R_perp)
    img_comp2 = SLCchip*phase_comp2
    return img_comp2

def degradeSLC(SLCchip, transpose, azimuthres, rangeres, groundrangepixsize, crossrangepixsize,oprangeres,oprangepixsize,opazimuthres,opazimuthpixsize,win_typein,win_typeout,nofftshift = False,norm_fac = False,numsubaperture = 1,numsubband = 1):
    """Degrades and returns the SLC chip. Since the degradation function assumes range is row and azimuth is column, set transpose to True in order for
    range and azimuth to be swapped. Note that the desired output resolution and pixelsize might differ from the actual output resolution and pixelsize.

    Args:
        SLCchip: the SLC image of the chip
        transpose: boolean value to transpose the image before degradation
        azimuthres: azimuth resolution of the chip
        rangeres: range resolution of the chip
        groundrangepixsize: ground range pixel size of the chip 
        crossrangepixsize: cross range pixel size of the chip 
        oprangeres: the desired output range resolution
        oprangepixsize: the desired output range pixel size
        opazimuthres: the desired output azimuth resolution 
        opazimuthpixsize: the desired output azimuth pixel size 
        win_typein: the window function of the input image 
        win_typeout: the desired window function of the output image 
        nofftshift: boolean value to fftshift the image. Defaults to False.
        norm_fac: boolean value for the norm factor. Defaults to False.
        numsubaperture: the number of subapertures. Defaults to 1.
        numsubband: the number of subbands. Defaults to 1.

    Returns:
        img_outlist[0][0],oprowres,oprowpixsize,opcolres,opcolpixsize,norm_fac_value
    """
    if transpose:
        ipcolres = azimuthres
        ipcolpixsize = crossrangepixsize
        iprowres = rangeres
        iprowpixsize = groundrangepixsize
        opcolres = opazimuthres
        opcolpixsize = opazimuthres
        oprowres = oprangeres
        oprowpixsize = oprangepixsize

    else:
        ipcolres = rangeres
        ipcolpixsize = groundrangepixsize
        iprowres = azimuthres
        iprowpixsize = crossrangepixsize
        opcolres = oprangeres
        opcolpixsize = oprangepixsize
        oprowres = opazimuthres
        oprowpixsize = opazimuthres
    img_outlist,oprowres,oprowpixsize,opcolres,opcolpixsize ,norm_fac_value = dr_coherent_multiple_inwtoutwt_preserveintensity(SLCchip,iprowres,iprowpixsize,ipcolres,ipcolpixsize,oprowres,oprowpixsize,opazimuthres,opazimuthpixsize,win_typein,win_typeout,nofftshift,norm_fac,numsubaperture,numsubband)
    return img_outlist[0][0],oprowres,oprowpixsize,opcolres,opcolpixsize ,norm_fac_value

if __name__ == '__main__':
    file = nc.Dataset('data/SAR_CPLX_20190823071330_9.6G_HH_12_pres_2_fdc_246.sar.rgo.sig.nc')
    orbitlatitude = file.variables['OrbitLatitude'][:]
    orbitlongitude = file.variables['OrbitLongitude'][:]
    orbitheight = file.variables['OrbitHeight'][:]
    
    orblatimage = file.variables['OrbLatImage'][:]
    orblonimage = file.variables['OrbLonImage'][:]
    orbheightimage = file.variables['OrbHeightImage'][:]
    
    modeltransformationtag = file.variables['ModelTransformationTag'][:]
    gbpgridinfo = file.variables['GBPGridInfo'][:]
    centralfreq = file.variables['CentralFreq'][:]
    
    transmittedpolarization = str(file.variables['TxPolarization'][:], 'utf-8')
    receivedpolarization = str(file.variables['RxPolarization'][:], 'utf-8')
    polarization = transmittedpolarization + receivedpolarization
    transmittedbandwidth = file.variables['TransmittedBandWidth'][:]
    c = 299792458 # Speed of light

    latimage = file.variables['LatImage'][:]
    lonimage = file.variables['LonImage'][:]
    
    realimg = file.variables['SigmaImageSingleLookRealPart'][:]
    imagimg = file.variables['SigmaImageSingleLookImaginaryPart'][:]
    fullimg = realimg + 1j*imagimg
    window = file.variables['WindowFunction'][:]
    match window:
        case 0:
            window = 'rect' 
            k = 0.89
        case 1:
            window = 'tri'
        case 2:
            window = 'hamming'
            k = 1.3
        case 3:
            window = 'hanning'
            k = 1.4466
    
    window = 'rect' ##hardcode
    k = 0.89
    
    minrow, maxrow, mincol, maxcol = getValidBoundaries(fullimg)
    validimage, validmpdeltransformationtag = getValidImage(fullimg, modeltransformationtag)
    validimage = fullimg[minrow:maxrow+1, mincol:maxcol-2]
    compressedimg = compressImg(validimage)
    convertToGeotiff_affine(compressedimg, validmpdeltransformationtag, gbpgridinfo, 'test.tif')
    
    #calculate azimuth resolution and range resolution 
    row, col = validimage.shape
    #find centre of full image
    row = int(row/2)
    col = int(col/2)
    teta = getIntegrationAngle(orbitlatitude, orbitlongitude, orbitheight, latimage[row][col], lonimage[row][col], gbpgridinfo)  
    azimuthres = k * (c / (2*transmittedbandwidth))
    print(k, transmittedbandwidth, azimuthres)
    wavelength = c/centralfreq
    rangeres = k * ( wavelength/ (2 * radians(teta))) 
    
    #calculate groudn range and cross range pixel size        
    groundrange = file.variables['GroundRange'][:]
    groundrangepixsize = abs(groundrange[0] - groundrange[-1])/ len(groundrange)
    crossrange = file.variables['CrossRange'][:]
    crossrangepixsize = abs(crossrange[0] - crossrange[-1]) / len(crossrange)
        
    validimage = fullimg[minrow:maxrow+1, mincol:maxcol-2]
    adjustedimage = adjustPhase(validimage, latimage[minrow:maxrow+1, mincol:maxcol-2], lonimage[minrow:maxrow+1, mincol:maxcol-2], orblatimage[minrow:maxrow+1, mincol:maxcol-2], orblonimage[minrow:maxrow+1, mincol:maxcol-2], orbheightimage[minrow:maxrow+1, mincol:maxcol-2], wavelength)
    img_outlist,oprowres,oprowpixsize,opcolres,opcolpixsize ,norm_fac_value = degradeSLC(adjustedimage, True, azimuthres, rangeres, groundrangepixsize, crossrangepixsize, 0.8, 0.25, 0.8, 0.25, window, 'hanning')
    compressedadjustedimg = compressImg(img_outlist[0][0])
    convertToGeotiff_affine(compressedadjustedimg, validmpdeltransformationtag, gbpgridinfo, 'testdegradedrect.tif')




