import netCDF4 as nc
import numpy as np
from math import asin, sin, pi, cos, atan2, sqrt, radians, degrees, acos
import rasterio as rio
from rasterio.windows import Window
import os, csv, shutil, sys, yaml


import sys
sys.path.append(os.path.join(sys.path[0],".."))
from PixelLatLongConversion.LatLon_pixel import getPixelfromGeotiff, getPixelfromLatLonSLC
from anglecalculations.anglecalculation import getGrazingAngle, getIntegrationAngle
from iceyeutils.utils import dr_coherent_multiple_inwtoutwt_preserveintensity
from degradeimage.slcdegrade import degradeSLC

        

def findCornerLatLonGivenDistance(lat1, lon1, Ad, b):
    """Find the latlon coordinates of the destination given the source latlon, angular distance and bearing

    Args:
        lat1: source latitude
        lon1: source longitude
        Ad: angular distance calculated from distance to destination/radius of Earth
        b: bearing of the destination with respect to the source
        
    Returns:
        lat2: destination latitude
        lon2: destination longitude
    """
    lat1_rad = radians(lat1)
    lon1_rad = radians(lon1)
    b_rad = radians(b)
    lat2_rad = asin(sin(lat1_rad) * cos(Ad) + cos(lat1_rad) * sin(Ad) * cos(b_rad))
    lon2_rad = lon1_rad + atan2(sin(b_rad) * sin(Ad) * cos(lat1_rad), cos(Ad) - sin(lat1_rad) * sin(lat2_rad))
    
    lat2 = degrees(lat2_rad)
    lon2 = degrees(lon2_rad)
    return lat2, lon2

def cropGeotiff(geotiffPath, col_tl, col_tr, row_tl, row_bl, outputPath):
    """Crops the geotiff at the given row and columns. Note outputPath does not need to include \'.tif\'.

    Args:
        geotiffPath: path to the original full sized geotiff
        col_tl: the top left column of the window
        col_tr: the top right column of the window
        row_tl: the top left row of the window
        row_bl: the bottom left row of the window 
        outputPath: the path to save the cropped geotiff 
    """
    with rio.open(geotiffPath) as src:
        #size of the window
        csize = abs(col_tl - col_tr)
        rsize = abs(row_tl - row_bl)
        window = Window(col_tl, row_tl, csize, rsize)
        transform = src.window_transform(window)
        window_profile = src.profile
        window_profile.update({
            'height' : rsize,
            'width' : csize,
            'transform' : transform
        })
        with rio.open(outputPath + '.tif', 'w', **window_profile) as dst:
            dst.write(src.read(window = window))

# def cropGeotiff(geotiffPath, x_tl, x_tr, y_tl, y_bl, outputPath):
    # with rio.open(geotiffPath) as src:
    #     #size of the window
    #     xsize = abs(x_tl - x_tr)
    #     ysize = abs(y_tl - y_bl)
        
    #     #prevent cropping of pixels beyond image dimensions
    #     xmax =  src.width - xsize 
    #     ymax =  src.height - ysize
    #     if(x_tl > xmax):
    #         print("x-coordinate of cropped geotiff overflows from the edge of the image, cropped image may be smaller than actual size specified.")
    #     if(y_tl > ymax):
    #         print("y-coordinate of cropped geotiff overflows from the edge of the image, cropped image may be smaller than actual size specified.")
    #     window = Window(x_tl, y_tl, xsize, ysize)
    #     transform = src.window_transform(window)
    #     window_profile = src.profile
    #     window_profile.update({
    #         'height' : xsize,
    #         'width' : ysize,
    #         'transform' : transform
    #     })
    #     with rio.open(outputPath + '/croppedgeotiff.tif', 'w', **window_profile) as dst:
    #         dst.write(src.read(window = window))

def getDist(lat1, lon1, lat2, lon2):
    """Gets the approximate distance between 2 latlon coordinates

    Args:
        lat1: latitude of point 1
        lon1: longitude of point 1
        lat2: latitude of point 2
        lon2: longitude of point 2

    Returns:
        distance: the approximate distance between the two latlon points
    """
    lat1_rad = radians(lat1)
    lon1_rad = radians(lon1)
    lat2_rad = radians(lat2)
    lon2_rad = radians(lon2)
    dotproduct = cos(lat1_rad) * cos(lat2_rad) * cos(lon1_rad-lon2_rad) + sin(lat1_rad) * sin(lat2_rad)
    angle = acos(dotproduct)
    distance = 6378 * angle #radius of earth * angle
    return distance

# def cropSLC(slcimage, lat_tl, lon_tl, lat_br, lon_br, modeltransformationtag, gbpgridinfo, outputPath):
    """Crops out a SLCimage and saves as a npy file

    Args:
        slcimage: the SLC numpy array
        lat_tl: latitude of the topleft
        lon_tl: longitude of the topleft
        lat_br: latitude of the bottom right
        lon_br: longitude of the bottom right
        modeltransformationtag: the modeltransformationtag of the SAR image
        gbpgridinfo: the gbpgridinfo of the SAR image
        outputPath: the destination path of the cropped SLC (no need to append filetype at the end)
    """
    x_tl, y_tl = getPixelfromLatLonSLC(lat_tl, lon_tl,modeltransformationtag, gbpgridinfo)
    x_br, y_br = getPixelfromLatLonSLC(lat_br, lon_br, modeltransformationtag, gbpgridinfo)
    # print(x_tl, y_tl, x_br, y_br)
    x_tl, x_br = min(x_tl, x_br), max(x_tl, x_br)
    y_tl, y_br = min(y_br, y_tl), max(y_tl, y_br)
    # print(x_tl, y_tl, x_br, y_br)
    croppedslc = slcimage[y_tl:y_br, x_tl:x_br]
    np.save(outputPath + '/croppedslc', croppedslc.data)

def getWindowBoundaries(slcimage, row_tl, col_tl, rowpixsize, colpixsize, chipsize):
    """Obtains the boundary rows and columns of the full sized slcimage to cut the chip.

    Args:
        slcimage: the original slcimage
        row_tl: the top left row of the window
        col_tl: the top left column of the window
        rowpixsize: the row pixel size of the slc image
        colpixsize: the column pixel size of the slc image 
        chipsize: the chip size in metres

    Returns:
        row_tl, row_bl,col_tl, col_tr : the boundaries of the window
    """
    rowpix = round(chipsize/rowpixsize)
    colpix = round(chipsize/colpixsize)

    maxrow, maxcol = slcimage.shape
    print('SLCimage shape:',maxrow, maxcol)
    if row_tl + rowpix > maxrow:
        print("row of cropped SLC overflows from the edge of the array, cropped SLC may be smaller than actual size specified.")
    if col_tl + colpix > maxcol:
        print("column of cropped SLC overflows from the edge of the array, cropped SLC may be smaller than actual size specified.")
    row_bl = min(maxrow, row_tl + rowpix)
    col_tr = min(maxcol, col_tl + colpix)
    print('bottommost window row:',row_tl + rowpix)
    print('rightmost window column:',col_tl + colpix)
    print('window boundaries by row and columns: ',row_tl, row_bl,col_tl, col_tr)

    return row_tl, row_bl,col_tl, col_tr
  
def cropSLC(slcimage, row_tl, row_bl, col_tl, col_tr):
    """Crops the SLC image at the given row and columns.

    Args:
        slcimage: the original SLC image
        row_tl: the top left row of the window
        row_bl: the bottom left row of the window 
        col_tl: the top left column of the window
        col_tr: the top right column of the window

    Returns:
        croppedslc: returns the cropped SLC image array
    """
    croppedslc = slcimage[row_tl:row_bl, col_tl:col_tr]
    return croppedslc

def saveSLC(slcimage, outputPath):
    np.save(outputPath, slcimage)
# def cutChip(lat, lon, size, geotiffPath, slcimage, gbpgridinfo, orbitlatitude, orbitlongitude, orbitheight, colpixsize, rowpixsize, transpose, outputPath):
#     if os.path.exists(outputPath):
#         shutil.rmtree(outputPath)
#     os.makedirs(outputPath)
    
    
#     R = 6371 #radius of earth
#     size_km = size/1000
#     distance = sqrt((size_km/2)**2 + (size_km/2)**2)
#     Ad = distance/R #angular distance

#     b_tl = 315 # tl bearing is north west of center
#     lat_tl, lon_tl = findCornerLatLonGivenDistance(lat, lon, Ad, b_tl)
#     # print(f'top left latlon is {lat_tl, lon_tl}')

#     # b_tr = 45 # tr bearing is north east of center
#     # lat_tr, lon_tr = findCornerLatLonGivenDistance(lat, lon, Ad, b_tr)
#     # print(f'top right latlon is {lat_tr, lon_tr}')
    
#     # b_bl = 225 # bl bearing is south west of center
#     # lat_bl, lon_bl = findCornerLatLonGivenDistance(lat, lon, Ad, b_bl)
#     # print(f'bottom left latlon is {lat_bl, lon_bl}')
    
#     # b_br = 135 # br bearing is south east of center
#     # lat_br, lon_br = findCornerLatLonGivenDistance(lat, lon, Ad, b_br)
#     # print(f'bottom right bearing is {lat_br, lon_br}')
#     # getDist(lat_tl, lon_tl, lat_bl, lon_bl)
#     # getDist(lat_br, lon_br, lat_bl, lon_bl)
#     # getDist(lat_br, lon_br, lat_tr, lon_tr)
#     # getDist(lat_tl, lon_tl, lat_tr, lon_tr)

#     row_tl, col_tl = getPixelfromGeotiff(lat_tl, lon_tl, geotiffPath)
#     if transpose:
#         row_tl, col_tl = col_tl, row_tl
#     col_tl, col_tr,row_tl, row_bl = getWindowBoundaries(slcimage, row_tl, col_tl, rowpixsize, colpixsize, size, transpose)
#     print(col_tl, col_tr,row_tl, row_bl)
#     cropSLC(slcimage,col_tl, col_tr,row_tl, row_bl, outputPath)
#     integrationangle = getIntegrationAngle(orbitlatitude, orbitlongitude, orbitheight, lat, lon, gbpgridinfo)
#     cropGeotiff(geotiffPath, col_tl, col_tr,row_tl, row_bl, outputPath)
#     with open(outputPath + '/chipinfo.csv', 'w', newline = '') as file:
#         writer = csv.writer(file)
#         writer.writerow(['center coordinates', lat, lon])
#         writer.writerow(['IntegrationAngle', integrationangle])
            
# def getChips(lat, lon, size, geotiffPathlist, slcimagelist, gbpgridinfolist, orbitlatitudelist, orbitlongitudelist, orbitheightlist,groundrangepixsizelist, crossrangepixsizelist, transposelist, outputPath):
#     if os.path.exists(outputPath):
#         shutil.rmtree(outputPath)
#     os.makedirs(outputPath)
#     for index, (geotiffPath, slcimage, gbpgridinfo, orbitlatitude, orbitlongitude, orbitheight, groundrangepixsize, crossrangepixsize, transpose) in \
#     enumerate(zip(geotiffPathlist, slcimagelist, gbpgridinfolist, orbitlatitudelist, orbitlongitudelist, orbitheightlist,groundrangepixsizelist, crossrangepixsizelist, transposelist)):
#         print(index)
#         cutChip(lat, lon, size, geotiffPath, slcimage, gbpgridinfo, orbitlatitude, orbitlongitude, orbitheight,groundrangepixsize, crossrangepixsize, transpose, outputPath, index)
        
    