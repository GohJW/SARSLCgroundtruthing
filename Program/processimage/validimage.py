import numpy as np
import os, sys
sys.path.append(os.path.join(sys.path[0],".."))
import copy

def recalculateAffineTransformationMatrix(modeltransformationtag, minrow, mincol):
    """Adjusts the modeltransformation tag to fit the new image based on its minimum row and column.

    Args:
        modeltransformationtag: the modeltransformationtag of the original image
        minrow: the minimum row specified
        mincol: the minimum column specified

    Returns:
        modeltransformationtag_new: the new adjusted modeltransformationtag to use in creating the geotiff of the cropped SLC image
    """
    modeltransformationtag_new = copy.deepcopy(modeltransformationtag)
    modeltransformationtag_new[3][0] = modeltransformationtag[3][0] + (mincol * modeltransformationtag[0][0]) + (minrow * modeltransformationtag[1][0])   
    modeltransformationtag_new[3][1] = modeltransformationtag[3][1] + (mincol * modeltransformationtag[0][1]) + (minrow * modeltransformationtag[1][1])
    # print(modeltransformationtag)
    # print(modeltransformationtag_new)
    return modeltransformationtag_new

def getValidBoundaries(fullimg):
    [row, col] = np.where(fullimg.real != -9999)
    minrow = row.min()
    maxrow = row.max()
    mincol = col.min()
    maxcol = col.max()
    return minrow, maxrow, mincol, maxcol

def getValidImage(fullimg, modeltransformationtag):
    """Removes the border of the SLC where there is dummy data. Note this function chooses the minimum and maximum rows and columns where there is a valid data pixel, uneven
    images with jagged borders would still contain dummy data that needs to be filtered out in the geotiff generation as nodata.

    Args:
        fullimg: the raw SLC image with borders
        modeltransformationtag: the original modeltransformationtag of the raw SLC image

    Returns:
        validimg: the cropped SLC with rows and columns of dummy data removed
        validmodeltransformationtag: the updated modeltransformationtag for the cropped SLC image dimensions
    """
    #Crop the image based on the min and max columns that contain data
    minrow, maxrow, mincol, maxcol = getValidBoundaries(fullimg)
    validimg = fullimg[minrow:maxrow+1, mincol:maxcol+1]
    validmodeltransformationtag = recalculateAffineTransformationMatrix(modeltransformationtag, minrow, mincol)
    return validimg, validmodeltransformationtag