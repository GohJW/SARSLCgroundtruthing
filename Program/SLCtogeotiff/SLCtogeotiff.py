import rasterio as rio
from pyproj import CRS
from rasterio.transform import from_gcps, Affine
from rasterio.control import GroundControlPoint
import os

def convertToGeotiff_gcp(img, latimage, lonimage, outputpath):
    """Converts the image into a geotiff at the outputpath using ground control points

    Args:
        img: The SLC image
        latimage: The matrix containing latitude for each pixel
        lonimage: The matrix containing longitude for each pixel
        outputpath: the output path for the generated geotiff image
        
    """
    max_x, max_y = img.shape[1], img.shape[0]
    tl = GroundControlPoint(0,0, lonimage[0][0], latimage[0][0],0)
    tr = GroundControlPoint(max_y-1, 0, lonimage[-1][0], latimage[-1][0], 0)
    br = GroundControlPoint(max_y-1, max_x-1, lonimage[-1][-1], latimage[-1][-1], 0)
    bl = GroundControlPoint(0, max_x-1, lonimage[0][-1], latimage[0][-1], 0)
    gcp_vectors = [br, bl, tl, tr]
    transform = from_gcps(gcp_vectors)
    if os.path.exists(outputpath):
        os.remove(outputpath)

    with rio.open(outputpath, 'w', crs = "EPSG:4326", driver= 'GTiff', height= img.shape[0], width= img.shape[1], count= 1, dtype=img.dtype, transform=transform, nodata = -1) as dst:
        dst.write(img, 1)
   
   #Using Affine Transformation matrix
def convertToGeotiff_affine(img, modeltransformationtag, gbpgridinfo, outputpath):
    """Converts the image into a geotiff at the outputpath using affine transformation matrix at the specified minimum and maximum rows and columns.

    Args:
        img: The SLC image
        modeltransformationtag: the modeltransformationtag of the given image
        gbpgridinfo: the gbpgridinfo of the given image
        outputpath: the output path of the geotiff
    """
    crs = CRS.from_dict({'proj': 'utm', 'zone': int(gbpgridinfo[3]), 'south': bool(gbpgridinfo[4])})
    transform = Affine(a = modeltransformationtag[0][0], b = modeltransformationtag[1][0],c = modeltransformationtag[3][0],d = modeltransformationtag[0][1], e = modeltransformationtag[1][1], f = modeltransformationtag[3][1])
        
    with rio.open(outputpath, 'w', crs = crs, driver= 'GTiff', height= img.shape[0], width= img.shape[1], count= 1, dtype=img.dtype, transform=transform, nodata=-1) as dst:
        dst.write(img, 1)   
