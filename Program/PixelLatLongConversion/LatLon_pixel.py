from utm import from_latlon
import sympy as sym
from pyproj import Transformer
import rasterio


def getPixelfromLatLonSLC(lat, lon, ModelTransformationTag, GBPGridInfo):
    """returns the UTM conversion of each pixel in the image

    Args:
        latitude
        longitude
        zone: UTM zone of the image
        south: UTM hemisphere
        ModelTransformationTag: transformation matrix
        GBPGridInfo: array containing grid information
        
    Returns:
        x: the x-coord in the SLC image
        y: the y-coord in the SLC image 
    """
    negativeDxcosrot = ModelTransformationTag[0][0]
    Dysinrot = ModelTransformationTag[1][0]
    Dxsinrot = ModelTransformationTag[0][1]
    Dycosrot = ModelTransformationTag[1][1]
    x_offset = ModelTransformationTag[3][0]
    y_offset = ModelTransformationTag[3][1]
    
    zone = GBPGridInfo[3]
    south = GBPGridInfo[4]
    
    UTM = from_latlon(lat, lon, zone, 'S' if south == 1 else 'N')
    UTM_x = UTM[0]
    UTM_y = UTM[1]
    x,y = sym.symbols('x,y')
    eq1 = sym.Eq(negativeDxcosrot*x + Dysinrot*y + x_offset, UTM_x)
    eq2 = sym.Eq(Dxsinrot*x + Dycosrot*y + y_offset, UTM_y)
    result = sym.solve([eq1,eq2],(x,y))
    x = round(result[x])
    y = round(result[y])
    return x, y
    
def getPixelfromGeotiff(lat, lon, geotiffPath):
    """returns the row and column coordinates on the geotiff of the given latlon coordinates.
    Args:
        lat: Latitude
        lon: Longitude
        geotiffPath: The path to the given geotiff image
        
    Returns:
        r: row-coordinate in the geotiff
        c: column-coordinate in the geotiff
    """
    geotiff = rasterio.open(geotiffPath)
    img_crs = geotiff.crs
    img_gt = geotiff.transform 
    #transform long and lat into image coordinate system
    transformer = Transformer.from_crs("EPSG:4326", img_crs, always_xy=True) 
    UTM_x, UTM_y = transformer.transform(lon,lat)
    x,y = sym.symbols('x,y')
    eq1 = sym.Eq(img_gt[0] * x + img_gt[1] * y + img_gt[2], UTM_x)
    eq2 = sym.Eq(img_gt[3] * x + img_gt[4] * y + img_gt[5], UTM_y)
    result = sym.solve([eq1,eq2], (x,y))
    c = round(result[x])
    r = round(result[y])
    return r,c
