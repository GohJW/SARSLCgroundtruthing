import netCDF4 as nc
from utm import to_latlon
import rasterio
from pyproj import Transformer

def getLatLonFromPixelSLC(x, y, ModelTransformationTag, GBPGridInfo):
    """returns the lat and lon of each pixel in the image. Assumes y is row, x is column.

    Args:
        x: x-coord of image
        y: y-coord of image
        ModelTransformationTag: transformation matrix
        GBPGridInfo: array containing grid information
        
    Returns:
        Latitude: The latitude
        Longitude: The longitude
    """
    negativeDxcosrot = ModelTransformationTag[0][0]
    Dysinrot = ModelTransformationTag[1][0]
    Dxsinrot = ModelTransformationTag[0][1]
    Dycosrot = ModelTransformationTag[1][1]
    x_offset = ModelTransformationTag[3][0]
    y_offset = ModelTransformationTag[3][1]
    
    zone = GBPGridInfo[3]
    south = GBPGridInfo[4]
    
    UTM_x = negativeDxcosrot*x + Dysinrot*y + x_offset
    UTM_y = Dxsinrot*x + Dycosrot*y + y_offset
    lat, long = to_latlon(UTM_x, UTM_y, zone, 'S' if south == 1 else 'N')
    return lat, long

def getLatLonFromGeotiff(imgPath, x, y):
    """Calculates and returns the Longitude and Latitude coordinates of a pixel(x,y) in image. Assumes y is row, x is column.

    Args:
        imgPath: the path to the source image
        x: x coordinate
        y: y coordinate

    Returns:
        Lat, Lon: the Longitude and Latitude if x and y are within the bounds
    """
    img = rasterio.open(imgPath)
    
    #check if x and y pixel coordinate is within the bounds of the img raster
    max_x = img.width - 1
    max_y = img.height - 1
    if(x < 0 or x > max_x or y < 0 or y > max_y):
        print(f"x and y coords {x,y} are out of bounds for the given geotiff")
        return

    # #obtain geotransform information
    img_gt = img.transform
    
    #convert geotransform data to long and lat
    lon = img_gt[0] * x + img_gt[1] * y + img_gt[2]
    lat = img_gt[3] * x + img_gt[4] * y + img_gt[5]
    
    #transform current obtained long and lat into EPSG:4326 coordinate system
    transformer = Transformer.from_crs(img.crs,"EPSG:4326",always_xy=True) 
    lon, lat = transformer.transform(lon,lat)
    return lat, lon
