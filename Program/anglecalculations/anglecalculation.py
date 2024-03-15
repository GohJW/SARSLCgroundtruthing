import numpy as np
from pyproj import Geod
from math import degrees, cos, sin, pi, atan2, radians,acos, sqrt
from utm import from_latlon

import sys, os
sys.path.append(os.path.join(sys.path[0],".."))

def getAspectAngle(lat1, lon1, lat2, lon2):
    """Obtain the aspect angle of the object with respect to True North from lat1, lon1 to lat2, lon2. Value changes depending on order of latlon input.

    Args:
        lat1: Latitude of point 1
        lon1: Longitude of point 1
        lat2: Latitude of point 2
        lon2: Longitude of point 2
        
    Returns:
        aspect angle: the angle from point 1 to point 2 with respect to the True North
    """
    lat1_rad = radians(lat1)
    lon1_rad = radians(lon1)
    lat2_rad = radians(lat2)
    lon2_rad = radians(lon2)
    x = cos(lat2_rad) * sin(lon2_rad - lon1_rad)
    y = cos(lat1_rad) * sin(lat2_rad) - sin(lat1_rad) * cos(lat2_rad) * cos(lon2_rad-lon1_rad)
    aspectangle_rad = atan2(x, y)
    aspectangle_deg = degrees(aspectangle_rad)
    return aspectangle_deg

def getAspectAngle_r(heading, aspectangle_n, lookdirection):
    """Obtain the aspect angle of the object with respect to Line of Sight.

    Args:
        heading: heading of the aircraft
        aspectangle_n: aspect angle of the target relative to true north
        lookdirection: lookdirection of the aircraft
        
    Returns:
        aspectangle_r: aspect angle of the object with respect to Line of Sight.
    """
    lineofsight = heading + 90 if lookdirection == 'R' else heading - 90
    if lineofsight > 360:
        lineofsight = lineofsight - 360
    aspectangle_r = aspectangle_n - lineofsight
    return aspectangle_r +360 if aspectangle_r < 0 else aspectangle_r

# Current squint angle calculation is wrong
# def getSquintAngle(orbitlatitude, orbitlongitude, latimage, lonimage, x, y):
#     latlon1 = [latimage[y][x], lonimage[y][x]]
#     latlon2 = [np.mean(orbitlatitude), np.mean(orbitlongitude)]
#     g = Geod(ellps = 'WGS84')
#     #order of object latlons matter, first object represents latlon on ground, second represents satellite
#     #refer to linuxtut.com/en/84d49e18e9dd3373ce0f/ for diagram
#     result = g.inv(latlon1[1], latlon1[0], latlon2[1], latlon2[0]) 
#     squintangle = result[1]
#     print(f'azimuth {result[0]}, back_azimuth {result[1]}')
#     return squintangle

def getGrazingAngle(orbitlatitude, orbitlongitude, orbitheight, lat, lon):   
    """Returns the grazing angle of the satellite with respect to the latlon coordinates of the target.

    Args:
        orbitlatitude: orbitlatitude image of the satellite
        orbitlongitude: orbitlongitude image of the satellite
        orbitheight: orbitheight image of the satellite
        lat: the latitude of the target
        lon: the longitude of the target

    Returns:
       grazingangle_deg: the grazing angle of the satellite with respect to the latlon coordinates of the target in degrees
    """
    latlon1 = [lat, lon]
    
    #average latlon of satellite
    latlon2 = [np.mean(orbitlatitude), np.mean(orbitlongitude)]
    
    g = Geod(ellps = 'WGS84')
    result = g.inv(latlon1[1], latlon1[0], latlon2[1], latlon2[0])
    dist_2d = result[2]
    # print(f'azimuth {result[0]}, back_azimuth {result[1]}')
    grazingangle = np.arctan(np.mean(orbitheight)/dist_2d)
    grazingangle_deg = grazingangle * 180/pi
    # print(f'the distance from the satellite to the object is {dist_2d}m, its grazing angle is {grazingangle} radians or {grazingangle_deg} deg by GPSTIime')
    return grazingangle_deg

def getGPSTime(y, x, orbitimage):
    """Returns the gpstime of the satellite at pixel (x,y). Assumes y is row, x is column.

    Args:
        y: y-pixel
        x: x-pixel
        orbitimage: the orbit image of the satellite

    Returns:
        orbitimage[y][x]: the gps time of the satellite when pixel (x,y) is focused.
    """
    print(orbitimage[y][x])
    return orbitimage[y][x]

def getIntegrationAngle(orbitlatitude, orbitlongitude, orbitheight, lat, lon, gbpgridinfo):
    """Obtains the integration angle of the satellite relative to the latlon coordinates of the target. Angle returned is in degrees

    Args:
        orbitlatitude: the orbit latitude array of the satellite
        orbitlongitude: the orbit longitude array of the satellite
        orbitheight: the orbit height array of the satellite 
        lat: the latitude of the target
        lon: the longitude of the target
        gbpgridinfo: the gbpgridinfo of the given SLC image

    Returns:
        integrationangle: the integration angle with respect to the flight path of the satellite
    """
    #get the two start and end points of the satellite
    latlonstart = [orbitlatitude[0], orbitlongitude[0]]
    latlonend = [orbitlatitude[-1], orbitlongitude[-1]]
    latloncentre = [lat, lon]
    # print(latlonstart, latloncentre, latlonend)
    
    #Using 3D euclidean space
    zone = gbpgridinfo[3]
    zone_letter = 'S' if gbpgridinfo[4] == 1 else 'N'
    latlonstart_east, latlonstart_north, _, _ = from_latlon(latlonstart[0], latlonstart[1], zone, zone_letter)
    latlonend_east, latlonend_north, _, _ = from_latlon(latlonend[0], latlonend[1], zone, zone_letter)
    latloncentre_east, latloncentre_north,_,_ = from_latlon(latloncentre[0], latloncentre[1], zone, zone_letter)
    AB = sqrt((latlonstart_east - latlonend_east)**2 + (latlonstart_north - latlonend_north)**2 + (orbitheight[0] - orbitheight[-1])**2)
    AC = sqrt((latlonstart_east - latloncentre_east)**2 + (latlonstart_north - latloncentre_north)**2 + (orbitheight[0]- 0)**2) #assumes height of point on image is zero
    BC = sqrt((latlonend_east - latloncentre_east)**2 + (latlonend_north - latloncentre_north)**2 + (orbitheight[-1]- 0)**2) #assumes height of point on image is zero
    integrationangle = acos((BC**2 + AC**2 - AB**2)/(2 * BC * AC))
    return integrationangle