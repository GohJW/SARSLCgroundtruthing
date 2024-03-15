import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import math
from pyproj import Geod
import os

import sys
sys.path.append(os.path.join(sys.path[0],".."))
from anglecalculation import getGrazingAngle, getGPSTime

#plotting of satellite trajectory graphs for viewing

filePath = 'data/SAR_CPLX_20190823071330_9.6G_VH_11_pres_2_fdc_246.sar.rgo.sig.nc'
file = nc.Dataset(filePath, 'r')
modeltransformationtag = file.variables['ModelTransformationTag'][:]
gbpgridinfo = file.variables['GBPGridInfo'][:]

orbitlatitude = file.variables['OrbitLatitude'][:]
orbitlongitude = file.variables['OrbitLongitude'][:]
orbitheight = file.variables['OrbitHeight'][:]
orbitheading = file.variables['OrbitHeading'][:]
gpstime = file.variables['GPSTime'][:]
data = np.column_stack((orbitlatitude, orbitlongitude))
orbitheight = file.variables['OrbitHeight'][:]
latimage = file.variables['LatImage'][:]
lonimage = file.variables['LonImage'][:]
[y,x] = latimage.shape
y = int(y/2)
x = int(x/2)
plt.plot(lonimage[y][x], latimage[y][x], 'x')
#plotting of graphs
plt.plot(orbitlongitude, orbitlatitude)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Plot of Latitude against Longitude')
plt.axhline(np.mean(orbitlatitude), color = 'black', ls = '--')
plt.ticklabel_format(useOffset=False, style='plain')
plt.tight_layout()
plt.savefig('images/satellitetrajectorywithcentreimagelatlon.png')
plt.show()

fig = plt.figure()
gs = fig.add_gridspec(4, hspace = 0)
axs = gs.subplots(sharex=True, sharey=False)
fig.suptitle('Against GPSTime')
axs[0].plot(gpstime, orbitlatitude) 
axs[0].set_ylabel('Latitude')
axs[0].axhline(np.mean(orbitlatitude), color = 'black', ls = '--')
axs[0].ticklabel_format(useOffset=False, style='plain')
axs[1].plot(gpstime, orbitlongitude)
axs[1].set_ylabel('Longitude')
axs[1].axhline(np.mean(orbitlongitude), color = 'black', ls = '--')
axs[2].plot(gpstime, orbitheight)
axs[2].set_ylabel('Height')
axs[2].axhline(np.mean(orbitheight), color = 'black', ls = '--')
axs[3].plot(gpstime, orbitheading)
axs[3].set_ylabel('Heading')
axs[3].set_xlabel('GPSTime')
axs[3].axhline(np.mean(orbitheading), color = 'black', ls = '--')
fig.tight_layout()

plt.savefig('images/againsttime.png')
plt.show()

#original grazing angle calculation
grazingangle_mean = getGrazingAngle(orbitlatitude, orbitlongitude, orbitheight, latimage, lonimage)

#calc by instantaneous satellite position
#get centre of image
x, y = latimage.shape[1], latimage.shape[0]
x = int(x/2)
y= int(y/2)
orbitimage = file.variables['OrbitImage'][:]
currgpstime_index = orbitimage[y][x]
print(currgpstime_index)
index_size = gpstime.size
olat = np.interp(currgpstime_index, range(0, index_size), orbitlatitude)
olon = np.interp(currgpstime_index, range(0, index_size), orbitlongitude)
oheight = np.interp(currgpstime_index, range(0, index_size), orbitheight)
latlon1 = [latimage[y][x], lonimage[y][x]]
latlon2 = [olat, olon]
g = Geod(ellps = 'WGS84')
result = g.inv(latlon1[1], latlon1[0], latlon2[1], latlon2[0])
dist_2d = result[2]
grazingangle = np.arctan(oheight/dist_2d)
grazingangle_instant = grazingangle * 180/math.pi
print(grazingangle_mean, grazingangle_instant)