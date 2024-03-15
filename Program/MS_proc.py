## Script created to process MS data to a state where phase history can be manipulated.
# To note: Guesswork and luck in getting the signal model, which is pending confirmation by vendor.

import netCDF4 as nc
import numpy as np
import os, sys

# from scipy.io import savemat

file = nc.Dataset('.\SAR_CPLX_20190823071330_9.6G_HH_12_pres_2_fdc_246.sar.rgo.sig.nc')

#Reading all variables into a metadata structure
metadata = {}
for var in file.variables:
    metadata[var] = file.variables[var]

#Image Extraction
realimg = metadata['SigmaImageSingleLookRealPart'][:]
imagimg = metadata['SigmaImageSingleLookImaginaryPart'][:]
fullimg = realimg + 1j*imagimg

#Extract other important variables
LatImage = metadata['LatImage']
LonImage = metadata['LonImage']
Alt = 0 #Assume ground plane is at 0 height
OrbLatImage = metadata['OrbLatImage']
OrbLonImage = metadata['OrbLonImage']
OrbAlt = np.mean(metadata['OrbitHeight'])

#Image cropping (to extract only valid portions)
[row,col] = np.where(realimg != -9999)
minrow = min(row)
maxrow = max(row)
mincol = min(col)
maxcol = max(col)-3
fullimg = fullimg[minrow:maxrow+1, mincol:maxcol +1]
LatImage = LatImage[minrow:maxrow+1, mincol:maxcol +1]
LonImage = LonImage[minrow:maxrow+1, mincol:maxcol +1]
OrbLatImage = OrbLatImage[minrow:maxrow+1, mincol:maxcol +1]
OrbLonImage = OrbLonImage[minrow:maxrow+1, mincol:maxcol +1]

# Compute range of each pixel to Orbit center
# Note: Error in using sin and cos given that input should be radians.
# But correct result is obtained. To be checking with vendor with regards to signal model
r = 6378137.0
r1 = r + Alt
r2 = r + OrbAlt

x1 = r1 * np.sin(LonImage) * np.cos(LatImage)
y1 = r1 * np.sin(LonImage) * np.sin(LatImage)
z1 = r1 * np.cos(LonImage)

x2 = r2 * np.sin(OrbLonImage) * np.cos(OrbLatImage)
y2 = r2 * np.sin(OrbLonImage) * np.sin(OrbLatImage)
z2 = r2 * np.cos(OrbLonImage)

# Due to uncertainties in the signal model, both R_pix using 3D distance and R_perp using 2D perpendicular distance were used.
R_pix = np.linalg.norm([x1 - x2, y1 - y2, z1 - z2], axis = 0)
R_perp = np.linalg.norm([y1 - y2, z1 - z2], axis = 0)

# Phase compensation
# Note: Tested two different phase compensations based on different range definitions, and generally prefer phase_comp2 at the moment
c = 299792458 # Speed of light
wl = c/np.array(metadata['CentralFreq'])
phase_comp = np.exp(-1j * 4 * np.pi / wl * R_pix)
phase_comp2 = np.exp(-1j * 4 * np.pi / wl * R_perp)

img_comp1 = fullimg*phase_comp
img_comp2 = fullimg*phase_comp2

# Saving parameters for viewing in Matlab
# data = {'phase_comp2': phase_comp2, 'img_comp2': img_comp2, 'fullimg2': fullimg, 'R_perp2': R_perp, 'lambda2': wl}
# savemat('test_data.mat', mdict = data, appendmat = False)

exit(1)

