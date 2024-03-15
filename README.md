# SAR SLC groundtruthing program
A program used to process SLC images from netCDF4 format and generate geotiff, as well as groundtruth and generation of chips from SNAP.

## How to use
### config.yaml
The config.yaml file is used to specify the netCDF file locations, shpfile location, output folder location and chip parameters.
```
netcdfPathlist:
  - 'data/SAR_CPLX_20190823071330_9.6G_HH_12_pres_2_fdc_246.sar.rgo.sig.nc'
  - 'data/SAR_CPLX_20190823071330_9.6G_VH_11_pres_2_fdc_246.sar.rgo.sig.nc'

shpPath: 'chiptest.shp'
transpose: True
outputFolder: for_training
range_res: 0.7
az_res: 0.7
range_px: 0.1
az_px: 0.1
size: 20
winout: 'hanning'
sat: 2e-2
```
The first time the program is run, there will be no shpfile detected as the geotiff images need to be generated first. A new folder with 4 subfolders will be created at `outputFolder`.
The 4 subfolders:
- `cropped-image`: Cropped chips in geotiff and npy formats. Chips numbered in generation order.
- `cropped-label`: Row and Column coordinates of the mask used to crop the chip relative to the rows and columns of the corresponding geotiff images in `large-label`. Saved as csv.
- `large-image`: Full geotiff images from the processing of netCDF data. Naming convention follows naming convention of original netCDF files.
- `large-label`: Information of each chip:
  - chipnumber
  - latitude and longitude of chip centre
  - size of chip
  - aspect_n
  - aspect_r
  - grazingangle
  - window
  - rangeres and azimuthres
  - rangepixelsize and azimuthpixelsize
  - 4 corner chip coordinates
### Snap
Once geotiffs have been generated, using snap, create a new vector container under `Vector -> New Vector Data Container`. Create polygons to groundtruth targets, saving them
to the vector container created.
Within the vector container, fill in the given columns,
