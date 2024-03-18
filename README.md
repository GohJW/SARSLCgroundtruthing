# SAR SLC groundtruthing program
A program used to process SLC images from netCDF4 format and generate geotiff, as well as groundtruth and generation of chips from SNAP.
This repo also contains functions that may not be utilized in the actual ground truthing program.
## How to use
The main python script to run when using this program is `processnetCDF.py` under `processimage` folder.
### config.yaml
The config.yaml file is used to specify the netCDF file locations, shpfile location, output folder location and chip parameters in the following format.
> Note: We assume that the program takes the image as range for row and azimuth for column, if the image is not in this orientation set `transpose` to ```True```
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
- `large-label`: Information of each chip from the corresponding geotiff:
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
> Note: A sample `chiptest.shp` file has been provided in order to test the program. Replace this file with the actual shpfile of the cropped chips when using the program.
### Snap
Once geotiffs have been generated, using snap, create a new vector container under `Vector -> New Vector Data Container`. Create polygons to groundtruth targets, saving them
to the vector container created.
Within the vector container with masks, fill out the individual chip details.
> Note: currently the program only uses the `geometry` and `aspect_n` column, the other parameters are either calculated by the program or set in `config.yaml`. However, the program can be modified to use these inputs by modifying the dictionary in `processnetCDF.py`.

Once the chips are all masked, export the vector container under `Vector -> Export -> Geometry as Shape file`. Save the shape file and change `shpPath` to the location of the shp
file. Running the program with the new shp file will process the images and create a new yaml file `chips.yaml`. This file will contain the parameters of each chip, which will then be read by the program to crop out the corresponding chips.
> Note: Pixel size and resolution are currently fixed for all chips, that means that `netcdfPathlist` needs to be adjusted to crop chips of different sizes and resolutions. E.g.
> cropping for X-band chips first, followed by cropping of KU-band chips with seperate resolution and pixel size.

