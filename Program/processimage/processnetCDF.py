import netCDF4 as nc
import os, sys, csv, glob, yaml, geopandas
import numpy as np
from math import sqrt, radians
from utm import to_latlon
sys.path.append(os.path.join(sys.path[0],".."))
from processimage.validimage import getValidImage, getValidBoundaries, recalculateAffineTransformationMatrix
from SLCtogeotiff.adjustSLC import compressImg
from SLCtogeotiff.SLCtogeotiff import convertToGeotiff_affine, convertToGeotiff_gcp
from chips.cutchip import findCornerLatLonGivenDistance, getWindowBoundaries, cropSLC, saveSLC
from PixelLatLongConversion.LatLon_pixel import getPixelfromGeotiff
from anglecalculations.anglecalculation import getIntegrationAngle, getGrazingAngle, getAspectAngle_r
from degradeimage.slcdegrade import degradeSLC, adjustPhase


def getPolygonCoords(shpPath):
    df = geopandas.read_file(shpPath)
    # print(df)
    polygonlist = []
    for polygon in df.geometry:
        polygonlist.append(polygon)
    l = []
    for polygon in polygonlist:
        coords = polygon.exterior.coords
        l.append(list(list(coord) for coord in coords)) 
    return l

def createYamlFromShpFile(configPath):
    data = readYaml(configPath) 
    try:
        shpPath = data['shpPath']
        df = geopandas.read_file(shpPath).to_dict('list')
        dictionary = {'size': data['size'], 'winout': data['winout'], 'sat': float(data['sat']), 'range_res':float(data['range_res']), 'az_res':float(data['az_res']), 'range_px':float(data['range_px']), 'az_px':float(data['az_px'])}
        polygoncoordlist = getPolygonCoords(shpPath)
        centrecoordlist = [np.mean(polycoord[2:], axis= 0).tolist() for polycoord in polygoncoordlist]
        dictionary.update({key:value for key,value in df.items() if key not in ['geometry', 'az_res', 'az_px', 'range_res', 'range_px']})   
        dictionary.update({'centrecoordlist': centrecoordlist})
        dictionary.update({'polygoncoordlist': polygoncoordlist})
        saveLoc = os.path.join(str(data['outputFolder']), 'chips.yaml')
        with open(saveLoc, 'w') as dst:
            dump = yaml.dump(dictionary, default_flow_style= False, sort_keys= False)
            dst.write(dump)
    except Exception as error:
        print(error)
        return None

        
def extractImageName(inputPath):
    """Extracts and returns the name of the geotiff, keeping important information stated in the file name, aquisition date, band, polarization.

    Args:
        inputPath: the path to the netCDF file

    Returns:
        name: the name that the saved images would use
    """
    variables  = inputPath.split('_')
    name = str('')
    name = name.join(variables[2] + '-' + variables[3] + '-' + variables[4])
    if not str.isdigit(variables[5]):
        name = name.join('-' + variables[5])
    return name

def readYaml(yamlPath):
    """reads and returns data from the yaml path

    Args:
        yamlPath: path of the yaml file

    Returns:
        data: the data contained in the yaml file
    """
    with open(yamlPath, 'r') as file:
        data = yaml.safe_load(file)
    return data
        
def processFileAndGetChips(configyamlPath):
    """Processes all images and allow for automatic chip cropping from the config yaml file. If image fed
        in is ground range for col, azimuth/cross range for row, set transpose to True in config.

    Args:
        configyamlPath: the path to the config.yaml file
    """
    
    
    configdata = readYaml(configyamlPath)
    OutputPath = configdata['outputFolder']
    
    geotiffnamelist = dict()
    slcimagelist = dict()
    polarizationlist = dict()
    windowlist = dict()
    gbpgridinfolist = dict()
    orblatimagelist = dict()
    orblonimagelist = dict()
    orbheightimagelist = dict()
    orbitlatitudelist = dict()
    orbitlongitudelist = dict()
    orbitheightlist = dict()
    latimagelist = dict()
    lonimagelist = dict()
    groundrangepixsizelist = dict()
    crossrangepixsizelist = dict()
    azimuthreslist = dict()
    rangereslist = dict()
    wavelengthlist = dict()
    headinglist = dict()
    lookdirectionlist = dict()
    geotiffimageFolder = os.path.join(OutputPath,'large-image')
    geotifflabelFolder = os.path.join(OutputPath, 'large-label')
    chipimageFolder =os.path.join(OutputPath, 'cropped-image')
    chiplabelFolder = os.path.join(OutputPath, 'cropped-label')
    
    #check if the folder already exists
    if not os.path.exists(OutputPath):
        os.makedirs(OutputPath)
        os.makedirs(geotiffimageFolder)
        os.makedirs(geotifflabelFolder)
        os.makedirs(chipimageFolder)
        os.makedirs(chiplabelFolder)
        
    #process each image    
    for inputPath in configdata['geotiffPathlist']:
        geotiffname = extractImageName(inputPath)
        print(f'processing image {geotiffname}')
        file = nc.Dataset(inputPath)

        realimg = file.variables['SigmaImageSingleLookRealPart'][:]
        imagimg = file.variables['SigmaImageSingleLookImaginaryPart'][:]
        fullimg = realimg + 1j*imagimg   
        modeltransformationtag = file.variables['ModelTransformationTag'][:]
        
        #remove dummy values
        validimage, validmpdeltransformationtag = getValidImage(fullimg, modeltransformationtag)
        minrow, maxrow, mincol, maxcol = getValidBoundaries(fullimg)
        
        orblatimage = file.variables['OrbLatImage'][:][minrow:maxrow+1, mincol:maxcol+1]
        orblonimage = file.variables['OrbLonImage'][:][minrow:maxrow+1, mincol:maxcol+1]
        orbheightimage = file.variables['OrbHeightImage'][:][minrow:maxrow+1, mincol:maxcol+1]
        
        orbitlatitude = file.variables['OrbitLatitude'][:]
        orbitlongitude = file.variables['OrbitLongitude'][:]
        orbitheight = file.variables['OrbitHeight'][:]
        orbitheading = file.variables['OrbitHeading'][:]
        
        heading = np.mean(orbitheading)
        c = 299792458 # Speed of light

        gbpgridinfo = file.variables['GBPGridInfo'][:]
        window = file.variables['WindowFunction'][:]
        match window:
            case 0:
                window = 'rect' 
                k = 0.89
            case 1:
                window = 'tri'
            case 2:
                window = 'hamming'
                k = 1.3
            case 3:
                window = 'hanning'
                k = 1.4466
        
        transmittedpolarization = str(file.variables['TxPolarization'][:], 'utf-8')
        receivedpolarization = str(file.variables['RxPolarization'][:], 'utf-8')
        polarization = transmittedpolarization + receivedpolarization
        lookdirection = str(file.variables['LookDirection'][:], 'utf-8')
        transmittedbandwidth = file.variables['TransmittedBandWidth'][:]
        centralfreq = file.variables['CentralFreq'][:]
        # print(transmittedbandwidth, centralfreq)
        
        latimage = file.variables['LatImage'][:][minrow:maxrow+1, mincol:maxcol+1]
        lonimage = file.variables['LonImage'][:][minrow:maxrow+1, mincol:maxcol+1]
        
        
        
        #check if the image has already been processed and saved as a geotiff
        if not os.path.exists(os.path.join(geotiffimageFolder, f'{geotiffname}.tif')):
            compressedimg = compressImg(validimage)
            convertToGeotiff_affine(compressedimg, validmpdeltransformationtag, gbpgridinfo, os.path.join(geotiffimageFolder, f'{geotiffname}.tif'))
            with open(os.path.join(geotifflabelFolder, geotiffname + '.csv'), 'w', newline = '') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(['chipnumber', 'latitude', 'longitude','size', 'aspect_n', 'aspect_r','grazingangle','window', 'rangeres','azimuthres','rangepixelsize', 'azimuthpixelsize', 'chip coordinates'])
        
        #calculate groudn range and cross range pixel size        
        groundrange = file.variables['GroundRange'][:]
        groundrangepixsize = abs(groundrange[0] - groundrange[-1])/ len(groundrange)
        crossrange = file.variables['CrossRange'][:]
        crossrangepixsize = abs(crossrange[0] - crossrange[-1]) / len(crossrange)
        
        #calculate azimuth resolution and range resolution 
        row, col = validimage.shape
        #find centre of full image
        row = int(row/2)
        col = int(col/2)
        teta = getIntegrationAngle(orbitlatitude, orbitlongitude, orbitheight, latimage[row][col], lonimage[row][col], gbpgridinfo)  
        azimuthres = k * (c / (2*transmittedbandwidth))
        # print(k, transmittedbandwidth, azimuthres)
        wavelength = c/centralfreq
        rangeres = k * ( wavelength/ (2 * radians(teta))) 
        
        geotiffnamelist.update({geotiffname:0})
        slcimagelist.update({geotiffname:validimage})
        polarizationlist.update({geotiffname:polarization})
        windowlist.update({geotiffname:window})
        gbpgridinfolist.update({geotiffname:gbpgridinfo})
        orblatimagelist.update({geotiffname:orblatimage})
        orblonimagelist.update({geotiffname:orblonimage})
        orbheightimagelist.update({geotiffname:orbheightimage})
        orbitlatitudelist.update({geotiffname:orbitlatitude})
        orbitlongitudelist.update({geotiffname:orbitlongitude})
        orbitheightlist.update({geotiffname:orbitheight})
        headinglist.update({geotiffname:heading})
        lookdirectionlist.update({geotiffname:lookdirection})
        latimagelist.update({geotiffname:latimage})
        lonimagelist.update({geotiffname:lonimage})
        groundrangepixsizelist.update({geotiffname:groundrangepixsize})
        crossrangepixsizelist.update({geotiffname:crossrangepixsize})
        azimuthreslist.update({geotiffname:azimuthres})
        rangereslist.update({geotiffname:rangeres})
        wavelengthlist.update({geotiffname:wavelength})
    
    #if yaml file exists use it  
    if os.path.exists(os.path.join(OutputPath, 'chips.yaml')):

        d = readYaml(os.path.join(OutputPath, 'chips.yaml'))
        for index, UTM_coords in enumerate(d['centrecoordlist']):
            #get last chip number
            filelist = glob.glob(chipimageFolder + '\\*')
            if filelist:
                latestfile = max(filelist, key=os.path.getmtime)
                name = latestfile.split('\\')[-1]
                chipnumber = int(name.split('-')[-1].split('.')[0]) + 1
            else:
                chipnumber = 0
            print('---------------')
            print('chipnumber: ',chipnumber)
                
            #go through each image
            for geotiffname in geotiffnamelist.keys():
                print('___')
                print('geotiff ',geotiffname)
                slcimage = slcimagelist[geotiffname]
                transpose = configdata['transpose']
                polarization = polarizationlist[geotiffname]
                window = windowlist[geotiffname]
                gbpgridinfo = gbpgridinfolist[geotiffname]
                orblatimage = orblatimagelist[geotiffname]
                orblonimage = orblonimagelist[geotiffname]
                orbheightimage = orbheightimagelist[geotiffname]
                orbitlatitude = orbitlatitudelist[geotiffname]
                orbitlongitude = orbitlongitudelist[geotiffname]
                orbitheight = orbitheightlist[geotiffname]
                latimage = latimagelist[geotiffname]
                lonimage = lonimagelist[geotiffname]
                groundrangepixsize = groundrangepixsizelist[geotiffname]
                crossrangepixsize = crossrangepixsizelist[geotiffname]
                azimuthres = azimuthreslist[geotiffname]
                rangeres = rangereslist[geotiffname]
                wavelength = wavelengthlist[geotiffname]
                lookdirection = lookdirectionlist[geotiffname]
                heading = headinglist[geotiffname]
                
                #convert UTM coords to lat lon coords in geotiff
                lat, lon = to_latlon(UTM_coords[0], UTM_coords[1], gbpgridinfo[3], 'S' if gbpgridinfo[4] == 1 else 'N')
                
                #calculate the position of the chip based on the radius of the earth and the given centre latitude and longitude coordinates
                R = 6371 #radius of earth
                size_km = float(d['size'])/1000
                distance = sqrt((size_km/2)**2 + (size_km/2)**2)
                Ad = distance/R #angular distance
                
                #from centre latlon, we find the top left to crop the window        
                bearing_tl = 315
                lat_tl, lon_tl = findCornerLatLonGivenDistance(float(lat), float(lon), Ad, bearing_tl)
                row_tl, col_tl = getPixelfromGeotiff(lat_tl, lon_tl, os.path.join(geotiffimageFolder, geotiffname + '.tif'))
                row_tl, row_bl, col_tl, col_tr = getWindowBoundaries(slcimage, row_tl, col_tl, crossrangepixsize, groundrangepixsize, float(d['size']))
                lat_tl, lon_tl = latimage[row_tl][col_tl], lonimage[row_tl][col_tl]
                lat_tr, lon_tr = latimage[row_tl][col_tr], lonimage[row_tl][col_tr]
                lat_bl, lon_bl = latimage[row_bl][col_tl], lonimage[row_bl][col_tl]
                lat_br, lon_br = latimage[row_bl][col_tr], lonimage[row_bl][col_tr]

                croppedslc = cropSLC(slcimage,row_tl, row_bl,col_tl, col_tr)
                print('imageshape: ', croppedslc.shape)
                
                
                
                #adjust the phase of the image before performing degradation
                croppedslc_phaseadjusted = adjustPhase(croppedslc, latimage[row_tl:row_bl, col_tl:col_tr],lonimage[row_tl:row_bl, col_tl:col_tr], orblatimage[row_tl:row_bl, col_tl:col_tr], orblonimage[row_tl:row_bl, col_tl:col_tr], orbheightimage[row_tl:row_bl, col_tl:col_tr], wavelength)
                degradedslc,_,_,_,_ ,_ = degradeSLC(croppedslc_phaseadjusted,transpose, azimuthres, rangeres, groundrangepixsize, crossrangepixsize, d['range_res'], d['range_px'], d['az_res'], d['az_px'], window, d['winout'])
                print(degradedslc.shape)
                saveSLC(degradedslc, os.path.join(chipimageFolder, geotiffname + f'-{chipnumber}'))
                degradedslc_compressed = compressImg(degradedslc, sat = float(d['sat']))
                # use either affine transform or ground control points to generate the geotiff image
                # convertToGeotiff_affine(degradedslc_compressed, croppedmodeltransformationtag, gbpgridinfo, os.path.join(chipimageFolder, geotiffname + f'-{chipnumber}.tif'))
                convertToGeotiff_gcp(degradedslc_compressed, latimage[row_tl:row_bl, col_tl: col_tr],lonimage[row_tl:row_bl, col_tl: col_tr], os.path.join(chipimageFolder, geotiffname + f'-{chipnumber}.tif'))
                grazingangle = getGrazingAngle(orbitlatitude, orbitlongitude, orbitheight, lat, lon)
                aspectangle_r = getAspectAngle_r(heading, d['aspect_n'][index], lookdirection)   
                
                #save chip info largelabel
                with open(os.path.join(geotifflabelFolder, geotiffname + '.csv'), 'a', newline = '') as file:
                    writer = csv.writer(file)
                    writer.writerow([int(chipnumber),lat,lon,d['size'],d['aspect_n'][index], aspectangle_r, grazingangle,window,d['range_res'],d['az_res'],d['range_px'],d['az_px'],[[lat_tl, lon_tl], [lat_tr, lon_tr], [lat_br, lon_br], [lat_bl, lon_bl]]])
                #save mask rows and columns in chiplabel
                maskarray = []
                for coordinate in d['polygoncoordlist'][index][2:]:
                    lat, lon = to_latlon(coordinate[0], coordinate[1], gbpgridinfo[3], 'S' if gbpgridinfo[4] == 1 else 'N')
                    row, col = getPixelfromGeotiff(lat, lon, os.path.join(geotiffimageFolder, geotiffname + '.tif'))
                    maskarray.append([row,col])
                with open(os.path.join(chiplabelFolder, geotiffname + f'-{chipnumber}.csv'), 'w', newline='') as file:
                    writer = csv.writer(file)
                    writer.writerow(maskarray)
        print('all chips sucessfully cropped')
        
        
if __name__ == '__main__':
    configpath = 'config.yaml'
    createYamlFromShpFile(configpath)
    processFileAndGetChips(configpath)    
