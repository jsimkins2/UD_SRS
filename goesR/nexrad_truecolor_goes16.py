import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from siphon.catalog import TDSCatalog
import urllib
from netCDF4 import Dataset, num2date
from matplotlib import ticker
import matplotlib as mpl
import pyart
from dateutil import tz
import time
from time import mktime
import os.path
import numpy as np
import matplotlib.image as image
from datetime import datetime, timedelta
from os import listdir
from os.path import isfile, join
import calendar
from pyproj import Proj     
from matplotlib.patches import Rectangle

############# Initial Set Up ##################
workdir = "/home/sat_ops/goesR/radar/"
datadir = "/home/sat_ops/goesR/data/mcmipc/"
imgdir = "/home/sat_ops/goesR/radar/tcconus/"
ltngdir = "/home/sat_ops/goesR/data/glm/"
site = 'KDOX'
################ Grab the Lat/Lon of the site we want ####################
loc = pyart.io.nexrad_common.get_nexrad_location(site)
lon0 = loc[1] ; lat0 = loc[0]
############# Declare Functions ############
def nearest(items, pivot):
    return min(items, key=lambda x: abs(x - pivot))

def contrast_correction(color, contrast):
    F = (259*(contrast + 255))/(255.*259-contrast)
    COLOR = F*(color-.5)+.5
    COLOR = np.minimum(COLOR, 1)
    COLOR = np.maximum(COLOR, 0)
    return COLOR

def JulianDate_to_MMDDYYY(y,jd):
    month = 1
    day = 0
    while jd - calendar.monthrange(y,month)[1] > 0 and month <= 12:
        jd = jd - calendar.monthrange(y,month)[1]
        month = month + 1
    return month,jd,y
    
# Define projections
mH = Basemap(resolution='i', projection='lcc', area_thresh=1500, \
            width=1800*3000, height=1060*3100, \
            lat_1=38.5, lat_2=38.5, \
            lat_0=38.5, lon_0=-97.5)

DH = Basemap(projection='lcc',lon_0=lon0,lat_0=lat0,
    llcrnrlat=lat0-4.5,llcrnrlon=lon0-5.5,
    urcrnrlat=lat0+5,urcrnrlon=lon0+5.5,resolution='h') 

######################### NEXRAD #############################
# Go to the Unidata Thredds Server for the Current Day
nowdate = datetime.utcnow()
cat = TDSCatalog('http://thredds-jumbo.unidata.ucar.edu/thredds/catalog/grib/nexrad/composite/unidata/NEXRAD_Unidata_Reflectivity-' + \
                  str(nowdate.year) + str("%02d"%nowdate.month) + str("%02d"%nowdate.day) + '/catalog.xml')

# run through the last 5 files for now and we'll see whether or not we've already created them or not
raw_list = []
raw_list = list(cat.catalog_refs)[-5:]
raw_list.append(str(cat.catalog_refs[-1]))

dataset_list = []
for r in raw_list:
    if os.path.isfile(workdir + 'tcconus/' + str(r) + ".png") == False:
        dataset_list.append(r)

######################### Lightning Files #############################
GLM_files = [f for f in listdir(ltngdir) if isfile(join(ltngdir, f))]
GLM_names = []

for i in range(0,len(GLM_files)):
    fname = str(GLM_files[i].split("_", 4)[3])
    fname = str(fname.split(".", 2)[0])
    fname = fname[1:]
    GLM_names.append(fname)

lnamelist = sorted(GLM_names, key=int)[-500:]
ldatetime = []
for t in lnamelist:
    jday = t[4:7]
    year = t[0:4]
    mdy = JulianDate_to_MMDDYYY(int(year),int(jday))
    hms = t[7:13]
    t = str(mdy[2]) + str(mdy[0]) + str(mdy[1]) + hms
    ldatetime.append(datetime.strptime(t, '%Y%m%d%H%M%S'))

######################### TRUE COLOR GOES16 #############################
mcmipc_list = sorted([f for f in listdir(datadir) if isfile(join(datadir, f))])

gdatetime = []
for x in mcmipc_list:
    t = x.split('s')[1]
    jday = t[4:7]
    year = t[0:4]
    mdy = JulianDate_to_MMDDYYY(int(year),int(jday))
    hms = t[7:13]
    t = str(mdy[2]) + str(mdy[0]) + str(mdy[1]) + hms
    gdatetime.append(datetime.strptime(t, '%Y%m%d%H%M%S'))

########################### Begin the loop for the matched data and plot image ################################
for i in range(0,len(dataset_list)):
    cat = TDSCatalog('http://thredds-jumbo.unidata.ucar.edu/thredds/catalog/grib/nexrad/composite/unidata/NEXRAD_Unidata_Reflectivity-' + \
                  str(nowdate.year) + str("%02d"%nowdate.month) + str("%02d"%nowdate.day) + '/' + dataset_list[i] + '/catalog.xml')
    dataset_name = sorted(cat.datasets.keys())[-1]
    nexrad_name = cat.datasets[dataset_name]
    nexrad = nexrad_name.remote_access()
    if i==0:
        refltime = len(nexrad.variables['reftime'][:]) - len(dataset_list) + i
        refltimelen = len(nexrad.variables['reftime'][:])
    else:
        refltime = refltime + 1
    
    geoy = np.array(nexrad.variables['y'][:]) * 1000.
    geox = np.array(nexrad.variables['x'][:]) * 1000.
    refl = np.array(nexrad.variables['Base_reflectivity_surface_layer'][refltime,:,:])
    proj_var = nexrad.variables['LambertConformal_Projection']
    time_var = nexrad.variables['time']
    timestamp = num2date(time_var[:].squeeze(), time_var.units)

    ############# Read the goes file ###############
    C_file = mcmipc_list[gdatetime.index(nearest(gdatetime, timestamp[refltime]))]
    # match the lightning files to the goes files
    tem = C_file.split('s')[1]
    jday = tem[4:7]
    year = tem[0:4]
    mdy = JulianDate_to_MMDDYYY(int(year),int(jday))
    hms = tem[7:13]
    mdy_str = str(mdy[2]) + str(mdy[0]) + str(mdy[1]) + hms
    goesdatetime=datetime.strptime(mdy_str, '%Y%m%d%H%M%S')
    
    ltng_index = ldatetime.index(nearest(ldatetime, goesdatetime))
    ltng_files = lnamelist[ltng_index - 15: ltng_index]
    
    Cnight = Dataset(datadir + C_file, 'r')
    
    # Load the RGB arrays
    R = Cnight.variables['CMI_C02'][:].data
    G = Cnight.variables['CMI_C03'][:].data
    B = Cnight.variables['CMI_C01'][:].data
    
    # Turn empty values into nans
    R[R==-1] = np.nan
    G[G==-1] = np.nan
    B[B==-1] = np.nan
    
    # Apply range limits for each channel becuase RGB values must be between 0 and 1
    R = np.maximum(R, 0)
    R = np.minimum(R, 1)
    G = np.maximum(G, 0)
    G = np.minimum(G, 1)
    B = np.maximum(B, 0)
    B = np.minimum(B, 1)
    
    # Apply the gamma correction
    gamma = 0.4
    R = np.power(R, gamma)
    G = np.power(G, gamma)
    B = np.power(B, gamma)
    
    # Calculate the "True" Green
    G_true = 0.48358168 * R + 0.45706946 * B + 0.06038137 * G
    G_true = np.maximum(G_true, 0)
    G_true = np.minimum(G_true, 1)
    
    # Grab the IR / Apply range limits for clean IR channel/Normalize the channel between a range/invert colors/lessen the brightness
    cleanIR = Cnight.variables['CMI_C13'][:].data
    cleanIR[cleanIR==-1] = np.nan
    cleanIR = np.maximum(cleanIR, 90)
    cleanIR = np.minimum(cleanIR, 313)
    cleanIR = (cleanIR-90)/(313-90)
    cleanIR = 1 - cleanIR
    cleanIR = cleanIR/1.5
    
    # Modify the RGB color contrast
    contrast = 125
    RGB_contrast = contrast_correction(np.dstack([R, G_true, B]), contrast)
    RGB_contrast_IR = np.dstack([np.maximum(RGB_contrast[:,:,0], cleanIR), np.maximum(RGB_contrast[:,:,1], cleanIR), np.maximum(RGB_contrast[:,:,2], cleanIR)])
    
    # Satellite Date, height, lon, sweep
    add_seconds = Cnight.variables['t'][0]
    DATE = datetime(2000, 1, 1, 12) + timedelta(seconds=add_seconds)
    sat_h = Cnight.variables['goes_imager_projection'].perspective_point_height
    sat_lon = Cnight.variables['goes_imager_projection'].longitude_of_projection_origin
    sat_sweep = Cnight.variables['goes_imager_projection'].sweep_angle_axis
    
    # Configure EST/EDT depending on time of year
    abi_time = DATE
    from_zone = tz.gettz('UTC')
    to_zone = tz.gettz('America/New_York')
    utc = abi_time.replace(tzinfo=from_zone)
    local = utc.astimezone(to_zone)
    
    lt = time.localtime()
    dst = lt.tm_isdst
    lt = time.localtime()
    dst = lt.tm_isdst
    if dst == 0:
        et = "EDT"
    else:
        et = "EST"
    
    # The projection x and y coordinates equals
    # the scanning angle (in radians) multiplied by the satellite height (http://proj4.org/projections/geos.html)
    X = Cnight.variables['x'][:] * sat_h
    Y = Cnight.variables['y'][:] * sat_h
    
    # map object with pyproj
    p = Proj(proj='geos', h=sat_h, lon_0=sat_lon, sweep=sat_sweep)
    # Convert map points to latitude and longitude with the magic provided by Pyproj
    XX, YY = np.meshgrid(X, Y)
    lons, lats = p(XX, YY, inverse=True)
    lats[np.isnan(R)] = np.nan
    lons[np.isnan(R)] = np.nan

    # Create a color tuple for pcolormesh
    rgb = RGB_contrast_IR[:,:-1,:] # Using one less column is very imporant, else your image will be scrambled! (This is the stange nature of pcolormesh)
    rgb = np.minimum(rgb, 1)
    colorTuple = rgb.reshape((rgb.shape[0] * rgb.shape[1]),3) # flatten array, becuase that's what pcolormesh wants.
    colorTuple = np.insert(colorTuple, 3, 1.0, axis=1) # adding an alpha channel will plot faster?? according to stackoverflow.
    
    # get the nexrad stuff right
    dBZ = refl
    dBZ[np.isnan(dBZ)] = 0
    dBZ = np.ma.array(dBZ)
    dBZ[dBZ < 10] = np.ma.masked
    lcc = Proj("+proj=lcc +lat_0=40 +lon_0=260 +lat_1=40 +lat_2=40 +units=km +no_defs +R=6371200.0")
    # create a meshgrid so that they are the same size and we can convert the projection coordinates to lat lon so basemap can understand them
    XX, YY = np.meshgrid(geox, geoy)
    nexlons, nexlats = lcc(XX, YY, inverse=True)
    cmap = 'pyart_NWSRef'
    levs = np.linspace(0,80,41,endpoint=True)
    norm = mpl.colors.BoundaryNorm(levs,256)

    ##################################### CONUS Plotting #################################
    rec_height = 120000
    rec_width = mH.xmax
    # Now we can plot the GOES data on the HRRR map domain and projection
    plt.figure(figsize=[16, 12], dpi=100)
    xH, yH = mH(lons, lats)
    # The values of R are ignored becuase we plot the color in colorTuple, but pcolormesh still needs its shape.
    newmap = mH.pcolormesh(xH, yH, R, color=colorTuple, linewidth=0)
    newmap.set_array(None) # without this line, the linewidth is set to zero, but the RGB colorTuple is ignored. I don't know why.
    mH.pcolormesh(nexlons, nexlats, dBZ, latlon=True,
          cmap=cmap,
          vmax=80, vmin=0)
    mH.drawstates(color='k')
    mH.drawcountries()
    mH.drawcoastlines(linewidth=0.7,color='k')
    
    cb = mH.colorbar(location='bottom', size = '2%', pad = '-1.95%', ticks=[5, 15, 25, 35, 45, 55, 65, 75])
    cb.ax.set_xticklabels(['5', '15', '25', '35', '45', '55', '65', '75'])
    cb.outline.set_visible(False) # Remove the colorbar outline
    cb.ax.tick_params(width = 0) # Remove the colorbar ticks 
    cb.ax.xaxis.set_tick_params(pad=-13.75) # Put the colobar labels inside the colorbar
    cb.ax.tick_params(axis='x', colors='black', labelsize=10) # Change the color and size of the colorbar labels
    
    clabeltext='Reflectivity [dBZ]'
    title = 'NOAA GOES-16 & NEXRAD II Reflectivity'
    timestr = local.strftime('%Y-%m-%d %H:%M ') + et
    currentAxis = plt.gca()
    currentAxis.add_patch(Rectangle((0, 0), 1000000000, rec_height * 1.3, alpha=1, zorder=3, facecolor='darkslateblue'))
    plt.text(9000, 90000,clabeltext,horizontalalignment='left', color = 'white', size=11)
    # btw, 5400000 comes from the width of mH basemap
    currentAxis.add_patch(Rectangle((0, mH.ymax - rec_height), rec_width, rec_height , alpha=1, zorder=3, edgecolor='black',facecolor='white'))
    plt.text(4440000, 3200000,timestr,horizontalalignment='left', color = 'black', size=14)
    plt.text(9000, 3200000,title,horizontalalignment='left', color = 'black', size=14)
    
    # add logo
    im1 = image.imread("/home/sat_ops/goesR/zfolder/combined38.png")
    #im2 = image.imread("/home/sat_ops/goes_r/nexrad/udel38.png")
    #plt.figimage(im1, 1210, 745, zorder=1)
    plt.figimage(im1, 15, 50, zorder=1)
    
    # save file
    output_file = workdir + 'tcconus/' + str(nexrad_name) + ".png"
    plt.savefig(output_file, dpi=100, bbox_inches='tight')
    plt.close()
    
    ########################################################################
    ################# NOW PLOT MIDATLANTIC DOMAIN ########################
    ########################################################################
        ########## LIGHTNING ##############
    # Now that we have all the satellite ABIs read in, we need to call the lightning files 
    # create dictionaries for the lat/lons
    ltng_lat = {}
    ltng_lon = {}
    
    for lt in range(0, len(ltng_files)):
        ltfile = ltng_files[lt]
        L_file = ltngdir + 'OR_GLM-L2-LCFA_G16_s' + str(ltfile) + '.nc'  # GOES16 East
        L = Dataset(L_file, 'r')
        ltng_lat[lt] = L.variables['flash_lat'][:]
        ltng_lon[lt] = L.variables['flash_lon'][:]
    
    
    for lt in xrange(0, len(ltng_lat)):
        subset = []
        llat=DH.latmin
        hlat=DH.latmax
        llon=DH.lonmin
        rlon=DH.lonmax
        if lt==0:
            midatl_flash_count = 0
        for v in range(0,len(ltng_lat[lt])):
            latval=ltng_lat[lt][v]
            lonval=ltng_lon[lt][v]
            if latval > llat and latval < hlat and lonval > llon and lonval < rlon:
                subset.append(latval)
        midatl_flash_count = midatl_flash_count + len(subset)
    
    rec_height = 40000
    rec_width = DH.xmax
    
    plt.figure(figsize=[8, 8], dpi=100)
    xH, yH = DH(lons, lats)
    # The values of R are ignored becuase we plot the color in colorTuple, but pcolormesh still needs its shape.
    newmap = DH.pcolormesh(xH, yH, R, color=colorTuple, linewidth=0)
    newmap.set_array(None) # without this line, the linewidth is set to zero, but the RGB colorTuple is ignored. I don't know why.
    
    DH.pcolormesh(nexlons, nexlats, dBZ, latlon=True,
          cmap=cmap,
          vmax=80, vmin=0)
    DH.drawcoastlines(linewidth=0.7, color = 'k')
    DH.drawstates(linewidth=0.7, color = 'k')
    cb = DH.colorbar(location='bottom', size = '2%', pad = '-1.95%', ticks=[5, 15, 25, 35, 45, 55, 65, 75])
    cb.ax.set_xticklabels(['5', '15', '25', '35', '45', '55', '65', '75'])
    cb.outline.set_visible(False) # Remove the colorbar outline
    cb.ax.tick_params(width = 0) # Remove the colorbar ticks 
    cb.ax.xaxis.set_tick_params(pad=-10.5) # Put the colobar labels inside the colorbar
    cb.ax.tick_params(axis='x', colors='black', labelsize=7) # Change the color and size of the colorbar labels
    
    clabeltext='Reflectivity [dBZ]'
    title = 'NOAA GOES-16 & NEXRAD II Reflectivity'
    timestr = local.strftime('%B %d, %Y %H:%M ') + et
    timestr = local.strftime('%Y-%m-%d %H:%M ') + et
    currentAxis = plt.gca()
    currentAxis.add_patch(Rectangle((0, 0), DH.xmax, rec_height * 1.2, alpha=1, zorder=3, facecolor='darkslateblue'))
    plt.text(5000, 28000,clabeltext,horizontalalignment='left', color = 'white', size=8)
    # btw, 5400000 comes from the width of mH basemap
    currentAxis.add_patch(Rectangle((0, DH.ymax - rec_height), rec_width, rec_height , alpha=1, zorder=3, edgecolor='black',facecolor='white'))
    plt.text(674000, 1025000,timestr,horizontalalignment='left', color = 'black', size=10)
    plt.text(7000, 1025000,title,horizontalalignment='left', color = 'black', size=10)


    # add logo
    im1 = image.imread("/home/sat_ops/goesR/zfolder/combined24.png")
    #im2 = image.imread("/home/sat_ops/goes_r/nexrad/udel38.png")
    #plt.figimage(im1, 1210, 745, zorder=1)
    plt.figimage(im1, 15, 42,zorder=1)
    # save file
    
    output_file = workdir + 'tcmid/' + str(nexrad_name) + ".png"
    plt.savefig(output_file, dpi=100, bbox_inches='tight')
    plt.close()
    
    ########################## NOW PLOT MID ATLANTIC LIGHTNING WITH THIS MAF #####################
    plt.figure(figsize=[8, 8], dpi=100)
    # The values of R are ignored becuase we plot the color in colorTuple, but pcolormesh still needs its shape.
    newmap = DH.pcolormesh(xH, yH, R, color=colorTuple, linewidth=0)
    newmap.set_array(None) # without this line, the linewidth is set to zero, but the RGB colorTuple is ignored. I don't know why.
    DH.drawcoastlines(linewidth=0.7, color = 'k')
    DH.drawcountries(linewidth=0.7, color = 'k')
    DH.drawstates(linewidth=0.7, color = 'k')
    DH.pcolormesh(nexlons, nexlats, dBZ, latlon=True,
              cmap=cmap,
              vmax=80, vmin=0)
    
    cb = DH.colorbar(location='bottom', size = '2%', pad = '-1.95%', ticks=[5, 15, 25, 35, 45, 55, 65, 75])
    cb.ax.set_xticklabels(['5', '15', '25', '35', '45', '55', '65', '75'])
    cb.outline.set_visible(False) # Remove the colorbar outline
    cb.ax.tick_params(width = 0) # Remove the colorbar ticks 
    cb.ax.xaxis.set_tick_params(pad=-10.5) # Put the colobar labels inside the colorbar
    cb.ax.tick_params(axis='x', colors='black', labelsize=7) # Change the color and size of the colorbar labels
    
    for g in range(0, len(ltng_files) -6):
        group_x, group_y = DH(ltng_lon[g], ltng_lat[g])
        DH.scatter(group_x, group_y, s=11, marker="D", c='white', edgecolor='black', lw=1.2)
    
    clabeltext='Reflectivity [dBZ]                                                                                                  Flash Count=' + str(midatl_flash_count)
    title = 'NOAA GOES-16 & NEXRAD II Reflectivity'
    timestr = local.strftime('%B %d, %Y %H:%M ') + et
    timestr = local.strftime('%Y-%m-%d %H:%M ') + et
    currentAxis = plt.gca()
    currentAxis.add_patch(Rectangle((0, 0), DH.xmax, rec_height * 1.2, alpha=1, zorder=3, facecolor='darkslateblue'))
    plt.text(5000, 28000,clabeltext,horizontalalignment='left', color = 'white', size=8)
    # btw, 5400000 comes from the width of mH basemap
    currentAxis.add_patch(Rectangle((0, DH.ymax - rec_height), rec_width, rec_height , alpha=1, zorder=3, edgecolor='black',facecolor='white'))
    plt.text(674000, 1025000,timestr,horizontalalignment='left', color = 'black', size=10)
    plt.text(7000, 1025000,title + ' & GLM Lightning Flashes',horizontalalignment='left', color = 'black', size=8)


    # add logo
    im1 = image.imread("/home/sat_ops/goesR/zfolder/combined24.png")
    plt.figimage(im1, 15, 42,zorder=1)
    output_file = workdir + "lrgmid/" + str(nexrad_name) + ".png"
    plt.savefig(output_file, dpi=100, bbox_inches='tight')
    Cnight.close()
    
    plt.figimage(im1, 15, 32,zorder=1)


######################## Make gifs for conus ########################
######################## ######################## ######################## 
import imageio
import numpy as np
img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
img_names = sorted(img_list)[-20:]

imglen = len(img_names)
images = []
dur_vals = []
for i in xrange(1,imglen):
    if i != imglen:
        dur_vals.append(.07)
dur_vals.append(2)

for i in img_names:
    input_file=imgdir + str(i)
    images.append(imageio.imread(input_file))
imageio.mimsave(workdir + 'radar_goes_conus.gif', images, duration=dur_vals)

# now for the midatlantic
imgdir = "/home/sat_ops/goesR/radar/tcmid/"

img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
img_names = sorted(img_list)[-20:]

imglen = len(img_names)
images = []
dur_vals = []
for i in xrange(1,imglen):
    if i != imglen:
        dur_vals.append(.07)
dur_vals.append(2)

for i in img_names:
    input_file=imgdir + str(i)
    images.append(imageio.imread(input_file))
imageio.mimsave(workdir + 'radar_goes_midatlantic.gif', images, duration=dur_vals)



# now for the special sauce
imgdir = "/home/sat_ops/goesR/radar/lrgmid/"

img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
img_names = sorted(img_list)[-20:]

imglen = len(img_names)
images = []
dur_vals = []
for i in xrange(1,imglen):
    if i != imglen:
        dur_vals.append(.07)
dur_vals.append(2)

for i in img_names:
    input_file=imgdir + str(i)
    images.append(imageio.imread(input_file))
imageio.mimsave(workdir + 'lightning_radar_goes_midatlantic.gif', images, duration=dur_vals)
