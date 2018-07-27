import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from siphon.catalog import TDSCatalog
import urllib
from netCDF4 import Dataset, num2date
from matplotlib import ticker
import matplotlib as mpl
from dateutil import tz
import time
from time import mktime
import os.path
from os import listdir
from os.path import isfile, join
import numpy as np
import matplotlib.image as image
from datetime import datetime, timedelta
from pyproj import Proj
from matplotlib.patches import Rectangle

##
############# Initial Set Up ##################
workdir = "/home/sat_ops/goesR/tpwc/"
imgdir = "/home/sat_ops/goesR/tpwc/conus/"
# Make a new map object for the HRRR model domain map projection
mH = Basemap(resolution='i', projection='lcc', area_thresh=1500, \
            width=1800*3000, height=1060*3100, \
            lat_1=38.5, lat_2=38.5, \
            lat_0=38.5, lon_0=-97.5)
            
# Go to the Unidata Thredds Server for the Current Day
nowdate = datetime.utcnow()
cat = TDSCatalog('http://thredds-jumbo.unidata.ucar.edu/thredds/catalog/satellite/goes16/GOES16/Products/TotalPrecipitableWater/CONUS/current/catalog.xml')
dataset_name = sorted(cat.datasets.keys())[-1]
print dataset_name

if os.path.isfile(imgdir + str(dataset_name) + ".png") == False:
    dataset = cat.datasets[dataset_name]
    
    nc = dataset.remote_access()
    
    geoy = np.array(nc.variables['y'][:]) * 1000.
    geox = np.array(nc.variables['x'][:]) * 1000.
    
    tpwc = np.array(nc.variables['TPW'][:,:]) * nc.variables['TPW'].scale_factor
    add_seconds = int(nc.variables['t'][0])
    DATE = datetime(2000, 1, 1, 12) + timedelta(seconds=add_seconds)
    sat_h = nc.variables['goes_imager_projection'].perspective_point_height
    sat_lon = nc.variables['goes_imager_projection'].longitude_of_projection_origin
    sat_sweep = nc.variables['goes_imager_projection'].sweep_angle_axis
    
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
    X = nc.variables['x'][:] * sat_h
    Y = nc.variables['y'][:] * sat_h
    
    # map object with pyproj
    p = Proj(proj='geos', h=sat_h, lon_0=sat_lon, sweep=sat_sweep)
    # Convert map points to latitude and longitude with the magic provided by Pyproj
    XX, YY = np.meshgrid(X, Y)
    lons, lats = p(XX, YY, inverse=True)
    lats[np.isnan(tpwc)] = np.nan
    lons[np.isnan(tpwc)] = np.nan
    xH, yH = mH(lons, lats)

    ##################################### CONUS Plotting #################################
    rec_height = 120000
    rec_width = mH.xmax
    # Now we can plot the GOES data on the HRRR map domain and projection
    plt.figure(figsize=[16, 12], dpi=100)
    vmin = 0
    vmax = 60
    cmap = "jet"
    # The values of R are ignored becuase we plot the color in colorTuple, but pcolormesh still needs its shape.
    newmap = mH.pcolormesh(xH, yH, tpwc, cmap=cmap, vmin=vmin, vmax=vmax)
    newmap.set_array(None) # without this line, the linewidth is set to zero, but the RGB colorTuple is ignored. I don't know why.
    mH.drawstates(color='k')
    mH.drawcountries()
    mH.drawcoastlines(linewidth=0.7,color='k')
    
    cb = mH.colorbar(location='bottom', size = '2%', pad = '-1.95%', ticks=[5, 15, 25, 35, 45, 55, 65, 75])
    cb.ax.set_xticklabels(['5', '15', '25', '35', '45', '55', '65', '75'])
    cb.outline.set_visible(False) # Remove the colorbar outline
    cb.ax.tick_params(width = 0) # Remove the colorbar ticks 
    cb.ax.xaxis.set_tick_params(pad=-13.75) # Put the colobar labels inside the colorbar
    cb.ax.tick_params(axis='x', colors='black', labelsize=10) # Change the color and size of the colorbar labels
    
    clabeltext='[mm]'
    title = 'NOAA GOES-16 Total Precipitable Water Content'
    timestr = local.strftime('%Y-%m-%d %H:%M ') + et
    currentAxis = plt.gca()
    currentAxis.add_patch(Rectangle((0, 0), 1000000000, rec_height * 1.3, alpha=1, zorder=3, facecolor='darkslateblue'))
    plt.text(9000, 90000,clabeltext,horizontalalignment='left', color = 'white', size=11)
    # btw, 5400000 comes from the width of mH basemap
    currentAxis.add_patch(Rectangle((0, mH.ymax - rec_height), rec_width, rec_height , alpha=1, zorder=3, edgecolor='black',facecolor='white'))
    plt.text(4440000, 3200000,timestr,horizontalalignment='left', color = 'black', size=14)
    plt.text(9000, 3200000,title,horizontalalignment='left', color = 'black', size=14)
    
    # add logo
    im1 = image.imread("/home/sat_ops/goesR/zfolder/udelcemagoes38.png")
    plt.figimage(im1, 15, 15, zorder=1)
    
    # save file
    output_file = imgdir + str(dataset_name) + ".png"
    plt.savefig(output_file, dpi=100, bbox_inches='tight')
    plt.close()
    