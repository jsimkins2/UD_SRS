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
import pyart

############# Initial Set Up ##################
workdir = "/home/sat_ops/goesR/radar/"
imgdir = "/home/sat_ops/goesR/radar/imgconus/"
#site = 'KDOX'
#loc = pyart.io.nexrad_common.get_nexrad_location(site)
lon0 = -75.44 ; lat0 = 38.82556
londix = -74.41111 ; latdix = 39.94694
lonlwx = -77.48751 ; latlwx = 38.97628
lonokx = -72.86444 ; latokx = 40.86556

# Go to the Unidata Thredds Server for the Current Day
nowdate = datetime.utcnow()
cat = TDSCatalog('http://thredds-jumbo.unidata.ucar.edu/thredds/catalog/grib/nexrad/composite/unidata/NEXRAD_Unidata_Reflectivity-' + \
                  str(nowdate.year) + str("%02d"%nowdate.month) + str("%02d"%nowdate.day) + '/catalog.xml')
# run through the last 5 files for now and we'll see whether or not we've already created them or not
raw_list = []
raw_list = list(cat.catalog_refs)[-2:]

dataset_list = []
for r in raw_list:
    if os.path.isfile(imgdir + str(r.split('.')[0]) + ".png") == False:
        dataset_list.append(r)

mH = Basemap(resolution='i', projection='lcc', area_thresh=1500, \
        width=1800*3000, height=1060*3100, \
        lat_1=38.5, lat_2=38.5, \
        lat_0=38.5, lon_0=-97.5)

DH = Basemap(projection='lcc',lon_0=lon0,lat_0=lat0,
    llcrnrlat=lat0-4.5,llcrnrlon=lon0-5.5,
    urcrnrlat=lat0+5,urcrnrlon=lon0+5.5,resolution='h') 

if len(dataset_list) > 0:    
    for i in range(0,len(dataset_list)):
        cat = TDSCatalog('http://thredds-jumbo.unidata.ucar.edu/thredds/catalog/grib/nexrad/composite/unidata/NEXRAD_Unidata_Reflectivity-' + \
                          str(nowdate.year) + str("%02d"%nowdate.month) + str("%02d"%nowdate.day) + '/' + dataset_list[i] + '/catalog.xml')
        
        dataset_name = sorted(cat.datasets.keys())[-1]
        dataset = cat.datasets[dataset_name]
        nc = dataset.remote_access('OPENDAP')
        
        # this is a patch as of 1/3/2019, they changed the structure of the ncdf and got rid of reftime...which is good because 
        # this is a not a forecast product, it's a composite observation product.
        refltime=0
        
        geoy = np.array(nc.variables['y'][:]) * 1000.
        geox = np.array(nc.variables['x'][:]) * 1000.
        refl = np.array(nc.variables['Base_reflectivity_surface_layer'][refltime,:,:])
        proj_var = nc.variables['LambertConformal_Projection']
        time_var = nc.variables['time']
        timestamp = num2date(time_var[:].squeeze(), time_var.units)
        
        
        from_zone = tz.gettz('UTC')
        to_zone = tz.gettz('America/New_York')
        utc = timestamp.replace(tzinfo=from_zone)
        local = utc.astimezone(to_zone)
        
        lt = time.localtime()
        dst = lt.tm_isdst
        lt = time.localtime()
        dst = lt.tm_isdst
        
        if dst == 0:
            et = "EDT"
        else:
            et = "EST"
        
        # need to convert the NaNs to 0 and then mask them in order for them to all be hidden in the image
        dBZ = refl
        dBZ[np.isnan(dBZ)] = 0
        dBZ = np.ma.array(dBZ)
        dBZ[dBZ < 10] = np.ma.masked
        
        # need to find the central projection coordinates give the central lat lon and add those to the geox, geoy arrays to match the new projection with original
        lcc = Proj("+proj=lcc +lat_0=40 +lon_0=260 +lat_1=40 +lat_2=40 +units=km +no_defs +R=6371200.0")
        # create a meshgrid so that they are the same size and we can convert the projection coordinates to lat lon so basemap can understand them
        XX, YY = np.meshgrid(geox, geoy)
        lons, lats = lcc(XX, YY, inverse=True)
        rec_height = 120000
        rec_width = mH.xmax
        # Begin the plotting 
        plt.figure(figsize=[16, 12], dpi=100)
        # Plot the GOES image
        cmap = 'pyart_NWSRef'
        levs = np.linspace(0,80,41,endpoint=True)
        norm = mpl.colors.BoundaryNorm(levs,256)
        # Plot the HRRR reflectivity
        mH.drawcounties(linewidth=0.3, color='lightgray')
        mH.pcolormesh(lons, lats, dBZ, latlon=True,
                      cmap=cmap,
                      vmax=80, vmin=0)
        
        mH.drawcountries(linewidth=0.7,color='k')
        mH.drawstates(linewidth=0.7,color='k')
        mH.drawcoastlines(linewidth=0.7,color='k')
        cb = mH.colorbar(location='bottom', size = '2%', pad = '-1.95%', ticks=[5, 15, 25, 35, 45, 55, 65, 75])
        cb.ax.set_xticklabels(['5', '15', '25', '35', '45', '55', '65', '75'])
        cb.outline.set_visible(False) # Remove the colorbar outline
        cb.ax.tick_params(width = 0) # Remove the colorbar ticks 
        cb.ax.xaxis.set_tick_params(pad=-13.75) # Put the colobar labels inside the colorbar
        cb.ax.tick_params(axis='x', colors='black', labelsize=10) # Change the color and size of the colorbar labels
        
        clabeltext='Reflectivity [dBZ]'
        title = 'NEXRAD II Composite Reflectivity'
        timestr = local.strftime('%Y-%m-%d %H:%M ') + et
        currentAxis = plt.gca()
        currentAxis.add_patch(Rectangle((0, 0), 1000000000, rec_height * 1.3, alpha=1, zorder=3, facecolor='darkslateblue'))
        plt.text(9000, 90000,clabeltext,horizontalalignment='left', color = 'white', size=11)
        # btw, 5400000 comes from the width of mH basemap
        currentAxis.add_patch(Rectangle((0, mH.ymax - rec_height), rec_width, rec_height , alpha=1, zorder=3, edgecolor='black',facecolor='white'))
        plt.text(4440000, 3200000,timestr,horizontalalignment='left', color = 'black', size=14)
        plt.text(9000, 3200000,title,horizontalalignment='left', color = 'black', size=14)
        
        # add logo
        im1 = image.imread("/home/sat_ops/goesR/zfolder/udelcemaunidata38.png")
        plt.figimage(im1, 15, 55, zorder=1)
        
        # save file
        output_file = imgdir + str(dataset_name.split('.')[0]) + ".png"
        plt.savefig(output_file, dpi=100, bbox_inches='tight')
        
        ########################################################################
        ################# NOW PLOT MIDATLANTIC DOMAIN ########################
        ########################################################################
        rec_height = 40000
        rec_width = DH.xmax
        
        plt.figure(figsize=[8, 8], dpi=100)
        # The values of R are ignored becuase we plot the color in colorTuple, but pcolormesh still needs its shape.
        DH.pcolormesh(lons, lats, dBZ, latlon=True,
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
        title = 'NEXRAD II Reflectivity'
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
        plt.figimage(im1, 15, 42,zorder=1)
        
        output_file = workdir + 'imgmid/' + str(dataset_name.split('.')[0]) + ".png"
        plt.savefig(output_file, dpi=100, bbox_inches='tight')
        plt.close()

# Now create the gif
import imageio
import numpy as np
img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
img_names = sorted(img_list)[-15:]

imglen = len(img_names)
images = []
dur_vals = []
for i in range(0,imglen -1):
    if i != imglen:
        dur_vals.append(.07)
dur_vals.append(2)

for i in img_names:
    input_file=imgdir + str(i)
    images.append(imageio.imread(input_file))
imageio.mimsave(workdir + 'radar_conus.gif', images, duration=dur_vals)


###################################################################
imgdir = workdir + 'imgmid/'
img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
img_names = sorted(img_list)[-15:]

imglen = len(img_names)
images = []
dur_vals = []
for i in range(0,imglen -1):
    if i != imglen:
        dur_vals.append(.07)
dur_vals.append(2)

for i in img_names:
    input_file=imgdir + str(i)
    images.append(imageio.imread(input_file))
imageio.mimsave(workdir + 'radar_midatlantic.gif', images, duration=dur_vals)

