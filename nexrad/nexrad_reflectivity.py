import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
#from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import numpy as np
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
import cartopy.crs as ccrs
import metpy
from metpy.plots import colortables as ctables
import xarray as xr

from siphon.catalog import TDSCatalog
import urllib
from netCDF4 import Dataset, num2date
from matplotlib import ticker
import matplotlib as mpl
#import pyart
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
import pandas as pd
############# Initial Set Up ##################
workdir = "/home/sat_ops/goesR/radar/"
imgdir = "/home/sat_ops/goesR/radar/imgconus/"
#site = 'KDOX'
#loc = pyart.io.nexrad_common.get_nexrad_location(site)
lon0 = -75.44 ; lat0 = 38.82556
londix = -74.41111 ; latdix = 39.94694
lonlwx = -77.48751 ; latlwx = 38.97628
lonokx = -72.86444 ; latokx = 40.86556

def nearest(items, pivot):
    return min(items, key=lambda x: abs(x - pivot))

def contrast_correction(color, contrast):
    F = (259*(contrast + 255))/(255.*259-contrast)
    COLOR = F*(color-.5)+.5
    COLOR = np.minimum(COLOR, 1)
    COLOR = np.maximum(COLOR, 0)
    return COLOR
    
######################### NEXRAD #############################
# Go to the Unidata Thredds Server for the Current Day
nowdate = datetime.utcnow()
# https://thredds.ucar.edu/thredds/catalog/grib/NCEP/MRMS/BaseRef/latest.html
cat = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/MRMS/CONUS/BaseRef/catalog.html')

nexrad_name = cat.datasets['Full Collection Dataset']
nexrad = nexrad_name.remote_access(use_xarray=True)
proj_var = nexrad.variables['LatLon_Projection']
created_plot = False

fileind = [-1]
for i in fileind:
    refltime=i
    try:
        refl = nexrad['MergedBaseReflectivityQC_altitude_above_msl'].isel(time=refltime, altitude_above_msl=0)
        timestamp = pd.Timestamp(refl.time.values).to_pydatetime()
    except:
        pass
        try:
            refl = nexrad['MergedBaseReflectivityQC_altitude_above_msl'].isel(time2=refltime, altitude_above_msl=0)
            timestamp = pd.Timestamp(refl.time2.values).to_pydatetime()
        except:
            pass
            try:               
                refl = nexrad['MergedBaseReflectivityQC_altitude_above_msl'].isel(time3=refltime, altitude_above_msl=0)
                timestamp = pd.Timestamp(refl.time3.values).to_pydatetime()
            except:
                pass
                try:               
                    refl = nexrad['MergedBaseReflectivityQC_altitude_above_msl'].isel(time4=refltime, altitude_above_msl=0)
                    timestamp = pd.Timestamp(refl.time4.values).to_pydatetime()
                except:
                    print("NOT SURE WHAT THE TIME DIMENSION IS CALLED")


    geoy = np.array(refl['lat'].values)
    geox = np.array(refl['lon'].values)
    # cf_datetimes kwarg - https://github.com/pvlib/pvlib-python/issues/944
    
    output_file = workdir + 'imgconus/' + str(timestamp.strftime('%Y%m%d_%H%M')) + "nexradC.png"
    if os.path.isfile(output_file) == False:
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
        dBZ = np.array(refl.values)
        dBZ[np.isnan(dBZ)] = 0
        dBZ = np.ma.array(dBZ)
        dBZ[dBZ < 10] = np.ma.masked
        
        #######################################################################
        #######################################################################
        ####################### CONUS Plotting ################################
        #######################################################################
        #######################################################################
        newproj = ccrs.Mercator()
        # Define Plotting Locations
        fs_x = 16
        fs_y = 12
        dpi = 100
        toptext = 0.794
        toptextleft = 0.13
        toptextright = 0.76
        bottomtextleft = 0.13
        bottomtextheight = 0.183
        toprecx = 0.125
        toprecy = 0.786
        bottomrecx = 0.125
        bottomrecy = 0.178
        symbol = u'$\u26A1$'
        
        fig = plt.figure(figsize=[fs_x, fs_y], dpi=dpi)
        ax = fig.add_subplot(1,1,1, projection=newproj)
        rad = ax.pcolormesh(geox, np.flipud(geoy), np.flipud(dBZ), 
          cmap=ctables.get_colortable('NWSReflectivity'),vmax=80, vmin=0,transform=ccrs.PlateCarree())
        ax.set_extent((-65, -128, 21, 50), crs=ccrs.PlateCarree())  
        ax.set_title("")
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                        edgecolor='black', facecolor='none',linewidth=0.5))
        ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                        edgecolor='black', facecolor='none',linewidth=0.5))
        ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                        edgecolor='black', facecolor='none',linewidth=0.5))
                                        
        
        # top rectangle
        fig.patches.extend([plt.Rectangle((toprecx,toprecy),0.775,0.025,
                                      fill=True, alpha=1, facecolor='darkslateblue', zorder=3,
                                      transform=fig.transFigure, figure=fig)])
        fig.patches.extend([plt.Rectangle((bottomrecx,bottomrecy),0.775,0.025,
                              fill=True, facecolor='darkslateblue', zorder=3, alpha=1,
                              transform=fig.transFigure, figure=fig)])
        title = 'NCEP Merged Base Reflectivity'
        timestr = local.strftime('%Y-%m-%d %H:%M ') + et
        
        fig.text(toptextright, toptext,timestr,horizontalalignment='left', color = 'white', size=14, zorder=2000)
        fig.text(bottomtextleft, toptext,title,horizontalalignment='left', color = 'white', size=14, zorder=2000)
        
        im1 = image.imread("/home/sat_ops/goesR/zfolder/udelcemagoes38.png")
        plt.figimage(im1, 23, 58, zorder=1)
        ax.outline_patch.set_visible(False)
        ax.background_patch.set_visible(False)
        # save file
        fig.savefig(output_file, dpi=dpi, bbox_inches='tight', transparent=False)
        plt.close()
    
        ########################################################################
        ################# NOW PLOT MIDATLANTIC DOMAIN ########################
        ########################################################################
        # slice down for faster plotting
        geoy = refl.sel(lat=slice(44, 34.5), lon=slice(271,291)).lat.values
        geox = refl.sel(lat=slice(44, 34.5), lon=slice(271,291)).lon.values
        dBZ = np.array(refl.sel(lat=slice(44, 34.5), lon=slice(271,291)).values)
        dBZ[np.isnan(dBZ)] = 0
        dBZ = np.ma.array(dBZ)
        dBZ[dBZ < 10] = np.ma.masked
        
        fs_x = 8
        fs_y = 8
        dpi = 100
        toptext = 0.861
        textleft = 0.137
        toptextright = 0.69
        bottomtextleft = 0.139
        bottomtextheight = 0.115
        toprecx = 0.1352
        toprecy = 0.854
        bottomrecx = 0.1352
        bottomrecy = 0.11
        
        fig = plt.figure(figsize=[fs_x, fs_y], dpi=dpi)
        ax = fig.add_subplot(1,1,1, projection=newproj)
        rad = ax.pcolormesh(geox, np.flipud(geoy), np.flipud(dBZ), 
          cmap=ctables.get_colortable('NWSReflectivity'),vmax=80, vmin=0,transform=ccrs.PlateCarree())
        ax.set_extent((-69, -81, 34.5, 43.7))
        ax.set_title("")
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                        edgecolor='black', facecolor='none',linewidth=0.5))
        ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                        edgecolor='black', facecolor='none',linewidth=0.5))
        ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                        edgecolor='black', facecolor='none',linewidth=0.5))
    
        # top rectangle
        fig.patches.extend([plt.Rectangle((toprecx,toprecy),0.7538,0.025,
                                      fill=True, alpha=1, facecolor='darkslateblue', zorder=3,
                                      transform=fig.transFigure, figure=fig)])
        fig.patches.extend([plt.Rectangle((bottomrecx,bottomrecy),0.7538,0.025,
                              fill=True, facecolor='darkslateblue', zorder=3, alpha=1,
                              transform=fig.transFigure, figure=fig)])
        title = 'NCEP Merged Base Reflectivity'
        timestr = local.strftime('%Y-%m-%d %H:%M ') + et
        
        fig.text(toptextright, toptext,timestr,horizontalalignment='left', color = 'white', size=10, zorder=2000)
        fig.text(textleft, toptext,title,horizontalalignment='left', color = 'white', size=9, zorder=2000)
        im1 = image.imread("/home/sat_ops/goesR/zfolder/udelcemagoes24.png")
        plt.figimage(im1, 15, 34,zorder=1)
    
        ax.outline_patch.set_visible(False)
        ax.background_patch.set_visible(False)
        output_file = workdir + 'imgmid/' + str(timestamp.strftime('%Y%m%d_%H%M')) + "nexradM.png"
        fig.savefig(output_file, dpi=dpi, bbox_inches='tight', transparent=False)
        plt.close()
        
        created_plot = True

# Now create the gif
import imageio
import numpy as np

if created_plot==True:
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
