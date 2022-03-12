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
datadir = "/home/sat_ops/goesR/data/mcmipc/"
imgdir = "/home/sat_ops/goesR/radar/tcconus/"
ltngdir = "/home/sat_ops/goesR/data/glm/"
site = 'KDOX'
################ Grab the Lat/Lon of the site we want ####################
#loc = pyart.io.nexrad_common.get_nexrad_location(site)
#lon0 = loc[1] ; lat0 = loc[0]
############# Declare Functions ############

def nearest(items, pivot):
    return min(items, key=lambda x: abs(x - pivot))

def contrast_correction(color, contrast):
    F = (259*(contrast + 255))/(255.*259-contrast)
    COLOR = F*(color-.5)+.5
    COLOR = np.minimum(COLOR, 1)
    COLOR = np.maximum(COLOR, 0)
    return COLOR
    
######################### Lightning Files #############################
GLM_files = [f for f in listdir(ltngdir) if isfile(join(ltngdir, f))]

lnamelist = sorted(GLM_files)[-500:]
ldatetime = []
for t in lnamelist:
    t = t.split('s')[1][:-4]
    ldatetime.append(datetime.strptime(t, '%Y%j%H%M%S'))

######################### TRUE COLOR GOES16 #############################
mcmipc_list = sorted([f for f in listdir(datadir) if isfile(join(datadir, f))])

gdatetime = []
for x in mcmipc_list:
    t = x.split('s')[1][:-4]
    gdatetime.append(datetime.strptime(t, '%Y%j%H%M%S'))

########################### Begin the loop for the matched data and plot image ################################

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
    output_file = workdir + 'tcconus/' + str(timestamp.strftime('%Y%m%d_%H%M')) + "nexradTCC.png"
    if os.path.isfile(output_file) == False:
        ############# Read the goes file ###############
        C_file = mcmipc_list[gdatetime.index(nearest(gdatetime, timestamp))]
            # match the lightning files to the goes files
        tem = C_file.split('s')[1][:-4]
        goesdatetime=datetime.strptime(tem, '%Y%j%H%M%S')
        
        ltng_index = ldatetime.index(nearest(ldatetime, goesdatetime))
        ltng_files = lnamelist[ltng_index - 15: ltng_index]
        if ltng_index ==0:
            ltng_files = lnamelist[ltng_index - 15:]
        
        # open our Combined (C) file for goes16
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
        
        ########## LIGHTNING ##############
        # Now that we have all the satellite ABIs read in, we need to call the lightning files 
        # create dictionaries for the lat/lons
        ltng_lat = {}
        ltng_lon = {}
        
        # open up the lightning files that have been matched with the goes16 file
        # record the lat/lon location of the lightning flash and place them in ltng_lat/lon
        for lt in range(0, len(ltng_files)):
            ltfile = ltng_files[lt]
            L_file = ltngdir + str(ltfile)  # GOES16 East
            L = Dataset(L_file, 'r')
            ltng_lat[lt] = L.variables['flash_lat'][:]
            ltng_lon[lt] = L.variables['flash_lon'][:]
            if lt==0:
                ltngxr = xr.open_dataset(L_file)
                ltngxr = ltngxr.metpy.parse_cf("flash_energy")
                ltngproj = ltngxr.metpy.cartopy_crs
        
        # tell flash count it's zero if there are no strikes 
        for lt in range(0, len(ltng_lat)):
            if lt==0:
                conus_flash_count = 0
            ylt = len(ltng_lat[lt])
            conus_flash_count = conus_flash_count + ylt
        
        # perform the same operations for mid atlantic flash counts
        for lt in range(0, len(ltng_lat)):
            subset = []
            llat=34.32556
            hlat=43.82556
            llon=-80.94
            rlon=-69.94
            if lt==0:
                midatl_flash_count = 0
            for v in range(0,len(ltng_lat[lt])):
                latval=ltng_lat[lt][v]
                lonval=ltng_lon[lt][v]
                if latval > llat and latval < hlat and lonval > llon and lonval < rlon:
                    subset.append(latval)
            midatl_flash_count = midatl_flash_count + len(subset)
    
        # Modify the RGB color contrast to enhance the image
        contrast = 125
        RGB_contrast = contrast_correction(np.dstack([R, G_true, B]), contrast)
        RGB_contrast_IR = np.dstack([np.maximum(RGB_contrast[:,:,0], cleanIR), np.maximum(RGB_contrast[:,:,1], cleanIR), np.maximum(RGB_contrast[:,:,2], cleanIR)])
        # Create a color tuple for pcolormesh
        rgb = RGB_contrast_IR[:,:,:] # key thing here...we are NOT using one less column like before, pcolormesh works now....
        rgb = np.minimum(rgb, 1)
        colorTuple = rgb.reshape((rgb.shape[0] * rgb.shape[1]),3) # flatten array, becuase that's what pcolormesh wants.
        colorTuple = np.insert(colorTuple, 3, 1.0, axis=1) # adding an alpha channel will plot faster?? according to stackoverflow.
        
        # load in the same file to an xarray object so we can gather information from it. This is a quick and easy way to do this
        Cnight2 = xr.open_dataset(datadir + C_file)
        dat = Cnight2.metpy.parse_cf("CMI_C01")
        goesproj = dat.metpy.cartopy_crs
        newproj = ccrs.Mercator()
    
        # Figure out the time
        ymd = Cnight2.time_coverage_end.split("T")[0]
        hms = Cnight2.time_coverage_end.split("T")[1][:-3]
        gtimestamp = ymd + " " + hms
        gtimestamp = datetime.strptime(gtimestamp, "%Y-%m-%d %H:%M:%S")
        from_zone = tz.gettz('UTC')
        to_zone = tz.gettz('America/New_York')
        utc = gtimestamp.replace(tzinfo=from_zone)
        local = utc.astimezone(to_zone)
        
        lt = time.localtime()
        dst = lt.tm_isdst
        lt = time.localtime()
        dst = lt.tm_isdst
        
        if dst == 0:
            et = "EST"
        else:
            et = "EDT"
    
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
        im = ax.pcolormesh(dat['x'], dat['y'], R, color=colorTuple, transform=goesproj)
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
        title = 'NOAA GOES16 True Color & NCEP Merged Base Reflectivity'
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
        im = ax.pcolormesh(dat['x'], dat['y'], R, color=colorTuple, transform=goesproj)
        rad = ax.pcolormesh(geox, np.flipud(geoy), np.flipud(dBZ), 
          cmap=ctables.get_colortable('NWSReflectivity'),vmax=80, vmin=0,transform=ccrs.PlateCarree())
        ax.set_extent((-69, -81, 34.5, 44))
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
        title = 'NOAA GOES16 True Color & NCEP Merged Base Reflectivity'
        timestr = local.strftime('%Y-%m-%d %H:%M ') + et
        
        fig.text(toptextright, toptext,timestr,horizontalalignment='left', color = 'white', size=10, zorder=2000)
        fig.text(textleft, toptext,title,horizontalalignment='left', color = 'white', size=9, zorder=2000)
        im1 = image.imread("/home/sat_ops/goesR/zfolder/udelcemagoes24.png")
        plt.figimage(im1, 15, 34,zorder=1)
    
        ax.outline_patch.set_visible(False)
        ax.background_patch.set_visible(False)
        output_file = workdir + 'tcmid/' + str(timestamp.strftime('%Y%m%d_%H%M')) + "nexradTCM.png"
        fig.savefig(output_file, dpi=dpi, bbox_inches='tight', transparent=False)
        plt.close()
    
        ########################## NOW PLOT MID ATLANTIC LIGHTNING #####################
        ########################## NOW PLOT MID ATLANTIC LIGHTNING #####################
        ########################## NOW PLOT MID ATLANTIC LIGHTNING #####################
    
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
        im = ax.pcolormesh(dat['x'], dat['y'], R, color=colorTuple, transform=goesproj)
        rad = ax.pcolormesh(geox, np.flipud(geoy), np.flipud(dBZ), 
          cmap=ctables.get_colortable('NWSReflectivity'),vmax=80, vmin=0,transform=ccrs.PlateCarree())
        for g in range(0, len(ltng_lat)):
            ax.scatter(ltng_lon[g], ltng_lat[g],s=11, marker="D", c='white', edgecolor='k', lw=0.5, transform=ltngproj)
        ax.set_extent((-69, -81, 34.5, 44))
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
        title = 'NOAA GOES16, NCEP MRMS, & NOAA GLM - Powered by CEMA'
        timestr = local.strftime('%Y-%m-%d %H:%M ') + et
        
        fig.text(toptextright, toptext,timestr,horizontalalignment='left', color = 'white', size=10, zorder=2000)
        fig.text(textleft, toptext,title,horizontalalignment='left', color = 'white', size=9, zorder=2000)
        try:
            midatl_flash_count
        except NameError:
            midatl_flash_count = 'Unavailable'
        clabeltext = 'Flash Count=' + str(midatl_flash_count)
        fig.text(bottomtextleft, bottomtextheight,clabeltext,horizontalalignment='left', color = 'white', size=10, zorder=2000)
        im1 = image.imread("/home/sat_ops/goesR/zfolder/udelcemagoes24.png")
        plt.figimage(im1, 15, 34,zorder=1)
    
        ax.outline_patch.set_visible(False)
        ax.background_patch.set_visible(False)
        output_file = workdir + 'lrgmid/' + str(timestamp.strftime('%Y%m%d_%H%M')) + "nexradTCL.png"
        fig.savefig(output_file, dpi=dpi, bbox_inches='tight', transparent=False)
        plt.close()

        created_plot = True


######################## Make gifs ########################
######################## ######################## ######################## 
import imageio
import numpy as np
if created_plot==True:
    img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
    img_names = sorted(img_list)[-20:]
    
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
    imageio.mimsave(workdir + 'radar_goes_conus.gif', images, duration=dur_vals)
    
    # now for the midatlantic
    imgdir = "/home/sat_ops/goesR/radar/tcmid/"
    
    img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
    img_names = sorted(img_list)[-20:]
    
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
    imageio.mimsave(workdir + 'radar_goes_midatlantic.gif', images, duration=dur_vals)
    
    
    
    # now for the special sauce
    imgdir = "/home/sat_ops/goesR/radar/lrgmid/"
    
    img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
    img_names = sorted(img_list)[-20:]
    
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
    imageio.mimsave(workdir + 'lightning_radar_goes_midatlantic.gif', images, duration=dur_vals)