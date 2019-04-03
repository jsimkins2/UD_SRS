# This script creates true color imagery for CONUS and Mid Atlantic & Lightning/TC over CONUS and Mid Atlantic
# 4 different images altogether 
# James Simkins
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
from matplotlib import patheffects, ticker
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.image as image

from siphon.catalog import TDSCatalog
import urllib
from netCDF4 import Dataset, num2date
import xarray as xr
from xarray.backends import NetCDF4DataStore
import numpy as np
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
import cartopy.crs as ccrs
import metpy
from metpy.plots import colortables
import pyart

from dateutil import tz
import time
from time import mktime
import os.path
from datetime import datetime, timedelta
from os import listdir
from os.path import isfile, join
import calendar
import imageio
#suppress deprecation warnings
import warnings
warnings.simplefilter("ignore", category=DeprecationWarning)


############# Initial Set Up ##################
workdir = "/home/sat_ops/goesR/"
datadir = "/home/sat_ops/goesR/data/mcmipc/"
ltngdir = "/home/sat_ops/goesR/data/glm/"
imgdir = "/home/sat_ops/goesR/ltng_conus/"
site = 'KDOX'

############## Define Functions that will be used later #########################
# initial block mean allows us to downsample the Red Band which is a higher resolution than the others
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


################ Grab the Lat/Lon of the site we want ####################
# Note, this is for the regional map
loc = pyart.io.nexrad_common.get_nexrad_location(site)
lon0 = loc[1] ; lat0 = loc[0]


# create a logfile with most recent files created in the shell script
file_names = [f for f in listdir(datadir) if isfile(join(datadir, f))]

fnamelist = []
for i in range(0,len(file_names)):
    fname = str(file_names[i].split("_", 4)[3])
    fname = str(fname.split(".", 2)[0])
    fname = fname[1:]
    fnamelist.append(fname)

fnamelist = sorted(fnamelist, key=int)
GLM_files = [f for f in listdir(ltngdir) if isfile(join(ltngdir, f))]
GLM_names = []

for i in range(0,len(GLM_files)):
    fname = str(GLM_files[i].split("_", 4)[3])
    fname = str(fname.split(".", 2)[0])
    fname = fname[1:]
    GLM_names.append(fname)

lnamelist = sorted(GLM_names, key=int)
# sort through the files and make sure that they don't exist before we process them
ldatetime = []
for t in lnamelist:
    jday = t[4:7]
    year = t[0:4]
    mdy = JulianDate_to_MMDDYYY(int(year),int(jday))
    hms = t[7:13]
    t = str(mdy[2]) + jday + hms
    ldatetime.append(datetime.strptime(t, '%Y%j%H%M%S'))

ABI_datetime = []
for i in fnamelist:    
    if os.path.isfile(imgdir + str(i) + ".png") == False:
        ABI_datetime.append(i)

# begin the loop that makes the images
if len(ABI_datetime) > 0:
    for n in range(0, len(ABI_datetime)):
        t = ABI_datetime[n][:-3]
        gdatetime=datetime.strptime(t, '%Y%j%H%M')
        ltng_index = ldatetime.index(nearest(ldatetime, gdatetime))
        if ltng_index < 15:
            ltng_files = lnamelist[ltng_index - ltng_index + 1: ltng_index]
        else:
            ltng_files = lnamelist[ltng_index - 15: ltng_index]
        
        C_file = datadir + 'OR_ABI-L2-MCMIPC-M3_G16_s' + str(ABI_datetime[n]) + '.nc'  # GOES16 East
        if os.path.isfile(C_file) == False:
            C_file = datadir + 'OR_ABI-L2-MCMIPC-M4_G16_s' + str(ABI_datetime[n]) + '.nc'
            if os.path.isfile(C_file) == False:
                C_file = datadir + 'OR_ABI-L2-MCMIPC-M6_G16_s' + str(ABI_datetime[n]) + '.nc'
                if os.path.isfile(C_file) == False:
                    import smtplib
                    dfile=open("/home/sat_ops/goes_r/noaa_format/udelsatinfo.txt", "r")
                    dw = dfile.readline()
                    server = smtplib.SMTP('smtp.gmail.com', 587)
                    server.starttls()
                    server.login("goessatelliteudel@gmail.com", dw)
                    msg = "TC CONUS IS BREAKING"
                    server.sendmail("goessatelliteudel@gmail.com", "simkins@udel.edu", msg)
                    server.quit()
        
        # below is code from Brian Blaylock of Univ Utah, thanks Brian!
        Cnight = Dataset(C_file, 'r')
        
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
        
        for lt in range(0, len(ltng_files)):
            ltfile = ltng_files[lt]
            L_file = ltngdir + 'OR_GLM-L2-LCFA_G16_s' + str(ltfile) + '.nc'  # GOES16 East
            L = Dataset(L_file, 'r')
            ltng_lat[lt] = L.variables['flash_lat'][:]
            ltng_lon[lt] = L.variables['flash_lon'][:]
            if lt==0:
                ltngxr = xr.open_dataset(L_file)
                ltngxr = ltngxr.metpy.parse_cf("flash_energy")
                ltngproj = ltngxr.metpy.cartopy_crs
        
        for lt in range(0, len(ltng_lat)):
            if lt==0:
                conus_flash_count = 0
            ylt = len(ltng_lat[lt])
            conus_flash_count = conus_flash_count + ylt
        
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

        # Modify the RGB color contrast
        contrast = 125
        RGB_contrast = contrast_correction(np.dstack([R, G_true, B]), contrast)
        RGB_contrast_IR = np.dstack([np.maximum(RGB_contrast[:,:,0], cleanIR), np.maximum(RGB_contrast[:,:,1], cleanIR), np.maximum(RGB_contrast[:,:,2], cleanIR)])
        # Create a color tuple for pcolormesh
        rgb = RGB_contrast_IR[:,:-1,:] # Using one less column is very imporant, else your image will be scrambled! (This is the stange nature of pcolormesh)
        rgb = np.minimum(rgb, 1)
        colorTuple = rgb.reshape((rgb.shape[0] * rgb.shape[1]),3) # flatten array, becuase that's what pcolormesh wants.
        colorTuple = np.insert(colorTuple, 3, 1.0, axis=1) # adding an alpha channel will plot faster?? according to stackoverflow.
        
        # load in the same file to an xarray object so we can gather information from it. This is a quick and easy way to do this
        Cnight2 = xr.open_dataset(C_file)
        dat = Cnight2.metpy.parse_cf("CMI_C01")
        proj = dat.metpy.cartopy_crs
        newproj = ccrs.Mercator()

        # Figure out the time
        ymd = Cnight2.time_coverage_end.split("T")[0]
        hms = Cnight2.time_coverage_end.split("T")[1][:-3]
        timestamp = ymd + " " + hms
        timestamp = datetime.strptime(timestamp, "%Y-%m-%d %H:%M:%S")
        from_zone = tz.gettz('UTC')
        to_zone = tz.gettz('America/New_York')
        utc = timestamp.replace(tzinfo=from_zone)
        local = utc.astimezone(to_zone)
        
        lt = time.localtime()
        dst = lt.tm_isdst
        lt = time.localtime()
        dst = lt.tm_isdst
        
        if dst == 0:
            et = "EST"
        else:
            et = "EDT"
    
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
        bottomtextheight = 0.212
        toprecx = 0.125
        toprecy = 0.786
        bottomrecx = 0.125
        bottomrecy = 0.19
        symbol = u'$\u26A1$'
        
        fig = plt.figure(figsize=[fs_x, fs_y], dpi=dpi)
        ax = fig.add_subplot(1,1,1, projection=newproj)
        im = ax.pcolormesh(dat['x'], dat['y'], R, color=colorTuple, transform=proj)
        for g in range(0, len(ltng_lat)):
            ax.scatter(ltng_lon[g], ltng_lat[g], s=18, marker=symbol, c='red', edgecolor='red', lw=0, transform=ltngproj)
        fig.patches.extend([plt.Rectangle((toprecx,toprecy),0.7745,0.025,
                              fill=True, color='black', alpha=1, zorder=1000,
                              transform=fig.transFigure, figure=fig)])
        ax.set_extent((-65, -128, 21, 47))
        ax.set_title("")
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                        edgecolor='black', facecolor='none',linewidth=0.5))
        ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                        edgecolor='black', facecolor='none',linewidth=0.5))
        ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                        edgecolor='black', facecolor='none',linewidth=0.5))
                                        
        
        # top rectangle
        fig.patches.extend([plt.Rectangle((toprecx,toprecy),0.7745,0.025,
                                      fill=True, color='black', alpha=1, zorder=1000,
                                      transform=fig.transFigure, figure=fig)])
        
        title = 'NOAA GOES16 True Color & Lightning Flashes - Powered By CEMA'
        timestr = local.strftime('%Y-%m-%d %H:%M ') + et
        
        fig.text(toptextright, toptext,timestr,horizontalalignment='left', color = 'white', size=14, zorder=2000)
        fig.text(bottomtextleft, toptext,title,horizontalalignment='left', color = 'white', size=14, zorder=2000)
        
        try:
            conus_flash_count
        except NameError:
            conus_flash_count = 'Unavailable'
        
        clabeltext = 'Flash Count=' + str(conus_flash_count)
        fig.text(.13, .19,clabeltext,horizontalalignment='left', color = 'red', size=12, zorder=2000)
        
        ax.outline_patch.set_visible(False)
        ax.background_patch.set_visible(False)
        output_file = workdir + "ltng_conus/" + ABI_datetime[n] + ".png"
        fig.savefig(output_file, dpi=dpi, bbox_inches='tight', transparent=False)
        plt.close()

        #######################################################################
        #######################################################################
        ####################### Mid Atlantic ##################################
        #######################################################################
        #######################################################################
        fs_x = 8
        fs_y = 8
        dpi = 100
        toptext = 0.863
        textleft = 0.137
        toptextright = 0.68
        bottomtextleft = 0.13
        bottomtextheight = 0.212
        toprecx = 0.1355
        toprecy = 0.859
        bottomrecx = 0.125
        bottomrecy = 0.205
        
        fig = plt.figure(figsize=[fs_x, fs_y], dpi=dpi)
        ax = fig.add_subplot(1,1,1, projection=newproj)
        im = ax.pcolormesh(dat['x'], dat['y'], R, color=colorTuple, transform=proj)
        for g in range(0, len(ltng_lat)):
            ax.scatter(ltng_lon[g], ltng_lat[g], s=18, marker=symbol, c='red', edgecolor='red', lw=0, transform=ltngproj)
        ax.set_extent((-69, -81, 34.5, 44))
        ax.set_title("")
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                        edgecolor='black', facecolor='none',linewidth=0.5))
        ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                        edgecolor='black', facecolor='none',linewidth=0.5))
        ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                        edgecolor='black', facecolor='none',linewidth=0.5))

        # top rectangle
        fig.patches.extend([plt.Rectangle((toprecx,toprecy),0.7535,0.025,
                              fill=True, color='black', alpha=1, zorder=1000,
                              transform=fig.transFigure, figure=fig)])
        title = 'NOAA GOES16 True Color & Lightning Flashes - Powered By CEMA'
        timestr = local.strftime('%Y-%m-%d %H:%M ') + et
        
        fig.text(toptextright, toptext,timestr,horizontalalignment='left', color = 'white', size=10, zorder=2000)
        fig.text(textleft, toptext,title,horizontalalignment='left', color = 'white', size=8, zorder=2000)
        try:
            midatl_flash_count
        except NameError:
            midatl_flash_count = 'Unavailable'
        clabeltext = 'Flash Count=' + str(midatl_flash_count)
        fig.text(.16, .13,clabeltext,horizontalalignment='left', color = 'red', size=10, zorder=2000)

        ax.outline_patch.set_visible(False)
        ax.background_patch.set_visible(False)
        output_file = workdir + "ltng_mid/" + ABI_datetime[n] + ".png"
        fig.savefig(output_file, dpi=dpi, bbox_inches='tight', transparent=True)
        plt.close()

######################## TRUE COLOR WITH LIGHTNING GIFS ########################
######################## ######################## ######################## 
workdir = "/home/sat_ops/goesR/"
imgdir = "/home/sat_ops/goesR/ltng_conus/"
img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
img_names = sorted(img_list)[-72:]

imglen = len(img_names)
images = []
dur_vals = []
for i in range(1,imglen):
    if i != imglen:
        dur_vals.append(.07)
        
dur_vals.append(2)

for i in img_names:
    input_file=imgdir + str(i)
    images.append(imageio.imread(input_file))

imageio.mimsave(workdir + 'lightning_truecolor_conus.gif', images, duration=dur_vals)
imageio.mimsave(workdir + 'truecolor_conus.gif', images, duration=dur_vals)

# now for the midatlantic
imgdir = "/home/sat_ops/goesR/ltng_mid/"
img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
img_names = sorted(img_list)[-72:]

imglen = len(img_names)
images = []
dur_vals = []
for i in range(1,imglen):
    if i != imglen:
        dur_vals.append(.07)

dur_vals.append(2)

for i in img_names:
    input_file=imgdir + str(i)
    images.append(imageio.imread(input_file))
imageio.mimsave(workdir + 'lightning_truecolor_midatlantic.gif', images, duration=dur_vals)
imageio.mimsave(workdir + 'truecolor_midatlantic.gif', images, duration=dur_vals)