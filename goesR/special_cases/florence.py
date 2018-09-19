# this script is for making the Florence Images into a gif
# This is a quick patch for the mesoscale scans
# we need to do true color and band 13 here
# first let's just test it out and then afterwards we can 
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
from matplotlib import patheffects, ticker
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.image as image
from matplotlib.colors import LinearSegmentedColormap, ListedColormap # Linear interpolation for color maps

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
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from dateutil import tz
import time
from time import mktime
import os.path
from datetime import datetime, timedelta
from os import listdir
from os.path import isfile, join
import calendar
from cpt_convert import loadCPT # Import the CPT convert function
import imageio
from PIL import Image
#suppress deprecation warnings
import warnings
warnings.simplefilter("ignore", category=DeprecationWarning)
############# Initial Set Up ##################
datadir = "/home/sat_ops/goesR/special_cases/florence/data/"
workdir = "/home/sat_ops/goesR/special_cases/florence/"

def contrast_correction(color, contrast):
    F = (259*(contrast + 255))/(255.*259-contrast)
    COLOR = F*(color-.5)+.5
    COLOR = np.minimum(COLOR, 1)
    COLOR = np.maximum(COLOR, 0)
    return COLOR


file_names = [f for f in listdir(datadir) if isfile(join(datadir, f))]
fnamelist = sorted(file_names)

for i in range(0,len(fnamelist)):
    if os.path.isfile("/home/sat_ops/goesR/special_cases/florence/tcimg/" + fnamelist[i] + ".png") == False:
        Cnight = Dataset(datadir + fnamelist[i], 'r')
        Cnight2 = xr.open_dataset(datadir + fnamelist[i])
        dat = Cnight2.metpy.parse_cf("CMI_C01")
        proj = dat.metpy.cartopy_crs
        newproj = ccrs.Mercator()
        ##### Need to do this by band because I don't think the mcmipc exists on thredds
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
        b13 = cleanIR
        cleanIR = np.maximum(cleanIR, 90)
        cleanIR = np.minimum(cleanIR, 313)
        cleanIR = (cleanIR-90)/(313-90)
        cleanIR = 1 - cleanIR
        cleanIR = cleanIR/1.5
    
        contrast = 125
        RGB_contrast = contrast_correction(np.dstack([R, G_true, B]), contrast)
        RGB_contrast_IR = np.dstack([np.maximum(RGB_contrast[:,:,0], cleanIR), np.maximum(RGB_contrast[:,:,1], cleanIR), np.maximum(RGB_contrast[:,:,2], cleanIR)])
        # Create a color tuple for pcolormesh
        rgb = RGB_contrast_IR[:,:-1,:] # Using one less column is very imporant, else your image will be scrambled! (This is the stange nature of pcolormesh)
        rgb = np.minimum(rgb, 1)
        colorTuple = rgb.reshape((rgb.shape[0] * rgb.shape[1]),3) # flatten array, becuase that's what pcolormesh wants.
        colorTuple = np.insert(colorTuple, 3, 1.0, axis=1) # adding an alpha channel will plot faster?? according to stackoverflow.
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
        
        title = 'Hurricane Florence - GOES True Color - Powered By CEMA'
        timestr = local.strftime('%Y-%m-%d %H:%M ') + et
        
        fig.text(toptextright, toptext,timestr,horizontalalignment='left', color = 'white', size=14, zorder=2000)
        fig.text(bottomtextleft, toptext,title,horizontalalignment='left', color = 'white', size=14, zorder=2000)
        
        #request = cimgt.GoogleTiles(url="https://cartodb-basemaps-d.global.ssl.fastly.net/dark_nolabels/{z}/{x}/{y}.png")
        #ax.add_image(request, 7, zorder=0, interpolation='none')
        im = Image.open("/home/sat_ops/goesR/zfolder/combinedlarge.png")
        im = np.array(im).astype(np.float) / 255
        plt.figimage(im, 15, 15, zorder=1)
        output_file = workdir + "tcimg/" + fnamelist[i] + ".png"
        fig.savefig(output_file, dpi=dpi, bbox_inches='tight')
        plt.close()
    
        print('plotting' + fnamelist[i])
        
        
        colorscheme = 'IR4AVHRR6.cpt'
        v_min = -103
        v_max = 84
        clabeltext = "Brightness Temperature [DegC]"
        cpt = loadCPT('/home/sat_ops/goesR/indbands/colortables/' + colorscheme)
        cpt_convert = LinearSegmentedColormap('cpt', cpt) 
        state_col = 'black'
        b13 = b13 - 273.15
        
        fig = plt.figure(figsize=[fs_x, fs_y], dpi=dpi)
        ax = fig.add_subplot(1,1,1, projection=newproj)
        im = ax.pcolormesh(dat['x'], dat['y'], b13, cmap=cpt_convert,vmin = v_min, vmax=v_max, transform=proj)
        ax.set_extent((-65, -128, 21, 47))
        ax.set_title("")
    
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                        edgecolor='darkslategray', facecolor='none',linewidth=1))
        ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                        edgecolor='darkslategray', facecolor='none',linewidth=1))
        ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                        edgecolor='darkslategray', facecolor='none',linewidth=1))
        # top rectangle
        fig.patches.extend([plt.Rectangle((toprecx,toprecy),0.7745,0.025,
                                      fill=True, color='black', alpha=1, zorder=1000,
                                      transform=fig.transFigure, figure=fig)])
                                      
        title = 'Hurricane Florence - GOES Band 13 - Powered By CEMA'
        timestr = local.strftime('%Y-%m-%d %H:%M ') + et
        
        fig.text(toptextright, toptext,timestr,horizontalalignment='left', color = 'white', size=14, zorder=2000)
        fig.text(bottomtextleft, toptext,title,horizontalalignment='left', color = 'white', size=14, zorder=2000)
        #request = cimgt.GoogleTiles(url="https://cartodb-basemaps-d.global.ssl.fastly.net/dark_nolabels/{z}/{x}/{y}.png")
        #ax.add_image(request, 7, zorder=0, interpolation='none')
        im = Image.open("/home/sat_ops/goesR/zfolder/combinedlarge.png")
        # We need a float array between 0-1, rather than
        # a uint8 array between 0-255
        im = np.array(im).astype(np.float) / 255
        plt.figimage(im, 15, 15, zorder=1)
        output_file = workdir + "b13img/" + fnamelist[i] + ".png"
        fig.savefig(output_file, dpi=dpi, bbox_inches='tight')
        plt.close()


giflist = ['tcimg']
gifnames = ['florence_truecolor']

for g in range(0,len(giflist)):
    imgdir = workdir + giflist[g] + '/'
    
    img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
    img_names = sorted(img_list)
    
    imglen = len(img_names)
    images = []
    dur_vals = []
    for i in range(0,imglen -1):
        if i != imglen:
            dur_vals.append(.05)
    dur_vals.append(2)
   
    writer = imageio.get_writer('test.mp4', fps=20) 
    for i in img_names:
        print(i)
        input_file=imgdir + str(i)
        writer.append_data(imageio.imread(input_file, format="PNG"))
    writer.close()

    
    
