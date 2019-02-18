# Precipitation Depiction 
# By James Simkins with assistance from Dan Moore
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
from matplotlib import patheffects, ticker
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.image as image
from matplotlib.colors import LinearSegmentedColormap

from siphon.catalog import TDSCatalog
from siphon.radarserver import RadarServer
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
from scipy.interpolate import griddata
import numpy.ma as ma
from pyproj import Proj, transform
from metpy.plots import USCOUNTIES

from dateutil import tz
import time
from time import mktime
import os.path
from datetime import datetime, timedelta
from os import listdir
from os.path import isfile, join
import os.path
import os
import sys
from PIL import Image

import calendar
import imageio

# http://rammb.cira.colostate.edu/training/visit/quick_guides/Day_Cloud_Phase_Distinction.pdf
############# Initial Set Up ##################
workdir = "/home/sat_ops/goesR/daycloud/"
datadir = "/home/sat_ops/goesR/data/mcmipc/"
imgdir = "/home/sat_ops/goesR/daycloud/img_conus/"


# create a logfile with most recent files created in the shell script
file_names = [f for f in listdir(datadir) if isfile(join(datadir, f))]

fnamelist = []

fnamelist = sorted(file_names)[-3:]

ABI_datetime = []
for i in fnamelist:    
    if os.path.isfile(imgdir + str(i) + ".png") == False:
        ABI_datetime.append(i)

if len(ABI_datetime) > 0:
    for n in range(0, len(ABI_datetime)):
        
        mcmipc_file = xr.open_dataset(datadir + str(ABI_datetime[n]))
        R = mcmipc_file.metpy.parse_cf('CMI_C13')
        proj = R.metpy.cartopy_crs
        G = mcmipc_file.metpy.parse_cf('CMI_C02')
        B = mcmipc_file.metpy.parse_cf('CMI_C05')

        # Apply range limits for each channel becuase RGB values must be between 0 and 1
        R = np.maximum(R, (273.15-53.5))
        R = np.minimum(R, (7.5+273.15))
        R = (R - (273.15-53.5))/(((273.15+7.5) - (273.15-53.5)))
        R = 1 - R 
        #cleanIR = (cleanIR-90)/(313-90)
        
        G = np.maximum(G, 0)
        G = np.minimum(G, .78)
        B = np.maximum(B, .01)
        B = np.minimum(B, 0.59)
        
        gamma = 1
        R = np.power(R, gamma)
        G = np.power(G, gamma)
        B = np.power(B, gamma)
        
        
        RGB = np.dstack([R, G, B])
        rgb = RGB[:,:-1,:] # Using one less column is very imporant, else your image will be scrambled! (This is the stange nature of pcolormesh)
        rgb = np.minimum(rgb, 1)
        colorTuple = rgb.reshape((rgb.shape[0] * rgb.shape[1]),3) # flatten array, becuase that's what pcolormesh wants.
        colorTuple = np.insert(colorTuple, 3, 1.0, axis=1) # adding an alpha channel will plot faster?? according to stackoverflow.
        
        ymd = mcmipc_file.time_coverage_end.split("T")[0]
        hms = mcmipc_file.time_coverage_end.split("T")[1][:-3]
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
            et = "EDT"
        else:
            et = "EST"
        fs_x = 8
        fs_y = 8
        dpi = 100
        toptext = 0.8
        textleft = 0.137
        toptextright = 0.69
        bottomtextleft = 0.13
        bottomtextheight = 0.212
        toprecx = 0.125
        toprecy = 0.795
        bottomrecx = 0.125
        bottomrecy = 0.205
        
        fig = plt.figure(figsize=[16, 12], dpi=70)
        ax = fig.add_subplot(1,1,1, projection=ccrs.Mercator())
        im = ax.pcolormesh(R['x'], R['y'], R, color=colorTuple, transform=proj)
        ax.set_extent((-65, -128, 21, 47))
        timestr = local.strftime('%Y-%m-%d %H:%M ') + et
        fig.text(0.5,0.9, 'GOES16 Day Cloud Phase - Powered By CEMA\n'
                + timestr,horizontalalignment='center',fontsize=16)
                
        request = cimgt.GoogleTiles(url='http://3.api.cartocdn.com/base-midnight/{Z}/{X}/{Y}.png')
        
        ax.add_image(request, 7, zorder=0, interpolation='none')
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                        edgecolor='black', facecolor='none',linewidth=0.5))
        ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                        edgecolor='black', facecolor='none',linewidth=0.5))
        ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                        edgecolor='black', facecolor='none',linewidth=0.5))
        # add logo
        im = Image.open("/home/sat_ops/goesR/zfolder/combinedsmall.png")
        # We need a float array between 0-1, rather than
        # a uint8 array between 0-255
        im = np.array(im).astype(np.float) / 255
        plt.figimage(im,15, 30, zorder=1, alpha=0.8)
        
        # save file
        output_file = workdir + "img_conus" + ABI_datetime[n] + ".png"
        plt.savefig(output_file, dpi=100, bbox_inches='tight')
        plt.close()










######################## TRUE COLOR GIFS ########################
######################## ######################## ######################## 
import imageio
img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
img_names = sorted(img_list)[-72:]

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
imageio.mimsave(workdir + 'daycloud_phase.gif', images, duration=dur_vals)
