# This is a quick patch for the fulldisk atlantic basin
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
datadir = "/home/sat_ops/goesR/special_cases/dorian/data/"
ltngdir = "/home/sat_ops/goesR/special_cases/dorian/ltng_data/"
workdir = "/home/sat_ops/goesR/special_cases/dorian/"

def contrast_correction(color, contrast):
    F = (259*(contrast + 255))/(255.*259-contrast)
    COLOR = F*(color-.5)+.5
    COLOR = np.minimum(COLOR, 1)
    COLOR = np.maximum(COLOR, 0)
    return COLOR


abi_names = [f for f in listdir(datadir) if isfile(join(datadir, f))]
ltng_names = [f for f in listdir(ltngdir) if isfile(join(ltngdir, f))]
abi_names = sorted(abi_names)
ltng_names = sorted(ltng_names)
ltng_names = ltng_names[1:]

if len(abi_names) > 0:
    for n in range(0, len(abi_names)):
        if os.path.isfile(workdir + 'tcimg/' + abi_names[n][:-3] + ".png") == False:
            t = abi_names[n][:-3]
            print(t)

            Cnight = Dataset(datadir + abi_names[n], 'r')
            Cnight2 = xr.open_dataset(datadir + abi_names[n])
            dat = Cnight2.metpy.parse_cf("CMI_C01")
            proj = dat.metpy.cartopy_crs
            newproj = ccrs.PlateCarree()
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
            
                    ########## LIGHTNING ##############
            # Now that we have all the satellite ABIs read in, we need to call the lightning files 
            # create dictionaries for the lat/lons
            ltng_lat = {}
            ltng_lon = {}
            
            ltfile = ltng_names[n]
            L_file = ltngdir + str(ltfile) # GOES16 East
            L = Dataset(L_file, 'r')
            ltng_lat = L.variables['flash_lat'][:]
            ltng_lon = L.variables['flash_lon'][:]

            ltngxr = xr.open_dataset(L_file)
            ltngxr = ltngxr.metpy.parse_cf("flash_energy")
            ltngproj = ltngxr.metpy.cartopy_crs

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
            symbol = u'$\u26A1$'
            
        
            fig = plt.figure(figsize=[16,9], dpi=100)
            ax = fig.add_subplot(1,1,1, projection=newproj)
            im = ax.pcolormesh(dat['x'], dat['y'], R, color=colorTuple, transform=proj)
            for g in range(0, len(ltng_lat)):
                ax.scatter(ltng_lon[g], ltng_lat[g], s=18, marker=symbol, c='red', edgecolor='red', lw=0, transform=ltngproj)

            ax.set_extent((-45, -105, 0, 45))
            timestr = local.strftime('%Y-%m-%d %H:%M ') + et
            fig.text(0.5,0.9, 'Hurricane Dorian TC + Ltng - Powered By CEMA\n'
                    + timestr,horizontalalignment='center',fontsize=16)
            ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                            edgecolor='black', facecolor='none',linewidth=0.3))
            ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                            edgecolor='black', facecolor='none',linewidth=0.3))
            ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                            edgecolor='black', facecolor='none',linewidth=0.3))
            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.5, color='darkgray', alpha=0.5, linestyle='--')
            gl.xlabels_top = False
            gl.ylabels_left = False
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            gl.xlabel_style = {'color': 'black', 'size':'small'}
            gl.ylabel_style = {'color': 'black', 'size':'small'}
        
            request = cimgt.GoogleTiles(url="https://cartodb-basemaps-d.global.ssl.fastly.net/dark_nolabels/{z}/{x}/{y}.png")
            ax.add_image(request, 7, zorder=0, interpolation='none')
            im = Image.open("/home/sat_ops/goesR/zfolder/combinedsmall.png")
            # We need a float array between 0-1, rather than
            # a uint8 array between 0-255
            im = np.array(im).astype(np.float) / 255
            plt.figimage(im,15, 30, zorder=1, alpha=0.8)
            output_file = workdir + "tcimg/" +  abi_names[0][:-3] + ".png"
            fig.savefig(output_file, dpi=dpi, bbox_inches='tight')
            plt.close()
            
            # NOW PLOT THE BAND 13 STUFF
            colorscheme = 'IR4AVHRR6.cpt'
            v_min = -103
            v_max = 84
            cpt = loadCPT('/home/sat_ops/goesR/indbands/colortables/' + colorscheme)
            cpt_convert = LinearSegmentedColormap('cpt', cpt) 
            state_col = 'black'
            b13 = b13 - 273.15
            
            fig = plt.figure(figsize=[16,9], dpi=100)
            ax = fig.add_subplot(1,1,1, projection=newproj)
            im = ax.pcolormesh(dat['x'], dat['y'], b13, cmap=cpt_convert,vmin = v_min, vmax=v_max, transform=proj)
            ax.set_extent((-45, -105, 0, 45))
            timestr = local.strftime('%Y-%m-%d %H:%M ') + et
            fig.text(0.5,0.9, 'Hurricane Dorian Band13 - Powered By CEMA\n'
                    + timestr,horizontalalignment='center',fontsize=16)
            ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                            edgecolor='black', facecolor='none',linewidth=0.5))
            ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                            edgecolor='black', facecolor='none',linewidth=0.5))
            ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                            edgecolor='black', facecolor='none',linewidth=0.5))
            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.5, color='darkgray', alpha=0.5, linestyle='--')
            gl.xlabels_top = False
            gl.ylabels_left = False
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            gl.xlabel_style = {'color': 'black', 'size':'small'}
            gl.ylabel_style = {'color': 'black', 'size':'small'}
            
            request = cimgt.GoogleTiles(url="https://cartodb-basemaps-d.global.ssl.fastly.net/dark_nolabels/{z}/{x}/{y}.png")
            ax.add_image(request, 7, zorder=0, interpolation='none')
            im = Image.open("/home/sat_ops/goesR/zfolder/combinedsmall.png")
            # We need a float array between 0-1, rather than
            # a uint8 array between 0-255
            im = np.array(im).astype(np.float) / 255
            fig.figimage(im,15, 30, zorder=1, alpha=0.8)
        
            output_file = workdir + "b13img/" +  abi_names[0][:-3] + ".png"
            fig.savefig(output_file, dpi=dpi, bbox_inches='tight')
            plt.close()
        
        
giflist = ["tcimg", "b13img"]
gifnames = ['dorian_truecolor', 'dorian_band13']

for g in range(0,len(giflist)):
    imgdir = workdir + giflist[g] + '/'
    
    img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
    img_names = sorted(img_list)
    
    imglen = len(img_names)
    images = []
    dur_vals = []
    for i in range(0,imglen -1):
        if i != imglen:
            dur_vals.append(.07)
    dur_vals.append(2)
    
    writer = imageio.get_writer(workdir + 'dorian.mp4', fps=20) 
    for i in img_names:
        print(i)
        input_file=imgdir + str(i)
        writer.append_data(imageio.imread(input_file, format="PNG"))
    writer.close()
    
    
    
