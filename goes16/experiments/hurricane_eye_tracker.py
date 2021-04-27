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
#from cpt_convert import loadCPT # Import the CPT convert function
import imageio
from PIL import Image
#suppress deprecation warnings
import warnings
warnings.simplefilter("ignore", category=DeprecationWarning)
############# Initial Set Up ##################
datadir = "/home/sat_ops/goesR/data/meso/"
workdir = "/home/sat_ops/goesR/"

def contrast_correction(color, contrast):
    F = (259*(contrast + 255))/(255.*259-contrast)
    COLOR = F*(color-.5)+.5
    COLOR = np.minimum(COLOR, 1)
    COLOR = np.maximum(COLOR, 0)
    return COLOR


#Cnight = Dataset("/Users/leathers/Documents/Delaware/goesR/data/OR_ABI-L2-MCMIPC-M3_G16_s20181711212251_e20181711215024_c20181711215133.nc", 'r')
#Cnight2 = xr.open_dataset("/Users/leathers/Documents/Delaware/goesR/data/OR_ABI-L2-MCMIPC-M3_G16_s20181711212251_e20181711215024_c20181711215133.nc")
filename = "Downloads/GOES16_FullDisk_20180914_001530_10.35_6km_0.0S_75.0W.nc4"
filename = "Downloads/OR_ABI-L2-DSIF-M3_G16_s20182570000305_e20182570011072_c20182570011400.nc"
Cnight = Dataset(filename, 'r')
Cnight2 = xr.open_dataset(filename)
dat = Cnight2.metpy.parse_cf("CAPE")
proj = dat.metpy.cartopy_crs
newproj = ccrs.Mercator()
plat = ccrs.PlateCarree()
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

fig = plt.figure(figsize=[12,8])
ax = fig.add_subplot(1,1,1, projection=plat)
ax.pcolormesh(dat['x'], dat['y'],dat, transform=proj)
ax.set_extent((-10, -105, -10, 50))
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
              linewidth=0.5, color='darkgray', alpha=0.5, linestyle='--')


from scipy import ndimage

test = dat.values
test = np.ma.masked_invalid(test)
x, y = ndimage.measurements.center_of_mass(test)

from pyproj import Proj, transform
import pyproj
p = Proj(proj='geos', h=Cnight2.attrs['satellite_altitude'], lon_0=Cnight2.attrs['satellite_longitude'], sweep='x')
wgs84=Proj("+init=EPSG:3395")

transform(p, wgs84, x, y)

newproj.transform_point(x,y,src_crs=proj)
t,r = plat.transform_point(x,y,src_crs=proj)

_x, _y = newproj.transform_point(t,r,src_crs=plat)