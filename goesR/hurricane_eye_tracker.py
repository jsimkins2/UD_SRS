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


Cnight = Dataset("/Users/james/Documents/Delaware/goesR/OR_ABI-L2-MCMIPM1-M3_G16_s20182561413206.nc", 'r')
Cnight2 = xr.open_dataset("/Users/james/Documents/Delaware/goesR/OR_ABI-L2-MCMIPM1-M3_G16_s20182561413206.nc")
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