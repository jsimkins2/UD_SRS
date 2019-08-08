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

from dateutil import tz
import time
from time import mktime
import os.path
from datetime import datetime, timedelta
from os import listdir
from os.path import isfile, join
import calendar
 
#suppress deprecation warnings
import warnings
warnings.simplefilter("ignore", category=DeprecationWarning)

fs_x = 16
fs_y = 12
dpi = 100
toptext = 0.8065
toptextleft = 0.13
toptextright = 0.8
bottomtextleft = 0.13
bottomtextheight = 0.212
toprecx = 0.125
toprecy = 0.799
bottomrecx = 0.125
bottomrecy = 0.19

newproj = ccrs.Mercator()
# Now plot the conus Lightning
fig = plt.figure(figsize=[fs_x, fs_y], dpi=dpi)
ax = fig.add_subplot(1,1,1, projection=newproj)
ax.set_extent((-65, -128, 21, 47))
request = cimgt.GoogleTiles(url="https://cartodb-basemaps-d.global.ssl.fastly.net/dark_nolabels/{z}/{x}/{y}.png")
ax.add_image(request, 7, zorder=0, interpolation='none')
fig.patches.extend([plt.Rectangle((bottomrecx,.18),0.7745,0.04,
                              fill=True, color='blue', alpha=1, zorder=1000,
                              transform=fig.transFigure, figure=fig)])
clabeltext = 'Flash Count=' 
fig.text(.8, 0.225,clabeltext,horizontalalignment='left', color = 'white', size=10, zorder=2000)
im1 = image.imread("/home/sat_ops/goesR/zfolder/udelcemagoes24.png")
plt.figimage(im1, 15, 65, zorder=1)
output_file = "/home/sat_ops/web/images/background/conus.png"
fig.savefig(output_file, dpi=dpi, bbox_inches='tight', transparent=True)

