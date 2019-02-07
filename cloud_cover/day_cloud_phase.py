import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
from siphon.catalog import TDSCatalog, get_latest_access_url
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
from matplotlib.patches import Rectangle
import pyart
from siphon.radarserver import RadarServer, get_radarserver_datasets
import numpy.ma as ma
import netCDF4
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import shapely.geometry as sgeom
import xarray
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.io.img_tiles as cimgt
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
import xarray as xr

# http://rammb.cira.colostate.edu/training/visit/quick_guides/Day_Cloud_Phase_Distinction.pdf

datadir="/Users/james/Downloads/"

mcmipc_file = xr.open_dataset(datadir + "ABI-L2-MCMIPC_2019_037_19_OR_ABI-L2-MCMIPC-M3_G16_s20190371902140_e20190371904513_c20190371905023.nc")
R = mcmipc_file.metpy.parse_cf('CMI_C13')
proj = R.metpy.cartopy_crs
G = mcmipc_file.metpy.parse_cf('CMI_C02')
B = mcmipc_file.metpy.parse_cf('CMI_C05')
# Turn empty values into nans


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

fig = plt.figure(figsize=[16, 12], dpi=70)
ax = fig.add_subplot(1,1,1, projection=ccrs.Mercator())
im = ax.pcolormesh(R['x'], R['y'], R, color=colorTuple, transform=proj)
ax.set_extent((-65, -128, 21, 47))
ax.set_title("")
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))








