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


import xarray as xr
import numpy as np
import metpy
from datetime import datetime, timedelta
import pandas as pd
# paths
outpath1 = "/data/GOES/GOES-R/1day/"
outpath2 = "/data/GOES/GOES-R/daily_composite/"
datelist = pd.date_range('2019-06-20', pd.datetime.today()).tolist()

ncep_nc = xr.open_dataset("http://thredds.demac.udel.edu/thredds/dodsC/ncep_stage_iv.nc?precipitation_flux[0:10]")
goes_nc = goes_nc.sel(time=datetime.strftime(
    datelist[d].date(), '%Y-%m-%d'))
goes_nc = goes_nc.drop('Band15')
newtimestamp = goes_nc.time.values[-1]
    
    
    