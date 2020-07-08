# this script is designed to test band 1 - blue band and band 3 - veggie band - as proxies for water / chlorophyll/oceancolor

import xarray as xr
# prepare goes16 sst data for dineof
from datetime import datetime, timedelta
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from matplotlib import patheffects
import metpy
from metpy.plots import colortables
import xarray as xr
from xarray.backends import NetCDF4DataStore
from netCDF4 import Dataset, num2date
from dateutil import tz
import time
from time import mktime
import os.path
import os
from os import listdir
from os.path import isfile, join
import cartopy.crs as ccrs
import pandas as pd
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh
import matplotlib.font_manager as font_manager
from scipy.ndimage.filters import gaussian_filter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib as mpl
import random
mpl.rcParams["contour.negative_linestyle"] = 'solid'



# load our area, in this case area 2
area = [36.5, 40, -75.25, -72.0]
area = [26.5, 29.5, -91.0, -86.5]
extArea = [area[2], area[3], area[0], area[1]]

# load SST
today = datetime.utcnow() - timedelta(days=7)
goes_thredds = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/GOESR_SST_DAILY.nc")
goes_thredds = goes_thredds.sel(latitude=slice(area[1],area[0]), longitude=slice(area[2], area[3]), time=datetime.strftime(today, '%Y-%m-%d'))
sst = goes_thredds.metpy.parse_cf('sst')[-1]
sst = sst.where(sst.values > 200)
sst = sst - 273.15
sst = sst*(9/5) + 32
sstproj = sst.metpy.cartopy_crs


viirs = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/viirs_1day_aggregate.nc')
viirs = viirs.sel(lat=slice(area[0],area[1]), lon=slice(area[2], area[3]), time="2012-04-06")

aqua = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/Aqua1DayAggregate.nc')
aqua = aqua.sel(lat=slice(area[0],area[1]), lon=slice(area[2], area[3]), time="2012-04-06")













# load chlorophyll
viirs = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/viirs_1day_aggregate.nc')
viirs = viirs.sel(lat=slice(area[0],area[1]), lon=slice(area[2], area[3]), time=datetime.strftime(today, '%Y-%m-%d'))
r486 = viirs.metpy.parse_cf('Rrs_486')[0]
r486 = r486.where(r486.values > 0) #blue
r862 = viirs.metpy.parse_cf('Rrs_862')[0] #green
r862 = r862.where(r862.values > 0) 
viirsproj = ccrs.PlateCarree()