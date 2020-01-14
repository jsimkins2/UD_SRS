import pandas as pd
import sys
from os.path import isfile, join
from os import listdir
import os
import os.path
from time import mktime
import time
from dateutil import tz
import cartopy.io.img_tiles as cimgt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from netCDF4 import Dataset, num2date
from xarray.backends import NetCDF4DataStore
from metpy.plots import colortables
from matplotlib import patheffects
import xarray as xr
from netCDF4 import Dataset
import numpy as np
import metpy
from datetime import datetime
import cartopy.feature as cfeature
import matplotlib.pyplot as plt

# now make a rolling 1 day aka last 24 hours
goes_nc = xr.open_dataset(
    "http://basin.ceoe.udel.edu/thredds/dodsC/GOESR_SST_DAILY.nc")
# grab the last 24 hours of sst dataset
goes_nc = goes_nc.isel(time=range(-8, 0))
newtimestamp = goes_nc.time.values[-1]

goes_nc = goes_nc.mean('time')
newtimestamp = (newtimestamp - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
x = goes_nc.assign_coords(time=newtimestamp)
goes_nc = x.expand_dims('time')
goes_nc.time.attrs['units'] = 'seconds since 1970-01-01 00:00:00'
outpath = "/home/sat_ops/goesR/sstClimatology/"
goes_nc.to_netcdf(path=outpath + 'GOES16_SST_8day.nc')
