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
goes_main = xr.open_dataset(
    "http://basin.ceoe.udel.edu/thredds/dodsC/GOESJPLSST.nc")
goes_nc = goes_main.sel(time=slice("2020-04-21","2020-04-28"))



goes_nc = goes_nc.drop('projection')
goes_nc = goes_nc.drop('dt_analysis')
goes_nc = goes_nc.drop('satellite_zenith_angle')
goes_nc = goes_nc.drop('sses_standard_deviation')
goes_nc = goes_nc.drop('wind_speed')
goes_nc['sea_surface_temperature'] = goes_nc['sea_surface_temperature'] - goes_nc['sses_bias']
goes_nc = goes_nc.drop('sses_bias')
newtimestamp = goes_nc.time.values[-1]
for t in range(len(goes_nc.time.values)):
    x = goes_nc['sea_surface_temperature'][t]
    goes_nc['sea_surface_temperature'][t] = x.where(goes_nc['quality_level'][t] == 5)

goes_nc['sst'] = goes_nc['sea_surface_temperature']
goes_nc = goes_nc.drop(['sea_surface_temperature'])
goes_nc = goes_nc.drop(['quality_level'])
datasets = []
for i in range(2,11):
    datasets.append(goes_nc.sst[18*i])
combined = xr.concat(datasets, dim='time')

newtimestamp = combined.time.values[-1]
new8day = combined.mean('time')
newtimestamp = (newtimestamp - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
x = new8day.assign_coords(time=newtimestamp)
new8day = x.expand_dims('time')
new8day.time.attrs['units'] = 'seconds since 1970-01-01 00:00:00'
new8day.to_netcdf("~/Downloads/new8day.nc")





# run anomaly code on basin to generate new anomaly
# open the anomaly files

all8day = xr.open_dataset("Downloads/SSTanomalyGoesAqua8dayCelsius119.nc")
new = xr.open_dataset("Downloads/random8day0428.nc")
new.sst.values = new.sst.values - 273.15

diff = all8day.sst[0] - new.sst[0]


bounds=(-80.2,20.3,-70.85, 45.3)

diff = diff.sel(lat=slice(bounds[3], bounds[1]), lon=slice(bounds[0],bounds[2]))

diff.plot(size=15, vmin=-2, vmax=2, cmap='bwr')
plt.title("Anomaly Difference: (192 files in the 8 day composite) vs. (8 14Z files in the 8 day composite")
