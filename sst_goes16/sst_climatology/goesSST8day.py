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
from datetime import datetime, timedelta, date
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import pandas as pd

start_date = date(2018, 1, 1)
end_date = datetime.utcnow()
daterange = pd.date_range(start_date, end_date)

workdir = '/data/aquaGoesSST/C'
# now make a rolling 1 day aka last 24 hours
goes_main = xr.open_dataset(
    "http://basin.ceoe.udel.edu/thredds/dodsC/GOESR_SST_DAILY.nc")

for d in daterange:
    # grab the last 24 hours of sst dataset
    goes_nc = goes_main.sel(time=slice(d.to_pydatetime() - timedelta(days=8), d.to_pydatetime()))
    newtimestamp = goes_nc.time.values[-1]
    newtimestamp = (newtimestamp - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
    
    if os.path.isfile("/data/aquaGoesSST/C/" + str(time.strftime('%Y', time.localtime(newtimestamp))) + "/SSTanomalyGoesAqua8dayCelsius" + str(time.strftime('%Y%m%d', time.localtime(newtimestamp))) + ".nc") == False:
        goes_nc = goes_nc.mean('time')
        print(str(time.strftime('%Y%m%d', time.localtime(newtimestamp))))
        x = goes_nc.assign_coords(time=newtimestamp)
        goes_nc = x.expand_dims('time')
        goes_nc.time.attrs['units'] = 'seconds since 1970-01-01 00:00:00'
        outpath = "/home/sat_ops/goesR/jpl_sst/sstClimatology/"
        goes_nc.to_netcdf(path=outpath + 'GOES16_SST_8day.nc')
        os.system("Rscript /home/sat_ops/goesR/github/UD_SRS/sst_goes16/sst_climatology/goesVsAquaClimatologyCelsius.R")
        os.system("Rscript /home/sat_ops/goesR/github/UD_SRS/sst_goes16/sst_climatology/goesVsAquaClimatologyCelsius.R")
        os.system("rm /home/sat_ops/goesR/jpl_sst/sstClimatology/GOES16_SST_8day.nc")

