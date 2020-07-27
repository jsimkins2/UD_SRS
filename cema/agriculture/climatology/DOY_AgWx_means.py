# deos geopandas to be run at home 
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.colors import ListedColormap
import xarray as xr
import matplotlib as mpl
import pandas as pd
import time
from datetime import datetime, timedelta, date
from netCDF4 import num2date
import matplotlib.pyplot as plt

    
# load in the datasets 
bounds=(-76.2,38.3,-74.85, 40.3)

# open agwx main dataset
agwx_main = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/DEOSAG.nc")
agwx_main = agwx_main.sel(latitude=slice(bounds[3], bounds[1]), longitude=slice(bounds[0],bounds[2]))

# open NCEP stage IV QC dataset
dsPrec = xr.open_dataset("http://thredds.demac.udel.edu/thredds/dodsC/NCEPIVQC.nc")
dsPrec = dsPrec.sel(lat=slice(bounds[1], bounds[3]), 
                    lon=slice(bounds[0],bounds[2]), 
                    time=slice(datetime.strptime("2010-01-01", "%Y-%m-%d"),
                              date.today()))

dsPrec = dsPrec.reindex(lat=list(reversed(dsPrec.lat)))
dsPrec = dsPrec.rename(name_dict= {'lat' : 'latitude'})
dsPrec = dsPrec.rename(name_dict= {'lon' : 'longitude'})
dsPrec = dsPrec.drop('crs')

# calculate 
agwx_doy = agwx_main.groupby("time.dayofyear").mean("time")
prec_doy = dsPrec.groupby("time.dayofyear").mean("time")

# add precipitation flux to the entire dataset
agwx_doy['NCEPstageIVPrecip']=(['dayofyear', 'latitude', 'longitude'],  prec_doy['Precipitation_Flux'])

agwx_doy.to_netcdf("/data/DEOS/doy_climatology/deos_doy_climatology.nc")


ncc = xr.open_dataset('')

