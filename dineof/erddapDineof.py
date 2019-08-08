# prepare goes16 sst data for dineof
# dineof removes lat/lon values & time demensions, so we need to add those back using the dimensions
# of the original file
import xarray as xr
import numpy as np
import metpy
from datetime import datetime, timedelta
import pandas as pd

for a in range(1,6):
    dineof_nc = xr.open_dataset("/data/dineof_in/output/dineofOUT_roffs_area" + str(a) + "_7day.nc")
    roffs_nc = xr.open_dataset("/home/james/roffs/roffs_area" + str(a) + "_7day.nc")
    roffs_nc = roffs_nc.isel(time=-1)
    oldtimestamp = roffs_nc.time.values
    roffs_nc.sst.values = dineof_nc.forecasted_sst[-1].values
    roffs_nc.sst.values = (roffs_nc.sst.values - 273.15)*(9/5) + 32
    roffs_nc['forecasted_sst'] = roffs_nc['sst']
    roffs_nc = roffs_nc.drop(['sst'])
    roffs_nc = roffs_nc.where(roffs_nc.forecasted_sst > 31)
    
    roffs_nc.attrs['units'] = 'Fahrenheit'
    newtimestamp = (oldtimestamp - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
    x = roffs_nc.assign_coords(time=newtimestamp)
    roffs_nc = x.expand_dims('time')
    roffs_nc.time.attrs['units'] = 'seconds since 1970-01-01 00:00:00'
    outpath = "/home/james/roffs/dineof_out/"
    
    # Create the Timestamp object
    pdTime = pd.to_datetime(oldtimestamp)
    roffs_nc.to_netcdf(path=outpath + 'roffs_area' + str(a) + '_' + str(pdTime.year) + str(pdTime.month).zfill(2)
                                    + str(pdTime.day).zfill(2) + str(pdTime.hour).zfill(2) + str(pdTime.minute).zfill(2)
                                    + '.nc',format='NETCDF3_CLASSIC')