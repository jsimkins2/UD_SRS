import xarray as xr
import numpy as np

# define area of interest
lat_bounds = [18.0, 31.0]
lon_bounds = [-98, -80]

# grab sst data from the last 3 days and use DQF == 0 data
goes_nc = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/goes_r_sst.nc")
goes_nc = goes_nc.sel(latitude=slice(18.0,31.0), longitude=slice(-98,-80), 
          time=slice('2018-10-30'))

goes_nc = goes_nc.resample(time='1D').mean('time')
#aggregate daily sst
sst1 = goes_nc['SST'][-1]
dqf1 = goes_nc['DQF'][-1]

sst2 = goes_nc['SST'][-2]
dqf2 = goes_nc['DQF'][-2]

sst1 = sst1.where(dqf1 == 0)
sst2 = sst2.where(dqf2 == 0)
# concatenate these two xarrays to the same dataset
sst1 = sst1.rename('sst')
sst2 = sst2.rename('sst')
sst3 = xr.concat([sst2, sst1], dim='time')

# try above but with resample
x = sst3.resample(time='1D').mean('time')
# save as uncompressed netcdf level 3