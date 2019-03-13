# goes16 quality control algorithm
lat_bounds = [23, 27]
lon_bounds = [-81, -77]
# Biscayne Bay 
# grab sst data from the last 3 days and use DQF == 0 data
goes_nc = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/goes_r_sst.nc")
#goes_nc = goes_nc.sel(longitude=[-90], method='nearest')
#sst = sst.sel(longitude=[-51], method='nearest')


import xarray as xr
import metpy
import datetime

badfiles2 = []
ds = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/goes_r_sst.nc")
for t in range(0,len(ds.time.values)):
    sst = ds.metpy.parse_cf("SST")[t]
    #sst = sst.sel(time=datetime.datetime.strptime('2019-03-08 16:28:18', '%Y-%m-%d %H:%M:%S'))
    sst = sst.sel(longitude=slice(-65, -50), latitude=slice(16,20))
    sst = sst.where(sst.values > 270, np.nan)
    sst = sst.fillna(-999)
    if np.percentile(sst.values, 60) == -999:
        badfiles2.append(sst.time.values)
    





#
from itertools import groupby
[len(list(v)) for k,v in groupby(sst.values) if k==-999]


