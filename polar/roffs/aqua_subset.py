import xarray as xr
import numpy as np

# define area of interest
area = [18, 31, -98, -80]
# Biscayne Bay 
# grab sst data from the last 3 days and use DQF == 0 data
aqua_nc = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/Aqua1DayAggregate.nc")
aqua_nc = aqua_nc.sel(lat=slice(area[0],area[1]), lon=slice(area[2], area[3]))
aqua_nc = aqua_nc.sel(time=datetime.strptime('2019-03-13', '%Y-%m-%d'), method='nearest')

aqua_nc.to_netcdf(path='Downloads/aqua_test.nc', format='NETCDF3_CLASSIC')
