import xarray as xr
from datetime import datetime, timedelta, date

# write the most recent year to local disk
ds_write = xr.open_dataset("http://thredds.demac.udel.edu/thredds/dodsC/NCEPIVQC.nc").sel(time=slice(datetime.strptime(str(str(date.today().year) + "-01-01"), "%Y-%m-%d"),date.today()))
ds_write.to_netcdf(str("/home/james/agriculture/ncep_stageIV/aggregate_quality/ncep_stageIV_quality24hr_" + str(date.today().year) + ".nc"), mode='w')
ds_write.close()
print('dont with ncep')

# Declare bounds of the data
bounds=(-76.2,38.3,-74.85, 40.3)

# read in the refET dataset
dsRefET = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/DEOSAG.nc")
dsRefET = dsRefET.sel(latitude=slice(bounds[3], bounds[1]), longitude=slice(bounds[0],bounds[2]))
dsRefET.to_netcdf("/home/james/agriculture/deos_data/dsRefET.nc", mode='w')

print('done with refet')

climo = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/deos_doy_climatology.nc")
climo = climo.sel(latitude=slice(bounds[3], bounds[1]), longitude=slice(bounds[0],bounds[2]))
climo.to_netcdf("/home/james/agriculture/deos_data/climo.nc", mode='w')

print('done with climo')

