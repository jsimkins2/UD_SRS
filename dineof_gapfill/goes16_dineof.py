import xarray as xr
import numpy as np

# define area of interest
#25.777711167924" N, 80.1865375041962" W
lat_bounds = [23, 27]
lon_bounds = [-81, -77]
# Biscayne Bay 
# grab sst data from the last 3 days and use DQF == 0 data
goes_nc = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/goes_r_sst.nc")
goes_nc = goes_nc.sel(latitude=slice(lat_bounds[0],lat_bounds[1]), longitude=slice(lon_bounds[0], lon_bounds[1]), 
          time=slice('2019-02-27', '2019-02-28'))

for t in range(len(goes_nc.time.values)):
    x = goes_nc['SST'][t]
    goes_nc['SST'][t] = x.where(goes_nc['DQF'][t] == 0)
    x = goes_nc['DQF'][t]
    goes_nc['DQF'][t] = x.where(goes_nc['DQF'][t] == 3)

#landmask = goes_nc['DQF'][0]
goes_nc = goes_nc.drop(['DQF'])
goes_nc = goes_nc.drop('Band15')
goes_nc['sst'] = goes_nc['SST']
goes_nc = goes_nc.drop(['SST'])


forecast_nc = goes_nc.isel(time=[-3,-2,-1])
forecast_nc.time.values = forecast_nc.time.values.astype('datetime64[s]') + (3600*3)
forecast_nc['sst'] = forecast_nc['sst'].where(forecast_nc['sst'] < 2)
goes_nc= xr.concat([forecast_nc, goes_nc], dim='time')
#goes_nc = goes_nc.fillna(-999)
#landmask.to_netcdf()
goes_nc.to_netcdf(path='Downloads/goes16_sst_biscayne_bay_02272019_02282019.nc', format='NETCDF3_CLASSIC')

'''
sst3 = xr.concat([sst2, sst1], dim='time')

# try above but with resample
x = sst3.resample(time='1D').mean('time')
# save as uncompressed netcdf level 3