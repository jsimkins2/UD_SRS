import xarray as xr
import numpy as np

# define area of interest
25.777711167924" N, 80.1865375041962" W
lat_bounds = [23, 27]
lon_bounds = [-80, -76]
# Biscayne Bay 
# grab sst data from the last 3 days and use DQF == 0 data
goes_nc = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/goes_r_sst.nc")
goes_nc = goes_nc.sel(latitude=slice(23, 27), longitude=slice(-80, -76), 
          time=slice('2019-02-27', '2019-02-28'))

#dineof_ssts = {}
#sst_list = [str('sst' + str(f)) for f in range(len(goes_nc.time.values))]
#dineof_ssts = {sst_list : goes_nc.time.values}
#dineof_ssts = dict(zip(sst_list,goes_nc.time.values))

dineof_ssts = {}
for t in range(len(goes_nc.time.values)):
    tem_nc = goes_nc.sel(time=str(goes_nc.time.values[t]))
    tem_sst = tem_nc['SST'].where(tem_nc['DQF'] == 0)
    tem_sst = tem_sst.rename('SST_qc')
    dineof_ssts[t] = tem_sst

for t in range(len(dineof_ssts)):
    x = dineof_ssts[t]
    x.to_netcdf(path='Downloads/' + 'biscayne_SST_qc_' + str(str(goes_nc.time.values[t]).split('.')[0]), format = 'NETCDF3_CLASSIC')

xr.DataArray.to_netcdf(dineof_ssts[0])

x.to_netcdf
goes_nc = goes_nc.resample(time='1D').mean('time')

dineof_dataset0

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