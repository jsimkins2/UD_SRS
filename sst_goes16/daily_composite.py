# prepare goes16 sst data for dineof
import xarray as xr
import numpy as np
import metpy
from datetime import datetime, timedelta
import pandas as pd
# paths
outpath = "/data/GOES/GOES-R/daily_composite/"

datelist = pd.date_range('2018-10-30', '2018-10-31')#,pd.datetime.today()).tolist()
for d in range(0,len(datelist)):
    goes_nc = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/goes_r_sst.nc")
    goes_nc = goes_nc.sel(time=datetime.strftime(datelist[d].date(), '%Y-%m-%d'))
    goes_nc = goes_nc.drop('Band15')
    for t in range(len(goes_nc.time.values)):
        x = goes_nc['SST'][t]
        goes_nc['SST'][t] = x.where(goes_nc['DQF'][t] == 0)
        x = goes_nc['DQF'][t]
        goes_nc['DQF'][t] = x.where(goes_nc['DQF'][t] == 3)
    
    landmask = goes_nc['DQF'][0]
    landmask = landmask.where(landmask.values == 3, 1)
    landmask = landmask.where(landmask.values == 1, 0)
    goes_nc = goes_nc.drop(['DQF'])
    goes_nc['sst'] = goes_nc['SST']
    goes_nc = goes_nc.drop(['SST'])
    # resample to a daily composite
    goes_nc = goes_nc.resample(time='1D').mean('time')
    #landmask.to_netcdf(path=outpath + 'landmask_roffs_' +  'area' + str(a) + '_' + str(goes_nc.time.values[0])[0:10] + '_' + str(goes_nc.time.values[-1])[0:10] + '.nc', format='NETCDF3_CLASSIC')
    goes_nc.to_netcdf(path=outpath + '/' + str(datelist[d].year) + '/GOES16_SST_dailycomposite_' + str(datelist[d].year) + str("{0:0=3d}".format(datelist[d].dayofyear)) + '.nc', format='NETCDF3_CLASSIC')


