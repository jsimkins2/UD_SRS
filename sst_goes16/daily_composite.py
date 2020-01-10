# backfill all of the 1day goes sst files that we missed
import xarray as xr
import numpy as np
import metpy
from datetime import datetime, timedelta
import pandas as pd
# paths
outpath1 = "/data/GOES/GOES-R/1day/"
outpath2 = "/data/GOES/GOES-R/daily_composite/"
datelist = pd.date_range('2020-01-01', pd.datetime.today()).tolist()
for d in range(0, len(datelist)):
    print(datelist[d])
    goes_nc = xr.open_dataset(
        "http://basin.ceoe.udel.edu/thredds/dodsC/goes_r_sst.nc")
    goes_nc = goes_nc.sel(time=datetime.strftime(
        datelist[d].date(), '%Y-%m-%d'))
    goes_nc = goes_nc.drop('Band15')
    newtimestamp = goes_nc.time.values[-1]

    for t in range(len(goes_nc.time.values)):
        x = goes_nc['SST'][t]
        goes_nc['SST'][t] = x.where(goes_nc['DQF'][t] == 0)
        #x = goes_nc['DQF'][t]
        #goes_nc['DQF'][t] = x.where(goes_nc['DQF'][t] == 3)
    
    
    goes_nc['sst'] = goes_nc['SST']
    goes_nc = goes_nc.drop(['SST'])
    goes_nc.to_netcdf(path=outpath1 + '/' + str(datelist[d].year) + '/GOES16_SST_1day_' + str(datelist[d].year) + str("{0:0=3d}".format(
        datelist[d].dayofyear)) + '_' + str("{0:0=2d}".format(datelist[d].month)) + str("{0:0=2d}".format(datelist[d].day)) + '.nc', format='NETCDF3_CLASSIC')

    # resample to a daily composite
    goes_nc = goes_nc.drop(['DQF'])
    goes_nc = goes_nc.resample(time='1D').mean('time')
    
    newtimestamp = (newtimestamp - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
    goes_nc.time.values = np.array([newtimestamp], dtype='float64')
    goes_nc.time.attrs['units'] = 'seconds since 1970-01-01 00:00:00'
    
    #landmask.to_netcdf(path=outpath2 + 'landmask_roffs_' +  'area' + str(a) + '_' + str(goes_nc.time.values[0])[0:10] + '_' + str(goes_nc.time.values[-1])[0:10] + '.nc', format='NETCDF3_CLASSIC')
    goes_nc.to_netcdf(path=outpath2 + '/' + str(datelist[d].year) + '/GOES16_SST_dailycomposite_' + str(datelist[d].year) + str("{0:0=3d}".format(
        datelist[d].dayofyear)) + '_' + str("{0:0=2d}".format(datelist[d].month)) + str("{0:0=2d}".format(datelist[d].day)) + '.nc', format='NETCDF3_CLASSIC')

    
