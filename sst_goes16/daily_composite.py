# create 1-day and daily composites as they come in
import xarray as xr
import numpy as np
from datetime import datetime, timedelta
import pandas as pd
# paths
outpath1 = "/data/GOES/GOES-R/1day/"
outpath2 = "/data/GOES/GOES-R/daily_composite/"
datelist = pd.date_range(pd.datetime.today() - timedelta(days=7), pd.datetime.today() - timedelta(days=1)).tolist()
for d in range(0, len(datelist)):
    print(datelist[d])
    goes_main = xr.open_dataset(
        "http://basin.ceoe.udel.edu/thredds/dodsC/GOESJPLSST.nc")
    goes_nc = goes_main.sel(time=datetime.strftime(
        datelist[d].date(), '%Y-%m-%d'))
    goes_nc = goes_nc.drop('projection')
    goes_nc = goes_nc.drop('dt_analysis')
    goes_nc = goes_nc.drop('satellite_zenith_angle')
    goes_nc = goes_nc.drop('sses_standard_deviation')
    goes_nc = goes_nc.drop('wind_speed')
    goes_nc['sea_surface_temperature'] = goes_nc['sea_surface_temperature'] - goes_nc['sses_bias']
    goes_nc = goes_nc.drop('sses_bias')
    newtimestamp = goes_nc.time.values[-1]
    for t in range(len(goes_nc.time.values)):
        x = goes_nc['sea_surface_temperature'][t]
        goes_nc['sea_surface_temperature'][t] = x.where(goes_nc['quality_level'][t] == 5)
        #x = goes_nc['DQF'][t]
        #goes_nc['DQF'][t] = x.where(goes_nc['DQF'][t] == 3)
    
    goes_nc['sst'] = goes_nc['sea_surface_temperature']
    goes_nc = goes_nc.drop(['sea_surface_temperature'])
    # have to add the following line becuase of a weird xarray netcdf4 error when writing the xarray to netcdf
    del goes_nc.attrs['_NCProperties']
    goes_nc.to_netcdf(path=outpath1 + '/' + str(datelist[d].year) + '/GOES16_SST_1day_' + str(datelist[d].year) + str("{0:0=3d}".format(
        datelist[d].dayofyear)) + '_' + str("{0:0=2d}".format(datelist[d].month)) + str("{0:0=2d}".format(datelist[d].day)) + '.nc', mode='w',format='NETCDF4')

    # resample to a daily composite
    goes_nc = goes_nc.drop(['quality_level'])
    goes_nc = goes_nc.resample(time='1D').mean('time')
    
    newtimestamp = (newtimestamp - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')
    goes_nc['time'] = np.array([newtimestamp], dtype='float64')
    goes_nc.time.attrs['units'] = 'seconds since 1970-01-01 00:00:00'
    
    #landmask.to_netcdf(path=outpath2 + 'landmask_roffs_' +  'area' + str(a) + '_' + str(goes_nc.time.values[0])[0:10] + '_' + str(goes_nc.time.values[-1])[0:10] + '.nc', format='NETCDF3_CLASSIC')
    goes_nc.to_netcdf(path=outpath2 + '/' + str(datelist[d].year) + '/GOES16_SST_dailycomposite_' + str(datelist[d].year) + str("{0:0=3d}".format(
        datelist[d].dayofyear)) + '_' + str("{0:0=2d}".format(datelist[d].month)) + str("{0:0=2d}".format(datelist[d].day)) + '.nc', mode='w',format='NETCDF4')



