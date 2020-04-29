# backfill all of the 1day goes sst files that we missed
import xarray as xr
import numpy as np
from datetime import datetime, timedelta
import pandas as pd
# paths
outpath1 = "/data/GOES/GOES-R/1day/"
outpath2 = "/data/GOES/GOES-R/daily_composite/"
datelist = pd.date_range('2018-01-01', pd.datetime.today()).tolist()

# for each day, process a file
for d in range(0, len(datelist)):
    print(datelist[d])
    #open the main THREDDS instance
    goes_main = xr.open_dataset(
        "http://basin.ceoe.udel.edu/thredds/dodsC/GOESJPLSST.nc")
    #subset to day we want to work with
    goes_nc = goes_main.sel(time=datetime.strftime(
        datelist[d].date(), '%Y-%m-%d'))
    #drop unnecessary variables
    goes_nc = goes_nc.drop('projection')
    goes_nc = goes_nc.drop('dt_analysis')
    goes_nc = goes_nc.drop('satellite_zenith_angle')
    goes_nc = goes_nc.drop('sses_standard_deviation')
    goes_nc = goes_nc.drop('wind_speed')
    #subtract the bias before we delete the bias variable
    goes_nc['sea_surface_temperature'] = goes_nc['sea_surface_temperature'] - goes_nc['sses_bias']
    goes_nc = goes_nc.drop('sses_bias')
    #mark the timestamp as the most recent timestamp
    newtimestamp = goes_nc.time.values[-1]
    #only accept SST value where we have clear sky
    for t in range(len(goes_nc.time.values)):
        x = goes_nc['sea_surface_temperature'][t]
        goes_nc['sea_surface_temperature'][t] = x.where(goes_nc['quality_level'][t] == 5)

    #rename sea_surface_temperature sst so things don't break downstream
    goes_nc['sst'] = goes_nc['sea_surface_temperature']
    goes_nc = goes_nc.drop(['sea_surface_temperature'])
    # have to add the following line becuase of a weird xarray netcdf4 error when writing the xarray to netcdf
    del goes_nc.attrs['_NCProperties']
    goes_nc.to_netcdf(path=outpath1 + '/' + str(datelist[d].year) + '/GOES16_SST_1day_' + str(datelist[d].year) + str("{0:0=3d}".format(
        datelist[d].dayofyear)) + '_' + str("{0:0=2d}".format(datelist[d].month)) + str("{0:0=2d}".format(datelist[d].day)) + '.nc', mode='w',format='NETCDF4')

    # resample to a daily composite
    goes_nc = goes_nc.drop(['quality_level'])
    goes_nc = goes_nc.resample(time='1D').mean('time')
    # convert to celsius
    goes_nc['sst'] = goes_nc['sst'] - 273.15
    # have to do weird things with the timestamp since we resampled to daily resolution
    newtimestamp = (newtimestamp - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')
    goes_nc['time'] = np.array([newtimestamp], dtype='float64')
    # set the units in the attributes
    goes_nc.time.attrs['units'] = 'seconds since 1970-01-01 00:00:00'
    goes_nc.sst.attrs['units'] = "Celsius"
    goes_nc.sst.attrs['data_summary'] = "Daily composite of SST with clear sky pixels (quality level = 5, or highest quality level data)"
    #landmask.to_netcdf(path=outpath2 + 'landmask_roffs_' +  'area' + str(a) + '_' + str(goes_nc.time.values[0])[0:10] + '_' + str(goes_nc.time.values[-1])[0:10] + '.nc', format='NETCDF3_CLASSIC')
    goes_nc.to_netcdf(path=outpath2 + '/' + str(datelist[d].year) + '/GOES16_SST_dailycomposite_' + str(datelist[d].year) + str("{0:0=3d}".format(
        datelist[d].dayofyear)) + '_' + str("{0:0=2d}".format(datelist[d].month)) + str("{0:0=2d}".format(datelist[d].day)) + '.nc', mode='w',format='NETCDF4')
