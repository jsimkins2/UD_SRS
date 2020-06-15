# prepare goes16 sst data for dineof
import xarray as xr
import numpy as np
#import metpy
from datetime import datetime, timedelta
import pandas as pd


# define area of interest
area1 = [32.5, 36, -78.25, -73.75]
area2 = [36.5, 40, -75.25, -72.0]
area3 = [27.0, 30.5, -91.5, -85.5]
area4 = [26.5, 29.5, -91.0, -86.5]
area5 = [34, 36, -76, -74]

# create a dictionary so these areas can be called via for loop
areas = {'area1': area1,
         'area2': area2,
         'area3': area3,
         'area4': area4,
         'area5': area5}

nowday = datetime.utcnow()
daysback = [7]#,14,21] #must keep in brackets to python recognizes it as a list

for a in range(2,3):
    for d in range(0,len(daysback)):
        #print(a)
        area = areas['area' + str(a)]
        if d == 0:
            dayOffset = 0
            addOffset = 0

        # grab sst data from the last d days and use DQF == 0 data
        print(d)
        goes_main = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/GOESJPLSST.nc")
        goes_nc = goes_main.sel(latitude=slice(area[1],area[0]), longitude=slice(area[2], area[3]), time=slice(datetime.strftime(nowday - timedelta(days=daysback[d]), '%Y-%m-%d %H:%M'), datetime.strftime(nowday, '%Y-%m-%d %H:%M')))
        # save multiple 1 week intervals
        #dayOffset = daysback[d] + 1
        #addOffset = 1

        # COMMENTING OUT BELOW BECAUSE THIS STEP IS ALREADY BEING DONE IN 1DAY
        for t in range(len(goes_nc.time.values)):
            x = goes_nc['sea_surface_temperature'][t]
            goes_nc['sea_surface_temperature'][t] = x.where(goes_nc['quality_level'][t] == 5)
            x = goes_nc['quality_level'][t]
            goes_nc['quality_level'][t] = x.where(goes_nc['quality_level'][t] == 0)

        #landmask = goes_nc['DQF'][0]
        #landmask = landmask.rename('landmask')
        #landmask = landmask.where(landmask.values == 3, 1)
        #landmask = landmask.where(landmask.values == 1, 0)
        goes_nc = goes_nc.drop(['quality_level'])
        goes_SST = goes_nc["sea_surface_temperature"]
        goes_SST.values = goes_SST.values - goes_nc['sses_bias'].values

        # Clean out files that are missing too much data
        '''
        print(len(goes_nc.time.values))
        badfiles=[]
        for tval in range(len(goes_nc.time.values)):
            sst = goes_nc.metpy.parse_cf("sst")[tval]
            sst = sst.where(sst.values > 270, np.nan)
            sst = sst.fillna(-999)
            if np.percentile(sst.values, 95) == -999:
                badfiles.append(sst.time.values)

        for t in badfiles:
            temArray = goes_nc.time.values
            temElem = np.where(temArray == t)
            print(temArray[temElem])
            print(goes_nc.time.values[temElem[0]])
            goes_nc = goes_nc.drop(goes_nc.time.values[temElem[0]], dim='time')

        print(len(goes_nc.time.values))
        '''

        # Add a forecast day
        forecast_nc = goes_SST.isel(time=[-1])
        forecast_nc = forecast_nc.assign_coords(time=forecast_nc['time'].values.astype('datetime64[s]') + (3600))
        forecast_nc = forecast_nc.where(forecast_nc < 0)
        goes_SST= xr.concat([goes_SST, forecast_nc], dim='time')
        
        # paths
        outpath = "/Users/james/Downloads/"
        #landmask.to_netcdf(path=outpath + 'landmask_roffs_area' + str(a) + '.nc', format='NETCDF3_CLASSIC')
        goes_SST.to_netcdf(path=outpath + 'roffs_' +  'area' + str(a) + '_' + str(daysback[d]) + "day"  + '.nc', format='NETCDF3_CLASSIC')
        '''
        if daysback[d] == 7:
            outpath = '/home/sat_ops/goesR/data/sst/roffs/area' + str(a) + '/'
            goes_nc.to_netcdf(path=outpath + 'roffs_' +  'area' + str(a) + '_'
            + datetime.strftime(pd.to_datetime(goes_nc.time.values[0]), '%Y%j%H%M') + '_' 
            + datetime.strftime(pd.to_datetime(goes_nc.time.values[-1]), '%Y%j%H%M') + '.nc', format='NETCDF3_CLASSIC')


	'''
