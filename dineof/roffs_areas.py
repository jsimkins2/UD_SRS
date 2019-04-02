# prepare goes16 sst data for dineof
import xarray as xr
import numpy as np
import metpy
from datetime import datetime, timedelta

# paths
outpath = "/home/sat_ops/goesR/data/sst/roffs/"

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
daysback = [7,14,21,28,35,42] #must keep in brackets to python recognizes it as a list

for a in range(1,2):
    for d in range(0,len(daysback)):
        #print(a)
        area = areas['area' + str(a)]
        if d == 0:
            dayOffset = 0
            addOffset = 0
        
        # grab sst data from the last 3 days and use DQF == 0 data
        goes_nc = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/goes_r_sst_daily.nc")
        goes_nc = goes_nc.sel(latitude=slice(area[0],area[1]), longitude=slice(area[2], area[3]), time=slice(datetime.strftime(nowday - timedelta(days=daysback[d] + addOffset), '%Y-%m-%d'), datetime.strftime(nowday - timedelta(days=dayOffset), '%Y-%m-%d')))

        # save multiple 1 week intervals
        #dayOffset = daysback[d] + 1
        #addOffset = 1

        
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
        #forecast_nc = goes_nc.isel(time=[-1])
        #forecast_nc.time.values = forecast_nc.time.values.astype('datetime64[s]') + (3600)
        #forecast_nc['sst'] = forecast_nc['sst'].where(forecast_nc['sst'] < 2)
        #goes_nc= xr.concat([goes_nc, forecast_nc], dim='time')

        landmask.to_netcdf(path=outpath + 'landmask_roffs_area' + str(a) '.nc', format='NETCDF3_CLASSIC')
        goes_nc.to_netcdf(path=outpath + 'roffs_' +  'area' + str(a) + '_' + str(goes_nc.time.values[0])[0:10] + '_' + str(goes_nc.time.values[-1])[0:10] + '.nc', format='NETCDF3_CLASSIC')
    
    
    
    
    
    
    
    

