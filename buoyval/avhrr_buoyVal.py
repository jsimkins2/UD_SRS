import xarray as xr
import numpy as np
import math
import pandas as pd
import metpy
import sklearn.metrics as metrics #standard deviation of the residuals
from sklearn.ensemble import RandomForestRegressor as rfr
from math import sqrt
import matplotlib.pyplot as plt
from datetime import datetime, timedelta


## set up initial data frame


def find_nearest(array, value):
    ''' Find nearest value is an array '''
    idx = (np.abs(array-value)).argmin()
    return idx

# from Maine down to florida
stations = {'44027' : 'Jonesport, ME',
  '44007' : 'Portland, ME',
  '44005' : 'Gulf of Maine',
  '44013' : 'Boston, MA',
  '44018' : 'Cape Cod, MA',
  '44020' : 'Nantucket Sound',
  '44011' : 'Georges Bank', 
  '44008' : 'Nantucket', 
  '44017' : 'Montauk Point',
  '44065' : 'New York Harbor Entrance',
  '44025' : 'Long Island',
  '44066' : 'Texas Tower',
  '44009' : 'Delaware Bay',
  '44089' : 'Wallops Island'
  
}

statsDF = pd.DataFrame(columns=list(stations.keys()), index=['Name', 'RMSE', 'MSE', 'MAE', 'Rsquared', 'median_absolute_error', 'explained_variance_score', 'max_error'])


# grab the buoy data and throw it into temporary data frame
for s in list(stations.keys()):
    print(s)
    year=2019
    statsDF[s]['Name'] = stations[s]
    # grab the buoy buoy data 
    
    buoy_nc = xr.open_dataset('https://dods.ndbc.noaa.gov/thredds/dodsC/data/stdmet/' + s + '/' + s + 'h' + str(year) + '.nc')
    buoy_nc = buoy_nc.sel(time = slice('2019-01-01', '2019-10-15'))
    buoy_nc = buoy_nc['sea_surface_temperature']
    # grab avhrr 16 data at the same location and throw it into dataframe
    lat_s = buoy_nc.latitude.values[0]
    lon_s = buoy_nc.longitude.values[0]
    
    avhrr_nc = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/avhrr_unfiltered_sst.nc')
    badtimes = avhrr_nc.time.values[6518:6545]
    for b in badtimes:
        avhrr_nc =avhrr_nc.drop(time=[np.datetime64(b)])
    
    avhrr_nc = avhrr_nc.sel(time = slice('2019-01-01', '2019-10-15'))
    avhrr_nc = avhrr_nc.sel(lat=lat_s, lon=lon_s, method='nearest')
    
    
    
    # match up the times for when there is good data for both
    if len(buoy_nc.time.values) > 1:
        for t in range(len(avhrr_nc.time.values)):
            x = avhrr_nc['mcsst'][t]
            avhrr_nc['mcsst'][t] = x.where(avhrr_nc['cloud_land_mask'][t] == 0)
        
        avhrr_nc = avhrr_nc.drop(['cloud_land_mask'])
        avhrr_nc = avhrr_nc.metpy.parse_cf("mcsst")
        avhrr_nc.values = avhrr_nc.values
        avhrr_nc = avhrr_nc.where(avhrr_nc.values > 0, np.nan)
        sst_time = []
        buoy_time = []
        sst_vals = []
        buoy_vals = []
        for val in range(len(avhrr_nc.values)):
            sstTime = avhrr_nc.time.values[[val]]
            buoyTimeInd = find_nearest(buoy_nc.time.values,sstTime[0])
            
            timediff = buoy_nc.time.values[buoyTimeInd] - sstTime[0]
            timediff = timediff.astype('timedelta64[m]')
            timediff = timediff / np.timedelta64(1, 'm')
            if timediff < 60*3:
                if math.isnan(avhrr_nc.values[val]) == False:
                    if math.isnan(buoy_nc.values[buoyTimeInd]) == False:
                        sst_time.append(avhrr_nc.time.values[val])
                        buoy_time.append(buoy_nc.time.values[buoyTimeInd])
                        sst_vals.append(avhrr_nc.values[val])
                        buoy_vals.append(buoy_nc.values[buoyTimeInd][0][0])
            
        if len(buoy_vals) > 1:
            statsDF[s]['RMSE'] = sqrt(metrics.mean_squared_error(buoy_vals, sst_vals))
            statsDF[s]['MSE'] = metrics.mean_squared_error(buoy_vals, sst_vals)
            statsDF[s]['MAE'] = metrics.mean_absolute_error(buoy_vals, sst_vals)
            statsDF[s]['Rsquared'] = metrics.r2_score(buoy_vals, sst_vals)
            statsDF[s]['explained_variance_score'] = metrics.explained_variance_score(buoy_vals, sst_vals)
            statsDF[s]['median_absolute_error'] = metrics.median_absolute_error(buoy_vals, sst_vals)
            statsDF[s]['max_error'] = metrics.max_error(buoy_vals, sst_vals)
            statsDF[s]['count'] = len(sst_vals)
            
            fig = plt.figure(figsize=(10,6))
            plt.scatter(sst_time, sst_vals, color='red', s=5, label = 'avhrr SST')
            plt.scatter(sst_time, buoy_vals, color='blue', s=5, label = 'Buoy Value')
            plt.legend()
            plt.title('Buoy ' + s + ' ' + stations[s])
            plt.savefig("/Users/james/Documents/buoy_val/avhrr_images/buoy" + s)


statsDF.to_csv("/Users/james/Documents/buoy_val/avhrr2019.csv")



