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

statsDF = pd.DataFrame(columns=list(stations.keys()), index=['Name', 'RMSE', 'MeanSquareError', 'MeanAbsoluteError', 'Rsquared', 'median_absolute_error', 'explained_variance_score', 'max_error','Bias(Sat-Buoy)', 'count'])

years_viirs= np.arange(2002,2020,1)
# grab the buoy data and throw it into temporary data frame
for s in list(stations.keys()):
    print(s)
    statsDF[s]['Name'] = stations[s]
    # grab the buoy buoy data 
    
    success = []
    datasets = []
    for year in years_viirs:
        try:
            ds = xr.open_dataset('https://dods.ndbc.noaa.gov/thredds/dodsC/data/stdmet/' + s + '/' + s + 'h' + str(year) + '.nc')
        except:
            pass
        else:
            # Add the following in case buoy location changes slightly by year
            success.append(year)
            if year == success[0]:
                LAT = ds.latitude.values[0]
                LON = ds.longitude.values[0]
            ds=ds.assign_coords(latitude=[LAT])
            ds=ds.assign_coords(longitude=[LON])
            datasets.append(ds)
    
    combined = xr.concat(datasets, dim='time')
    buoy_nc = combined.sel(time = slice('2012-01-01', '2019-10-30'))
    #wind_nc = buoy_nc['wind_spd']
    #wind_nc = wind_nc.resample(time='1D').mean('time')
    buoy_nc = buoy_nc['sea_surface_temperature']
    buoy_nc = buoy_nc.resample(time='1D').mean('time')
    buoy_nc = buoy_nc.where(buoy_nc > 2, drop=True) # 2 degrees celsius
    lat_s = buoy_nc.latitude.values[0]
    lon_s = buoy_nc.longitude.values[0]
    # grab modis viirs data at the same location and throw it into dataframe
    
    viirs_nc = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/viirs_1day_aggregate.nc')
    viirs_nc = viirs_nc = viirs_nc.metpy.parse_cf("sst")
    viirs_nc = viirs_nc.sel(time = slice('2012-01-01', '2019-10-30'))
    viirs_nc = viirs_nc.sel(lat=lat_s, lon=lon_s, method='nearest')
    viirs_nc = viirs_nc.where(viirs_nc > 2, drop=True) # 2 degrees celsius
    viirs_nc = viirs_nc.where(viirs_nc < 40, drop=True) # 2 degrees celsius
    
    # match up the times for when there is good data for both
    if len(buoy_nc.time.values) > 1:

        
        sst_time = []
        buoy_time = []
        sst_vals = []
        buoy_vals = []
        wind_vals = []
        
        for val in range(len(viirs_nc.values)):
            sstTime = viirs_nc.time.values[[val]]
            buoyTimeInd = find_nearest(buoy_nc.time.values,sstTime[0])
            #windTimeInd = find_nearest(wind_nc.time.values[buoyTimeInd], buoy_nc.time)
            
            timediff = buoy_nc.time.values[buoyTimeInd] - sstTime[0]
            timediff = timediff.astype('timedelta64[m]')
            timediff = np.abs(timediff / np.timedelta64(1, 'm'))
            if timediff < 1440:
                if math.isnan(viirs_nc.values[val]) == False:
                    if math.isnan(buoy_nc.values[buoyTimeInd]) == False:
                        sst_time.append(viirs_nc.time.values[val])
                        buoy_time.append(buoy_nc.time.values[buoyTimeInd])
                        sst_vals.append(viirs_nc.values[val])
                        buoy_vals.append(buoy_nc.values[buoyTimeInd][0][0])
                        #wind_vals.append(wind_nc.values[windTimeInd][0][0])
            
        if len(buoy_vals) > 1:
            statsDF[s]['RMSE'] = sqrt(metrics.mean_squared_error(buoy_vals, sst_vals))
            statsDF[s]['MeanSquareError'] = metrics.mean_squared_error(buoy_vals, sst_vals)
            statsDF[s]['MeanAbsoluteError'] = metrics.mean_absolute_error(buoy_vals, sst_vals)
            statsDF[s]['Rsquared'] = metrics.r2_score(buoy_vals, sst_vals)
            statsDF[s]['explained_variance_score'] = metrics.explained_variance_score(buoy_vals, sst_vals)
            statsDF[s]['median_absolute_error'] = metrics.median_absolute_error(buoy_vals, sst_vals)
            statsDF[s]['max_error'] = metrics.max_error(buoy_vals, sst_vals)
            statsDF[s]['Bias(Sat-Buoy)'] = ((np.sum(sst_vals) - np.sum(buoy_vals)) * (1.0/len(sst_vals)))
            statsDF[s]['count'] = len(sst_vals)
            
            fig = plt.figure(figsize=(16,8))
            plt.scatter(sst_time, sst_vals, color='red', s=5, label = 'NPP VIIRS SST')
            plt.scatter(sst_time, buoy_vals, color='blue', s=5, label = 'Buoy Value')
            plt.legend()
            plt.title('Buoy ' + s + ' ' + stations[s])
            plt.savefig("/Users/james/Documents/buoy_val/viirs_images/buoy" + s)
            plt.close()
        #fig = plt.figure(figsize=(16,8))
        #plt.scatter(sst_time, (np.array(sst_vals) - np.array(buoy_vals)), color='green', s=5, label = 'SST Difference')
        #plt.scatter(sst_time, wind_vals, color='black', s=5, label = 'Wind Speed (Buoy)')
        #plt.legend()
        #plt.title('Buoy ' + s + ' ' + stations[s])
        #plt.savefig("/Users/james/Documents/buoy_val/wind_viirs/buoy" + s)
        #plt.close()
statsDF.to_csv("/Users/james/Documents/buoy_val/viirs_buoyVal.csv")




