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
#stations = {'44027' : 'Jonesport, ME',
#  '44007' : 'Portland, ME',
#  '44005' : 'Gulf of Maine',
#  '44013' : 'Boston, MA',
#  '44018' : 'Cape Cod, MA',
#  '44020' : 'Nantucket Sound',
#  '44011' : 'Georges Bank', 
#  '44008' : 'Nantucket', 
#  '44017' : 'Montauk Point',
#  '44065' : 'New York Harbor Entrance',
#  '44025' : 'Long Island',
#  '44066' : 'Texas Tower',
#  '44009' : 'Delaware Bay',
#  '44089' : 'Wallops Island'
#}

stations = {'44009' : 'Delaware Bay'}

statsDF = pd.DataFrame(columns=list(stations.keys()), index=['Name', 'RMSE', 'MeanSquareError', 'MeanAbsoluteError', 'Rsquared', 'median_absolute_error', 'explained_variance_score', 'max_error','Bias(Sat-Buoy)', 'count'])

years_aqua = np.arange(2002,2020,1)
# grab the buoy data and throw it into temporary data frame
for s in list(stations.keys()):
    print(s)
    statsDF[s]['Name'] = stations[s]
    # grab the buoy buoy data 
    
    success = []
    datasets = []
    for year in years_aqua:
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
    buoy_nc = combined.sel(time = slice('2002-01-01', '2019-10-30'))
    #wind_nc = buoy_nc['wind_spd']
    #wind_nc = wind_nc.resample(time='1D').mean('time')
    buoy_nc = buoy_nc['sea_surface_temperature']
    buoy_nc = buoy_nc.resample(time='1D').mean('time')
    buoy_nc = buoy_nc.where(buoy_nc > 0, drop=True) # 2 degrees celsius
    lat_s = buoy_nc.latitude.values[0]
    lon_s = buoy_nc.longitude.values[0]
    # grab modis aqua data at the same location and throw it into dataframe
    
    aqua_nc = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/Aqua1DayAggregate.nc')
    aqua_nc = aqua_nc.metpy.parse_cf("sst")
    aqua_nc = aqua_nc.sel(time = slice('2002-01-01', '2019-10-30'))
    aqua_nc = aqua_nc.sel(lat=lat_s, lon=lon_s, method='nearest')
    aqua_nc = aqua_nc.where(aqua_nc > 0, drop=True) # 2 degrees celsius
    
    
    # match up the times for when there is good data for both
    if len(buoy_nc.time.values) > 1:

        
        sst_time = []
        buoy_time = []
        sst_vals = []
        buoy_vals = []
        wind_vals = []
        
        for val in range(len(aqua_nc.values)):
            sstTime = aqua_nc.time.values[[val]]
            buoyTimeInd = find_nearest(buoy_nc.time.values,sstTime[0])
            #windTimeInd = find_nearest(wind_nc.time.values[buoyTimeInd], buoy_nc.time)
            timediff = buoy_nc.time.values[buoyTimeInd] - sstTime[0]
            timediff = timediff.astype('timedelta64[m]')
            timediff = np.abs(timediff / np.timedelta64(1, 'm'))
            if timediff < 1440:
                if math.isnan(aqua_nc.values[val]) == False:
                    if math.isnan(buoy_nc.values[buoyTimeInd]) == False:
                        sst_time.append(aqua_nc.time.values[val])
                        buoy_time.append(buoy_nc.time.values[buoyTimeInd])
                        sst_vals.append(aqua_nc.values[val])
                        buoy_vals.append(buoy_nc.values[buoyTimeInd][0][0])
                        #wind_vals.append(wind_nc.values[windTimeInd][0][0])

            #fig = plt.figure(figsize=(16,8))
            #plt.scatter(sst_time, sst_vals, color='red', s=5, label = 'aqua SST')
            #plt.scatter(sst_time, buoy_vals, color='blue', s=5, label = 'Buoy Value')
            #plt.legend()
            #plt.title('Buoy ' + s + ' ' + stations[s])
            #plt.savefig("/Users/james/Documents/buoy_val/aqua_images/buoy" + s)
            #plt.close()
        
        # add trendline here and plot the bias over the length of the aqua dataset here
        fig = plt.figure(figsize=(16,8))
        plt.scatter(sst_time, (np.array(sst_vals) - np.array(buoy_vals)), color='green', s=5, label = 'SST Bias: (Sat-Buoy)')
        #plt.scatter(sst_time, wind_vals, color='black', s=5, label = 'Wind Speed (Buoy)')
        plt.legend()
        plt.title('Buoy ' + s + ' ' + stations[s])
        #plt.savefig("/Users/james/Documents/buoy_val/wind_aqua/buoy" + s)
        plt.close()

#statsDF.to_csv("/Users/james/Documents/buoy_val/aqua_buoyVal.csv")







