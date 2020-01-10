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


# grab the buoy data and throw it into temporary data frame
for s in list(stations.keys()):
    print(s)
    year=2019
    statsDF[s]['Name'] = stations[s]
    # grab the buoy buoy data 
    
    buoy_nc = xr.open_dataset('https://dods.ndbc.noaa.gov/thredds/dodsC/data/stdmet/' + s + '/' + s + 'h' + str(year) + '.nc')
    buoy_nc = buoy_nc.sel(time = slice('2019-01-01', '2019-10-30'))
    #wind_nc = buoy_nc['wind_spd']
    buoy_nc = buoy_nc['sea_surface_temperature']
    buoy_nc = buoy_nc.where(buoy_nc > 2, drop=True) # 2 degrees celsius
    buoy_nc.values = buoy_nc.values + 273.15
    # grab goes 16 data at the same location and throw it into dataframe
    
    lat_s = buoy_nc.latitude.values[0]
    lon_s = buoy_nc.longitude.values[0]

    goes_nc = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/goes_r_sst.nc')
    goes_nc =goes_nc.drop(time=[np.datetime64('2000-01-01T11:43:21.000000000')])
    goes_nc =goes_nc.drop(time=[np.datetime64('2009-08-19T19:56:53.000000000')])
    goes_nc =goes_nc.drop(time=[np.datetime64('2009-08-20T07:56:53.000000000')])
    goes_nc =goes_nc.drop(time=[np.datetime64('2009-08-20T19:56:51.000000000')])
    goes_nc =goes_nc.drop(time=[np.datetime64('2009-09-27T10:52:10.000000000')])
    goes_nc =goes_nc.drop(time=[np.datetime64('2009-11-06T14:51:48.000000000')])
    goes_nc =goes_nc.drop(time=[np.datetime64('2009-11-30T16:51:59.000000000')])
    goes_nc = goes_nc.drop('Band15')
    goes_nc = goes_nc.sel(time = slice('2019-01-01', '2019-10-30'))
    goes_nc = goes_nc.sel(latitude=lat_s, longitude=lon_s, method='nearest')

    

    
    # match up the times for when there is good data for both
    if len(buoy_nc.time.values) > 1:
        # commented out because it's slower actually? weird
        #x=datetime.datetime.utcnow()
        #goes_nc.SST = goes_nc.where(goes_nc['DQF'] == 0)
        #datetime.datetime.utcnow() - x
        
        for t in range(len(goes_nc.time.values)):
            x = goes_nc['SST'][t]
            goes_nc['SST'][t] = x.where(goes_nc['DQF'][t] == 0)

        goes_nc = goes_nc.drop(['DQF'])
        goes_nc = goes_nc.metpy.parse_cf("SST")
        goes_nc = goes_nc.where(goes_nc.values > 273.15, np.nan)
        sst_time = []
        buoy_time = []
        sst_vals = []
        buoy_vals = []
        wind_vals = []
        
        for val in range(len(goes_nc.values)):
            sstTime = goes_nc.time.values[[val]]
            buoyTimeInd = find_nearest(buoy_nc.time.values,sstTime[0])
            #windTimeInd = find_nearest(wind_nc.time.values[buoyTimeInd], buoy_nc.time)
            timediff = buoy_nc.time.values[buoyTimeInd] - sstTime[0]
            timediff = timediff.astype('timedelta64[m]')
            timediff = np.abs(timediff / np.timedelta64(1, 'm'))
            if timediff < 60:
                if math.isnan(goes_nc.values[val]) == False:
                    if math.isnan(buoy_nc.values[buoyTimeInd]) == False:
                        sst_time.append(goes_nc.time.values[val])
                        buoy_time.append(buoy_nc.time.values[buoyTimeInd])
                        sst_vals.append(goes_nc.values[val])
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
        
            fig = plt.figure(figsize=(10,6))
            plt.scatter(sst_time, sst_vals, color='red', s=5, label = 'GOES SST')
            plt.scatter(sst_time, buoy_vals, color='blue', s=5, label = 'Buoy Value')
            plt.legend()
            plt.title('Buoy ' + s + ' ' + stations[s])
            plt.savefig("/Users/james/Documents/buoy_val/goes_images/buoy" + s)
            plt.close()
        
        #fig = plt.figure(figsize=(16,8))
        #plt.scatter(sst_time, (np.array(sst_vals) - np.array(buoy_vals)), color='green', s=5, label = 'SST Difference')
        #plt.scatter(sst_time, wind_vals, color='black', s=5, label = 'Wind Speed (Buoy)')
        #plt.legend()
        #plt.title('Buoy ' + s + ' ' + stations[s])
        #plt.savefig("/Users/james/Documents/buoy_val/wind_goes/buoy" + s)
        #plt.close()
statsDF.to_csv("/Users/james/Documents/buoy_val/goes_buoyVal.csv")









'''

from erddapy import ERDDAP

s = list(stations.keys())[-1]
e = ERDDAP(
  server='https://coastwatch.pfeg.noaa.gov/erddap',
  protocol='tabledap',
)

e.response = 'csv'
e.dataset_id = 'cwwcNDBCMet'
e.constraints = {
    'time>=': '2017-01-01T00:00:00Z',
    'time<=': '2017-01-14T00:00:00Z',
    'station=': s,
}
e.variables = [
    'time',
    'latitude',
    'longitude',
    'wtmp',
]

df = e.to_pandas()
'''