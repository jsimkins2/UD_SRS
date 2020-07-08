from datetime import datetime, timedelta
import geopandas as gpd
import rioxarray
import xarray as xr
import pandas as pd
import time
from calendar import monthrange
import os
import numpy as np
#################################################################
# declare paths
#################################################################
temPath = "/home/sat_ops/deos/temp/"
outPathNC = "/data/DEOS/agriculture/"
#################################################################
# declare in functions
#################################################################
def linear_rbf(x, y, z, xi, yi):
    dist = distance_matrix(x,y, xi,yi)

    # Mutual pariwise distances between observations
    internal_dist = distance_matrix(x,y, x,y)

    # Now solve for the weights such that mistfit at the observations is minimized
    weights = np.linalg.solve(internal_dist, z)

    # Multiply the weights for each interpolated point by the distances
    zi =  np.dot(dist.T, weights)
    return zi

def distance_matrix(x0, y0, x1, y1):
    obs = np.vstack((x0, y0)).T
    interp = np.vstack((x1, y1)).T

    # Make a distance matrix between pairwise observations
    # Note: from <http://stackoverflow.com/questions/1871536>
    # (Yay for ufuncs!)
    d0 = np.subtract.outer(obs[:,0], interp[:,0])
    d1 = np.subtract.outer(obs[:,1], interp[:,1])

    return np.hypot(d0, d1)
    
#### custom colormap  
import matplotlib as mpl
startcolor = '#8B4513'
midcolor = '#FFFFFF'
endcolor = '#008000'
own_cmap1 = mpl.colors.LinearSegmentedColormap.from_list( 'own2', [startcolor, midcolor, endcolor] )

#################################################################
# begin script
#################################################################

# grab data from json files
deos_data = pd.read_json("http://128.175.28.202/deos_json/map_data2.json")
loc_deos = pd.read_json("http://128.175.28.202/deos_json/station_metadata.json")
# deos_data['2020-02-06 18:30:00'][2302]
# index is the station nunumbers 
station_id = list(deos_data.index)
# need to find time value
date_deos = deos_data.columns[0]
dst = "EST" if time.localtime().tm_isdst==0 else "EDT"
zuluDIFF = 5 if dst=='EST' else 4
deos_dateSTR = str(str(date_deos.month) + '/' + str(date_deos.day) + '/' + str(date_deos.year) + ' ' + str(date_deos.hour - zuluDIFF) + ':' + str(date_deos.minute) + ' ' + dst)
# create a dict of station IDs with the name of the station
station_dict = {}
for s in list(loc_deos.columns):
    station_dict[s] = loc_deos[s]['station_id']

bad_sites = ['DNEM', 'DFHM','DWBD', 'DWWK', 'DSCR', 'DBUK1', 'DWCH', 'DTDF', 'DHOC', 'DCLY', 'DTLY', 'DCHI', 'DBKB', 'DWCC',
             'DPPN', 'DMTC', 'DMCB', 'DSJR', 'DFRE', 'DSND', 'DVIO', 'DADV', 'DPAR', 'DBBB', 'DSBY', 'DDAG', 'DGUM', 'DELN',
             'DMIL', 'DJCR', 'DPMH', 'DLEW', 'DNAS', 'DRBH', 'DIRL', 'DLNK', 'DSLB', 'DCPH']

for bad in bad_sites:
    try:
        del station_dict[bad]
    except:
        pass
rev_station_dict = dict(zip(station_dict.values(),station_dict.keys()))
nameDict = dict(zip(['Mean Daily Temp.','Max Daily Temp.','Min Daily Temp.','Heating Degree Days','Cooling Degree Days',
                     'Mean Wind Speed','Gage Precipitation (Daily)', 'Mean Daily Dew Point', 'Energy Density',
                     'Reference Evapotrans.', 'Growing Degree Days', 'Daily Avg RH', 'Daily Max RH', 'Daily Min RH',
                     'Daily Avg ST', 'Daily Max ST', 'Daily Min ST', 'Daily Avg VWC', 'Daily Max VWC', 'Daily Min VWC',
                     'Daily Solar', 'Mean Wind Direction', 'Peak Wind Gust Speed (Daily)', 'Daily Min WC', 'Daily Max HI'], 
                     ['meanTemp', 'maxTemp', 'minTemp', 'HDD', 'CDD', 'meanWS', 'dailyprecip', 'meanDP','energyDens',
                      'refET', 'GDD', 'meanRH', 'maxRH', 'minRH', 'meanST', 'maxST', 'minST', 'meanVWC', 'maxVWC', 'minVWC',
                      'meanSolar', 'meanWD', 'dailyGust', 'dailyMinWC', 'maxHI']))
fancyDict = dict(zip(list(nameDict.values()), ['Kelvin', 'Kelvin', 'Kelvin', ' ', ' ', 'm.s-1', 'mm', 'Kelvin', 'J.m-2',
                                               'mm.day-1', ' ', '%', '%', '%', 'Kelvin', 'Kelvin', 'Kelvin', ' ', ' ', ' ', 'J.m-2',
                                               'Rad', 'm.s-1', 'Kelvin', 'Kelvin']))

# create a dictionary for months
monthDict = dict(zip([1,2,3,4,5,6,7,8,9,10,11,12], ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']))
yearList = [2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020]
monthList = [1,2,3,4,5,6,7,8,9,10,11,12]

for yr in yearList:
    for mn in monthList:
        dayList = range(1,monthrange(yr,mn)[1] + 1)
        for dy in dayList:
            nowtime = datetime.strptime(str("{:04d}".format(yr) + "-" + "{:02d}".format(mn) + "-" + "{:02d}".format(dy)), "%Y-%m-%d")
            daytime = str("{:04d}".format(nowtime.year) + "-" + "{:02d}".format(nowtime.month) + "-" + "{:02d}".format(nowtime.day))
            
            if os.path.isfile(outPathNC + "/" + str(nowtime.year) + "/" + str("DEOS_agri_" + "{:04d}".format(nowtime.year) + "{:02d}".format(nowtime.month) + "{:02d}".format(nowtime.day) + ".nc")) == False:
                for var in nameDict:
                    # begin file loop to check if it exists 
                    lats=list()
                    lons=list()
                    varData = list()
                    for key in rev_station_dict:
                        stat_path = "http://128.175.28.202/deos_json/daily_summary/" + rev_station_dict[key] + "_" + monthDict[nowtime.month] + "-" + str(nowtime.year) + ".json"
                        #print(stat_path)
                        try:
                            agJson = pd.read_json(stat_path)
                            # use this for when we are real-time et.append(int(float(et_data[rev_station_dict[key]][str(str(nowtime.year) + "-" + str("{0:0=2d}".format(nowtime.month)) + "-" + str("{0:0=2d}".format(nowtime.day)))]['Reference Evapotrans.']['Value'])))
                            varData.append(round(float(agJson[rev_station_dict[key]][daytime][var]['Value']),4))
                            lats.append(loc_deos[rev_station_dict[key]]['latitude'])
                            lons.append(loc_deos[rev_station_dict[key]]['longitude'])
                        except:
                            pass
                    
                    if len(varData) != 0:
                        # add in four corners to expand the interpolated grid
                        lons = lons + list([-76.65,-76.65, -74.28,  -74.68])
                        lats = lats + list([38.0, 40.9, 38.0, 40.6])
                        try:
                            t1 = float(deos_data[date_deos][2321][var])
                        except:
                            t1 = np.nanmean(varData)
                        
                        try:
                            t2 = float(deos_data[date_deos][2980][var])
                        except:
                            t2 = np.nanmean(varData)
                        
                        try:
                            t3 = float(deos_data[date_deos][2304][var])
                        except:
                            t3 = np.nanmean(varData)
                        
                        try:
                            t4 = float(deos_data[date_deos][2983][var])
                        except:
                            t4 = np.nanmean(varData)
                        
                        varData = varData + list([t1,t2,t3,t4])
                        
                        lons=np.array(lons)
                        lats=np.array(lats)
                        varData = np.array(varData)
                        varData = varData.round(2)
                        
                        x = np.linspace(min(lons), max(lons), 750)
                        y = np.linspace(min(lats), max(lats), 750)
                        xi,yi = np.meshgrid(x,y)
                        # interpolate
                        #zi = griddata((lons,lats),temp,(xi,yi),method='cubic')
                        # try the idw interpolation scheme
                        xi, yi = xi.flatten(), yi.flatten()
                        
                        # Calculate IDW
                        zi = linear_rbf(lons,lats,varData,xi,yi)
                        zi=zi.reshape((len(x), len(y)))
                        
                        da = xr.DataArray(zi,dims=['lat', 'lon'],coords={'lon': x, 'lat' :y})
                        da.rio.set_crs("epsg:4326",inplace=True)
                        da.rio.set_spatial_dims('lon', 'lat', inplace=True)
                        da.rio.to_raster(temPath + nameDict[var] + "_temp.tif", overwrite=True)

                    # now that all variables have been interpolated as a spatial dataset, send to R so we can regrid and save as netCDF4
                    else:
                        print("no data is present")
                # once all the variables have processed, call et_regrid.R to rasterize them and place them into 1 nc file
                dtime = str("{:04d}".format(nowtime.year) + "{:02d}".format(nowtime.month) + "{:02d}".format(nowtime.day))
                print(dtime)
                regrid = "Rscript /home/sat_ops/deos/scripts/et_regrid.R " + str(dtime)
                os.system(regrid)
                os.system("rm " + temPath + "*_temp.tif")
                os.system("rm " + temPath + "*_temp.nc")
            else:
                print(str("DEOS_agri_" + "{:04d}".format(nowtime.year) + "{:02d}".format(nowtime.month) + "{:02d}".format(nowtime.day) + ".nc") + " already exists")



