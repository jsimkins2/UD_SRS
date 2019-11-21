import pandas as pd
import numpy as np 
import scipy

def nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
        
deos_data = pd.read_json("http://128.175.28.202/deos_json/map_data.json")
loc_deos = pd.read_json("http://128.175.28.202/deos_json/station_metadata.json")

# index is the station nunumbers 
station_id = list(deos_data.index)

# need to find time value
date_deos = deos_data.columns[0]

# run through each met variable
varnames = list(deos_data[date_deos][2302].keys())

# create a dict of station IDs with the name of the station
station_dict = {}
for s in list(loc_deos.columns):
    station_dict[s] = loc_deos[s]['station_id']

temp = {}
# Begin the for loop to place a location on a map
for s in station_id:
    if str(s) in list(station_dict.values()):
        ID = str(s)
        try:
            temp[ID] = deos_data[date_deos][s]['Air Temperature']
        except:
            pass


# maybe create an empty numpy array with associated lat lon markers and place the values in the appropriate grid cell
# once we have a numpy with a bunch of blank values we can interpolate



min_lon = -78.24357604980469 #lons['data'].min() + 2.5
min_lat = 36.69588851928711 #lats['data'].min() + 2
max_lat = 40.95521545410156 #lats['data'].max() - 2
max_lon = -72.63585662841797 #lons['data'].max() - 2.5

boundinglat = [min_lat, max_lat]
boundinglon = [min_lon, max_lon]
nlon = 500; nlat = 500 #can be changed, but seems good.
#Grid Data using matplotlib
grid_lons = np.linspace(boundinglon[0],boundinglon[1],nlon, dtype=float)
grid_lats = np.linspace(boundinglat[0],boundinglat[1],nlat, dtype=float)
glon,glat = np.meshgrid(grid_lons,grid_lats)

var_array = np.zeros((nlon,nlat))
var_array.fill(np.nan)

rev_station_dict = dict(zip(station_dict.values(),station_dict.keys()))
for key,val in temp.items():
    # grag latitude and longitudes location of station
    lat = loc_deos[rev_station_dict[key]]['latitude']
    lon = loc_deos[rev_station_dict[key]]['longitude']
    # find nearest index value of lat/lon, plug in value into our empy var array
    var_array[nearest(grid_lons, lon),nearest(grid_lats, lat)] = float(val)
    
    
# Interpolate data onto grid using linear interpolation
gref = griddata((glon,glat),list(float(temp.values())),(glon,glat),method='linear')
from scipy import interpolate
gref = scipy.interpolate.interp2d(grid_lons, grid_lats,var_array)
