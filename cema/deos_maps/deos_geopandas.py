import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.colors import ListedColormap
import geopandas as gpd
# Load the box module from shapely to create box objects
from shapely.geometry import box
import earthpy as et
from earthpy import clip as cl
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

# target grid to interpolate to
temp = list()
# Begin the for loop to place a location on a map
for s in station_id:
    if str(s) in list(station_dict.values()):
        ID = str(s)
        try:
            temp.append(float(deos_data[date_deos][s]['Air Temperature']))
        except:
            temp.append(np.nan)

rev_station_dict = dict(zip(station_dict.values(),station_dict.keys()))

lats=list()
lons=list()
for key in station_id:
    # grag latitude and longitudes location of station
    lats.append(loc_deos[rev_station_dict[str(key)]]['latitude'])
    lons.append(loc_deos[rev_station_dict[str(key)]]['longitude'])

# add in four corners to expand the interpolated grid
lons = lons + list([-76.15,-76.15, -74.98,  -74.98])
lats = lats + list([38.3, 40.3, 38.3, 40.3])

# just in case one of the temperature gauges is out
try:
    t1 = float(deos_data[date_deos][2321]['Air Temperature'])
except:
    t1 = np.nanmean(temp)

try:
    t2 = float(deos_data[date_deos][2980]['Air Temperature'])
except:
    t2 = np.nanmean(temp)

try:
    t3 = float(deos_data[date_deos][2304]['Air Temperature'])
except:
    t3 = np.nanmean(temp)

try:
    t4 = float(deos_data[date_deos][2983]['Air Temperature'])
except:
    t4 = np.nanmean(temp)

temp = temp + list([t1,t2,t3,t4])

lons=np.array(lons)
lats=np.array(lats)
temp = np.array(temp)
temp = temp - 273.15
temp = (temp*(9/5)) + 32

lons,lats, temp = remove_nan_observations(lons,lats, temp)

x = np.linspace(min(lons), max(lons), 500)
y = np.linspace(min(lats), max(lats), 500)
xi,yi = np.meshgrid(x,y)
# interpolate
zi = griddata((lons,lats),temp,(xi,yi),method='cubic')

# Import all of your data at the top of your notebook to keep things organized.
deos_boundarys = gpd.read_file(
    'Downloads/mapLayers/deoscounties.shp')

fig = plt.figure(figsize=(12,12))
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
ax.set_extent([-76.1, -75.02, 38.35, 40.25])
ax.pcolormesh(xi,yi,zi, cmap=cmap, norm=norm,transform=ccrs.PlateCarree(), clip_path=(clip))
plt.plot(lons,lats,'k.')
#ax.add_feature(USCOUNTIES.with_scale('500k'))

for county in deos_boundarys['NAME']:
    ind = deos_boundarys[deos_boundarys['NAME'] == county].index[0]
    ax.add_geometries([deos_boundarys['geometry'][ind]], ccrs.PlateCarree(),facecolor='none', edgecolor='black')
else:
    ax.add_geometries(facecolor='gray')

ind = can[can['name'] == state].index[0]

ax.add_geometries([can['geometry'][ind]], ccrs.PlateCarree(),
                  facecolor=facecolor, edgecolor=edgecolor)