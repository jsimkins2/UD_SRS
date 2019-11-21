import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import BoundaryNorm
import matplotlib.pyplot as plt
import numpy as np

from metpy.cbook import get_test_data
from metpy.interpolate import (interpolate_to_grid, remove_nan_observations,
                               remove_repeat_coordinates)

from metpy.plots import ctables
from scipy.interpolate import griddata

import pandas as pd
import numpy as np 
import scipy
# mask out the field
from metpy.plots import USCOUNTIES

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


# plot
levels = list(range(-30, 120, 1))
cmap = ctables.registry.get_colortable('rainbow')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)


fig = plt.figure(figsize=(12,12))
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
ax.add_feature(USCOUNTIES.with_scale('500k'))
plt.pcolormesh(xi,yi,zi, cmap=cmap, norm=norm,transform=ccrs.PlateCarree())
plt.plot(lons,lats,'k.')
plt.xlabel('xi',fontsize=16)
plt.ylabel('yi',fontsize=16)



from cartopy.feature import NaturalEarthFeature, LAND, COASTLINE

ocean = NaturalEarthFeature(category='physical', name='ocean',
                            scale='50m', facecolor='white')


import cartopy.io.shapereader as shpreader

reader = shpreader.Reader('Downloads/us_counties_500k.shp')
counties = reader.records()

fig = plt.figure(figsize=(12,12))
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
ax.set_extent([-76.1, -75.02, 38.35, 40.25])
plt.pcolormesh(xi,yi,zi, cmap=cmap, norm=norm,transform=ccrs.PlateCarree())
plt.plot(lons,lats,'k.')
ax.add_feature(USCOUNTIES.with_scale('500k'))

for country in counties:
    if country.attributes['STATEFP'] == '10' or country.attributes['STATEFP'] == '24' or country.attributes['STATEFP'] == '42' or country.attributes['STATEFP'] == '34' or country.attributes['STATEFP'] == '51':
        if country.attributes['NAME'] == 'Kent' or country.attributes['NAME'] == 'Sussex' or country.attributes['NAME'] == 'New Castle' or country.attributes['NAME'] == 'Chester':
            ax.add_geometries(country.geometry, ccrs.PlateCarree(),
                              facecolor='none')
        else:
            ax.add_geometries(country.geometry, ccrs.PlateCarree(),
                              facecolor='gray')




                      
