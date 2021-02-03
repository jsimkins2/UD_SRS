import xarray as xr
import numpy as np
import metpy
from datetime import datetime
from xarray.backends import NetCDF4DataStore
from netCDF4 import Dataset, num2date
from siphon.catalog import TDSCatalog
from pyproj import Proj
import cartopy.crs as ccrs
from scipy.interpolate import griddata
import numpy.ma as ma

import os.path
from datetime import datetime, timedelta
from os import listdir
from os.path import isfile, join
import os.path
import os
import sys


def trim_data(lats, lons, ref, boundinglat, boundinglon):
    for i in range(lats.shape[0]):
        for j in range(lats.shape[1]):
            if (lats[i][j]>boundinglat[0] and lats[i][j]<boundinglat[1] and \
                lons[i][j]>boundinglon[0] and lons[i][j]<boundinglon[1]):
                pass
            else:
                # lats[i] = np.nan; lons[i] = np.nan
                ref[i][j] = np.nan
    return ref

min_lon = -76.35 #DEOS lons 
min_lat = 38.0 #DEOS lats
max_lat = 40.6 #DEOS lats 
max_lon =  -74.68 #DEOS lons 
boundinglat = [min_lat, max_lat]
boundinglon = [min_lon, max_lon]
nlon=500
nlat=500
    
datadir = "/Users/james/Downloads/hrrr_temp/"
filenames = [f for f in listdir(datadir) if isfile(join(datadir, f))]
hrrrdata = datadir + ''.join(([f for f in filenames if f[0:3]=='rep']))
print('starting the hrrr layer processing')

ds = xr.open_dataset(hrrrdata)
# parse the temperature at various heights
tsurf = ds.metpy.parse_cf('temperature_surface')
t850 = ds.metpy.parse_cf('temperature_850')
t925 = ds.metpy.parse_cf('temperature_925')
ht1000 = ds.metpy.parse_cf("height_1000")
ht500 = ds.metpy.parse_cf("height_500")
thick = ht500 - ht1000

hproj = ccrs.Geodetic()
hproj = Proj(hproj.proj4_init)

# trim the data to save space
lons, lats = np.meshgrid(tsurf['longitude'], tsurf['latitude'])
hrrr_t850 = trim_data(lats, lons, ma.getdata(t850), boundinglat, boundinglon)
hrrr_t925 = trim_data(lats, lons, ma.getdata(t925), boundinglat, boundinglon)
hrrr_tsurf = trim_data(lats, lons, ma.getdata(tsurf), boundinglat, boundinglon)
thick = trim_data(lats, lons, ma.getdata(thick), boundinglat, boundinglon)

# we have to ravel these for scipy interpolate
rav_lats = lats.ravel()
rav_lons = lons.ravel()
rav_t850 = hrrr_t850.ravel()
rav_t925 = hrrr_t925.ravel()
rav_tsurf = hrrr_tsurf.ravel()
rav_thick = thick.ravel()

#Grid Data using scipy interpolate
grid_lons = np.linspace(boundinglon[0],boundinglon[1],nlon)
grid_lats = np.linspace(boundinglat[0],boundinglat[1],nlat)
glon,glat = np.meshgrid(grid_lons,grid_lats)
grid850= griddata((rav_lons,rav_lats),rav_t850,(glon,glat),method='linear')
grid925 = griddata((rav_lons,rav_lats),rav_t925,(glon,glat),method='linear')
gridsurf = griddata((rav_lons,rav_lats),rav_tsurf,(glon,glat),method='linear')
gridthick = griddata((rav_lons,rav_lats),rav_thick,(glon,glat),method='linear')

np.save(datadir + 'grid850.npy', grid850)
np.save(datadir + 'grid925.npy', grid925)
np.save(datadir + 'gridsurf.npy', gridsurf)
np.save(datadir + 'gridthick.npy', gridthick)

print('finished making regridded numpy datasets')