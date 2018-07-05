import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from siphon.catalog import TDSCatalog
import urllib
from netCDF4 import Dataset, num2date
from matplotlib import ticker
import matplotlib as mpl
import pyart
from dateutil import tz
import time
from time import mktime
import os.path
import numpy as np
import matplotlib.image as image
from pyproj import Proj     
from datetime import datetime, timedelta
############# Initial Set Up ##################
workdir = "/home/sat_ops/goesR/sst/"
imgdir = "imgatlantic/"


# Go to the Unidata Thredds Server for the Current Day
nowdate = datetime.utcnow()
cat = TDSCatalog('http://thredds-jumbo.unidata.ucar.edu/thredds/catalog/satellite/goes16/GOES16/Products/SeaSurfaceTemperature/FullDisk/' + \
                  str(nowdate.year) + str("%02d"%nowdate.month) + str("%02d"%nowdate.day) + '/catalog.xml')
# run through the last 5 files for now and we'll see whether or not we've already created them or not
dataset_name = sorted(cat.datasets.keys())[-1]

print dataset_name
dataset = cat.datasets[dataset_name]
nc = dataset.remote_access()
list(nc.variables)
sst = nc.variables['SST'][:]
sst[sst==-1] = np.nan
add_seconds = nc.variables['t'][0]
DATE = datetime(2000, 1, 1, 12) + timedelta(seconds=int(add_seconds))
sat_h = nc.variables['goes_imager_projection'].perspective_point_height
sat_lon = nc.variables['goes_imager_projection'].longitude_of_projection_origin
sat_sweep = nc.variables['goes_imager_projection'].sweep_angle_axis
X = nc.variables['x'][:] * sat_h
Y = nc.variables['y'][:] * sat_h


# Configure EST/EDT depending on time of year
abi_time = DATE
from_zone = tz.gettz('UTC')
to_zone = tz.gettz('America/New_York')
utc = abi_time.replace(tzinfo=from_zone)
local = utc.astimezone(to_zone)

lt = time.localtime()
dst = lt.tm_isdst
lt = time.localtime()
dst = lt.tm_isdst
if dst == 0:
    et = "EDT"
else:
    et = "EST"


p = Proj(proj='geos', h=sat_h, lon_0=sat_lon, sweep=sat_sweep)
# Convert map points to latitude and longitude with the magic provided by Pyproj
XX, YY = np.meshgrid(X, Y)
lons, lats = p(XX, YY, inverse=True)
lats[np.isnan(sst)] = np.nan
lons[np.isnan(sst)] = np.nan



newmap = mH.pcolormesh(xH, yH, dBZ, cmap='jet')
newmap.set_array(None) # without this line, the linewidth is set to zero, but the RGB colorTuple is ignored. I don't know why.


mH = Basemap(resolution='i', projection='lcc', area_thresh=1500, \
            width=1800*3000, height=1060*3000, \
            lat_1=38.5, lat_2=38.5, \
            lat_0=38.5, lon_0=-97.5)
xH, yH = mH(lons, lats)

plt.figure(figsize=[16, 12], dpi=100)
newmap = mH.pcolormesh(xH, yH, sst, cmap='jet')
newmap.set_array(None)
mH.fillcontinents()
















# bbk 
# Go to the Unidata Thredds Server for the Current Day
nowdate = datetime.utcnow()
cat = TDSCatalog('http://thredds-jumbo.unidata.ucar.edu/thredds/catalog/satellite/goes16/GOES16/Products/SeaSurfaceTemperature/FullDisk/' + \
                  str(nowdate.year) + str("%02d"%nowdate.month) + str("%02d"%nowdate.day) + '/catalog.xml')
# run through the last 5 files for now and we'll see whether or not we've already created them or not
dataset_name = sorted(cat.datasets.keys())[-1]

print dataset_name
dataset = cat.datasets[dataset_name]
nc = dataset.remote_access()
list(nc.variables)
sst = np.array(nc.variables['SST'][:,:])
sst[np.isnan(sst)] = -1
sst = np.ma.array(sst)
sst[sst < 0] = np.ma.masked

    
    
add_seconds = nc.variables['t'][0]
DATE = datetime(2000, 1, 1, 12) + timedelta(seconds=int(add_seconds))
sat_h = nc.variables['goes_imager_projection'].perspective_point_height
sat_lon = nc.variables['goes_imager_projection'].longitude_of_projection_origin
sat_sweep = nc.variables['goes_imager_projection'].sweep_angle_axis
X = nc.variables['x'][:] * sat_h
Y = nc.variables['y'][:] * sat_h

mF = Basemap(projection='geos', lon_0='-89.5',
            llcrnrx=X.min(),llcrnry=Y.min(),
            urcrnrx=X.max(),urcrnry=Y.max())

XX, YY = np.meshgrid(X, Y)

p = Proj(proj='geos', h=sat_h, lon_0=sat_lon, sweep=sat_sweep)
lons, lats = p(XX, YY, inverse=True)

lats[np.isnan(sst)] = np.nan
lons[np.isnan(sst)] = np.nan

mH = Basemap(resolution='i', projection='lcc', \
            width=1800*3000, height=1060*3000, \
            lat_1=38.5, lat_2=38.5, \
            lat_0=38.5, lon_0=-97.5)

plt.figure(figsize=[16, 12], dpi=100)

# Plot the HRRR reflectivity
mH.pcolormesh(lons, lats,sst, latlon=True,
              cmap='jet',
              vmax=80, vmin=0)
mH.fillcontinents()

plt.figure(figsize=[10, 8])
mF.imshow(np.flipud(sst)) # Remember, "images" are upside down, so flip up/down
mF.drawcoastlines()
mF.drawcountries()
mF.fillcontinents()

















sst2 = sst

# have to add offset and scale factor 
sst2 = sst*0.00244163 + 180

X = nc.variables['x'][:] * sat_h
Y = nc.variables['y'][:] * sat_h
XX, YY = np.meshgrid(X, Y)
p = Proj(proj='geos', h=sat_h, lon_0=sat_lon, sweep=sat_sweep)
lons, lats = p(XX, YY, inverse=True)
lats[np.isnan(sst)] = np.nan
lons[np.isnan(sst)] = np.nan

m = Basemap(projection='mill', llcrnrlon=-120, urcrnrlon=0, llcrnrlat=0, urcrnrlat=60, resolution='c')

xH, yH = m(lons, lats)

fig = plt.figure()
cs = m.pcolormesh(xH, yH, sst2, shading='flat', cmap=plt.cm.jet)
m.drawcoastlines(linewidth=0.1)
m.fillcontinents(color='None')
m.drawmapboundary(fill_color='None')

cbar = m.colorbar(cs)
plt.show()






