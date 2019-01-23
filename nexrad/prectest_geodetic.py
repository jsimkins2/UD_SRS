import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import patheffects, ticker
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.image as image
from matplotlib.colors import LinearSegmentedColormap

from siphon.catalog import TDSCatalog
from siphon.radarserver import RadarServer
import urllib
from netCDF4 import Dataset, num2date
import xarray as xr
from xarray.backends import NetCDF4DataStore
import numpy as np
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
import cartopy.crs as ccrs
import metpy
from metpy.plots import colortables
import pyart
from scipy.interpolate import griddata
import numpy.ma as ma
from pyproj import Proj, transform
from metpy.plots import USCOUNTIES

from dateutil import tz
import time
from time import mktime
import os.path
from datetime import datetime, timedelta
from os import listdir
from os.path import isfile, join
import calendar


# grab the HRRR data 

# grab the radar data 
site = 'KCCX'
#site = 'KDOX'
cat = TDSCatalog('http://thredds-aws.unidata.ucar.edu/thredds/radarServer/catalog.xml')
rs = RadarServer(cat.catalog_refs['S3 NEXRAD Level II'].href)
rs.stations[site]
query = rs.query()
query.stations(site).time(datetime.utcnow())
rs.validate_query(query)
catalog = rs.get_catalog(query)
dataset = list(catalog.datasets.values())[0]
radar = pyart.io.read_nexrad_cdm(dataset.access_urls['OPENDAP'])

my_gf = pyart.filters.GateFilter(radar)
my_gf.exclude_below('reflectivity', 12)
my_ds_gf = pyart.correct.despeckle_field(radar, 'reflectivity', gatefilter=my_gf)


timestamp = radar.time['units'].split(' ')[-1].split('T')
timestamp = timestamp[0] + ' ' + timestamp[1][:-1]
timestamp = datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S') + timedelta(minutes=4)
from_zone = tz.gettz('UTC')
to_zone = tz.gettz('America/New_York')
utc = timestamp.replace(tzinfo=from_zone)
local = utc.astimezone(to_zone)
lt = time.localtime()
dst = lt.tm_isdst
lt = time.localtime()
dst = lt.tm_isdst
if dst == 0:
  et = "EST"
else:
  et = "EDT"

lats = radar.gate_latitude
lons = radar.gate_longitude
min_lon = -80#-78.24357604980469 #lons['data'].min() + 2.5
min_lat = 39#36.69588851928711 #lats['data'].min() + 2
max_lat = 43#40.95521545410156 #lats['data'].max() - 2
max_lon = -76#-72.63585662841797 #lons['data'].max() - 2.5
display = pyart.graph.RadarMapDisplayCartopy(radar)
lat0 = display.loc[0]
lon0 = display.loc[1]



def regrid_to_cartesian(radar, lon0, lat0):
    display = pyart.graph.RadarMapDisplayCartopy(radar)
    x,y = display._get_x_y(0,True,None)
    x = x*1000; y = y*1000
    lambert_aea = {'proj': 'laea',
      'lat_0':lat0, 
      'lon_0':lon0, 
      'x_0':0., 
      'y_0':0.,
      'ellps': 'WGS84',
      'datum': 'WGS84',
      'R':6378137.0}
    lons, lats = pyart.core.cartesian_to_geographic(x,y, projparams = lambert_aea)
    #Otherwise, lats and lons have 1 extra data point.
    lats = lats[0:720, 0:1832]
    lons = lons[0:720, 0:1832]
    return lons,lats

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

def trim_rad_data(lats, lons, ref, boundinglat=0, boundinglon=0):
    for i in range(lats.shape[0]):
        for j in range(lats.shape[1]):
            if (lats[i][j]>boundinglat[0] and lats[i][j]<boundinglat[1] and \
                lons[i][j]>boundinglon[0] and lons[i][j]<boundinglon[1]):
                pass
            else:
                # lats[i] = np.nan; lons[i] = np.nan
                ref[i][j] = np.nan
                
            if ref[i][j]>conv_thresh:
                pass
            else:
                ref[i][j] = np.nan
    return ref


boundinglat = [min_lat, max_lat]
boundinglon = [min_lon, max_lon]
conv_thresh = 10.0#dBZ
lons, lats = regrid_to_cartesian(radar, lon0, lat0)
my_ref = radar.get_field(0, 'reflectivity')
unmasked = ma.getdata(my_ref)
my_ref = trim_rad_data(lats, lons, unmasked, boundinglat, boundinglon)
rav_lats = lats.ravel()
rav_lons = lons.ravel()
rav_ref = my_ref.ravel()

nlon = 500; nlat = 500 #can be changed, but seems good.
#Grid Data using matplotlib
grid_lons = np.linspace(boundinglon[0],boundinglon[1],nlon)
grid_lats = np.linspace(boundinglat[0],boundinglat[1],nlat)
glon,glat = np.meshgrid(grid_lons,grid_lats)

# Interpolate data onto grid using linear interpolation
gref = griddata((rav_lons,rav_lats),rav_ref,(glon,glat),method='linear')


ds = xr.open_dataset('Downloads/HRRR_CONUS_2p5km_20190123_1400.ncraster.nc')
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

# create a masked array for each precipitation type
rain = (gridsurf > 273.15) &  (np.isfinite(gref)) & (np.isfinite(grid850)) & (np.isfinite(grid925)) & (np.isfinite(gridsurf))
rain = np.ma.masked_array(gref, ~rain)
ice = (grid850 > 273.15) & (grid925 > 273.15) & (gridsurf < 273.15) &  (np.isfinite(gref)) & (np.isfinite(grid850)) & (np.isfinite(grid925)) & (np.isfinite(gridsurf))
ice = np.ma.masked_array(gref, ~ice)
sleet = (grid850 > 273.15) & (grid925 < 273.15) & (gridsurf < 273.15) &  (np.isfinite(gref)) & (np.isfinite(grid850)) & (np.isfinite(grid925)) & (np.isfinite(gridsurf))
sleet = np.ma.masked_array(gref, ~sleet)
snow = (grid850 < 273.15) & (grid925 < 273.15) & (gridsurf < 273.15) &  (np.isfinite(gref)) & (np.isfinite(grid850)) & (np.isfinite(grid925)) & (np.isfinite(gridsurf))
snow = np.ma.masked_array(gref, ~snow) 

# create colormaps for each precip type
cmap_rain = LinearSegmentedColormap.from_list('mycmap', ['palegreen', 'springgreen','darkseagreen','mediumseagreen','seagreen', 'green', 'darkgreen'], N=20)
cmap_ice = LinearSegmentedColormap.from_list('mycmap', ['lightpink','Pink', 'HotPink', 'deeppink'], N=20)
cmap_sleet = LinearSegmentedColormap.from_list('mycmap', ['Lavender', 'Violet', 'DarkViolet', 'purple'], N=20)
#cmap_snow = mpl.cm.Blues
cmap_snow = LinearSegmentedColormap.from_list('mycmap', ['powderblue', 'deepskyblue', 'dodgerblue', 'blue', 'mediumblue','midnightblue'],N=20)

# Now let's plot it all up!
fig = plt.figure(figsize=[11,11])
ax = plt.subplot(1,1,1, projection=ccrs.Mercator())
ax.set_extent((min_lon, max_lon, min_lat, max_lat))
ax.plot(lon0, lat0,
         color='k', linewidth=4, marker='o',
         transform=ccrs.PlateCarree())
im1 = ax.pcolormesh(glon, glat,rain,cmap=cmap_rain, vmin=0, vmax=50, transform = ccrs.PlateCarree())
im2 = ax.pcolormesh(glon, glat,ice,cmap=cmap_ice, vmin=0, vmax=50,transform = ccrs.PlateCarree())
im3 = ax.pcolormesh(glon, glat,sleet,cmap=cmap_sleet, vmin=0, vmax=50,transform = ccrs.PlateCarree())
im4 = ax.pcolormesh(glon, glat,snow,cmap=cmap_snow, vmin=0, vmax=50,transform = ccrs.PlateCarree())
title = 'CEMA Precipitation Depiction & 1000-500mb Thickness Lines '
timestr = local.strftime('%Y-%m-%d %H:%M ') + et
plt.title(title + "\n" + timestr, fontsize=18)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                edgecolor='black', facecolor='none',linewidth=1.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                edgecolor='black', facecolor='none',linewidth=1.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                edgecolor='black', facecolor='none',linewidth=1.5))
ax.add_feature(USCOUNTIES.with_scale('5m'), linewidth=0.5)
im5 = ax.contour(glon, glat, gridthick,levels=[5460, 5520,5550,5600,5650, 5700,5750, 5800], colors='k',linestyles='--', transform = ccrs.PlateCarree())
CS = ax.contour(glon, glat, gridthick,levels = [5200, 5250, 5300, 5350,5400],
                 colors='blue',linestyles='--',linewidths=2, transform = ccrs.PlateCarree())
plt.clabel(im5, fmt='%1.0f', transform=tempiso.metpy.cartopy_crs)
plt.clabel(CS,fmt='%1.0f', transform=tempiso.metpy.cartopy_crs)

cbaxes = fig.add_axes([0.93, 0.15, 0.02, 0.15]) 
cb1 = plt.colorbar(im1, cax = cbaxes, )  
cb1.ax.get_yaxis().labelpad = 12
cb1.ax.set_ylabel('Rainfall  [dBZ]', fontsize=12)
cb1.set_ticks([0, 10,20, 30, 40, 50])

cbaxes = fig.add_axes([0.93, 0.333, 0.02, 0.15]) 
cb2 = plt.colorbar(im2, cax = cbaxes)  
cb2.ax.get_yaxis().labelpad = 12
cb2.ax.set_ylabel('Ice  [dBZ]', fontsize=12)
cb2.set_ticks([0, 10,20, 30, 40, 50])

cbaxes = fig.add_axes([0.93, 0.513, 0.02, 0.15]) 
cb3 = plt.colorbar(im3, cax = cbaxes)  
cb3.ax.get_yaxis().labelpad = 12
cb3.ax.set_ylabel('Sleet  [dBZ]', fontsize=12)
cb3.set_ticks([0, 10,20, 30, 40, 50])

cbaxes = fig.add_axes([0.93, 0.7, 0.02, 0.15]) 
cb4 = plt.colorbar(im4, cax = cbaxes)  
cb4.ax.get_yaxis().labelpad = 12
cb4.ax.set_ylabel('Snowfall  [dBZ]', fontsize=12)
cb4.set_ticks([0, 10,20, 30, 40, 50])

plt.savefig("Downloads/preliminary_precipitation_depiction.png", dpi=100, bbox_inches='tight')
