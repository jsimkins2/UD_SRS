# test script for snow sensor precip


# Precipitation Depiction with DEOS Sites
# By James Simkins with assistance from Dan Moore
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

import matplotlib as mpl
mpl.use('agg')
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
from metpy.plots import USCOUNTIES

from dateutil import tz
import time
import os.path
from datetime import datetime, timedelta
from os import listdir
from os.path import isfile, join
import os

import imageio

# monkey patch from nightmare that is january 10th
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh

# Some plotting work here with help from Dan Moore
def regrid_to_cartesian(radar, lon0, lat0):
    display = pyart.graph.RadarMapDisplay(radar)
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

def create_gif(workdir, imgdir, gifname):
    img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
    img_names = sorted(img_list)[-15:]
    imglen = len(img_names)
    images = []
    dur_vals = []
    for i in range(0,imglen -1):
        if i != imglen:
            dur_vals.append(.11)
    dur_vals.append(2)
    for i in img_names:
        input_file=imgdir + str(i)
        images.append(imageio.imread(input_file))
    imageio.mimsave(workdir + gifname, images, duration=dur_vals)








# try each radar location

try:
    site = 'KDOX'
    nowtime = datetime.utcnow().replace(second=0, microsecond=0)
    cat = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/nexrad/level2/' + site + '/' + str(nowtime.year) + str(nowtime.month).zfill(2) + str(nowtime.day).zfill(2) + '/catalog.xml')
    dataset = list(cat.datasets.values())[1]
    data_ymd = str(dataset).split('_')[2]
    data_hm = str(str(dataset).split('_')[3]).split('.')[0]
    data_datetime = datetime.strptime(data_ymd + data_hm,'%Y%m%d%H%M')
    time_diff = nowtime - data_datetime
    if time_diff.total_seconds() > 1600:
        raise ValueError('Datetimes too far apart, moving to next site')
except (HTTPError,ValueError):
    try:
        site = 'KLWX'
        nowtime = datetime.utcnow().replace(second=0, microsecond=0)
        cat = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/nexrad/level2/' + site + '/' + str(nowtime.year) + str(nowtime.month).zfill(2) + str(nowtime.day).zfill(2) + '/catalog.xml')
        dataset = list(cat.datasets.values())[1]
        data_ymd = str(dataset).split('_')[2]
        data_hm = str(str(dataset).split('_')[3]).split('.')[0]
        data_datetime = datetime.strptime(data_ymd + data_hm,'%Y%m%d%H%M')
        time_diff = nowtime - data_datetime
        if time_diff.total_seconds() > 1600:
            raise ValueError('Datetimes too far apart, moving to next site')
    except (HTTPError,ValueError):
        site = 'KDIX'
        nowtime = datetime.utcnow().replace(second=0, microsecond=0)
        cat = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/nexrad/level2/' + site + '/' + str(nowtime.year) + str(nowtime.month).zfill(2) + str(nowtime.day).zfill(2) + '/catalog.xml')
        dataset = list(cat.datasets.values())[1]
        data_ymd = str(dataset).split('_')[2]
        data_hm = str(str(dataset).split('_')[3]).split('.')[0]
        data_datetime = datetime.strptime(data_ymd + data_hm,'%Y%m%d%H%M')
        time_diff = nowtime - data_datetime
        if time_diff.total_seconds() > 1600:
            raise ValueError('Datetimes too far apart, moving to next site')




workdir = 'Users/james/downloads'
datadir = '/Users/james/Downloads/hrrr_temp/'
conv_thresh = 8 #dBZ
# create colormaps for each precip type

cmap_rain = LinearSegmentedColormap.from_list('mycmap', [(0,'honeydew'), (.2,'palegreen'),(.4,'darkseagreen'),(.6,'seagreen'), (.8,'green'), (1,'darkgreen')], N=50)
cmap_ice = LinearSegmentedColormap.from_list('mycmap', [(0,'mistyrose'), (.25,'pink'),(.5,'hotpink'), (.75,'deeppink'), (1,'mediumvioletred')], N=50)
cmap_sleet = LinearSegmentedColormap.from_list('mycmap', [(0,'Lavender'), (.33,'Violet'), (.66, 'DarkViolet'), (1,'purple')], N=50)
cmap_snow = LinearSegmentedColormap.from_list('mycmap', [(0,'powderblue'), (.2,'deepskyblue'), (.4,'dodgerblue'), (.6,'blue'), (.8,'mediumblue'),(1,'midnightblue')],N=50)


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

min_lon = -76.35 #DEOS lons 
min_lat = 38.0 #DEOS lats
max_lat = 40.6 #DEOS lats 
max_lon =  -74.68 #DEOS lons 
display = pyart.graph.RadarMapDisplay(radar)
lat0 = display.loc[0]
lon0 = display.loc[1]
boundinglat = [min_lat, max_lat]
boundinglon = [min_lon, max_lon]
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
# load in the temperatures data
grid850=np.load(datadir + 'grid850.npy')
grid925=np.load(datadir + 'grid925.npy')
gridsurf=np.load(datadir + 'gridsurf.npy')
gridthick=np.load(datadir + 'gridthick.npy')
# create a masked array for each precipitation type
rain = (gridsurf > 273.15) &  (np.isfinite(gref)) & (np.isfinite(grid850)) & (np.isfinite(grid925)) & (np.isfinite(gridsurf))
rain = np.ma.masked_array(gref, ~rain)
ice = (grid850 > 273.15) & (grid925 > 273.15) & (gridsurf < 273.15) &  (np.isfinite(gref)) & (np.isfinite(grid850)) & (np.isfinite(grid925)) & (np.isfinite(gridsurf))
ice = np.ma.masked_array(gref, ~ice)
sleet = (grid850 > 273.15) & (grid925 < 273.15) & (gridsurf < 273.15) &  (np.isfinite(gref)) & (np.isfinite(grid850)) & (np.isfinite(grid925)) & (np.isfinite(gridsurf))
sleet = np.ma.masked_array(gref, ~sleet)
snow = (grid850 < 273.15) & (grid925 < 273.15) & (gridsurf < 273.15) &  (np.isfinite(gref)) & (np.isfinite(grid850)) & (np.isfinite(grid925)) & (np.isfinite(gridsurf))
snow = np.ma.masked_array(gref, ~snow) 


snow_stations = ['DTLY','DCLY','DGRN','DPRC','DWCC','DHOC','DAGF','DGLW','DDMV','DBKB','DPPN','DDFS','DSMY','DWDS','DPAR',
                     'DFRE','DHAR','DDAG','DBNG','DSTK','DNAS','DELN','DLEW','DSEA','DBRG','DLAU','DNOT','DSGM','DTDF','DWBD',
                     'DWCH','DWHW','DWPK','DWWK']

import pandas as pd
import matplotlib.patheffects as path_effects
loc_deos = pd.read_json("http://128.175.28.202/deos_json/station_metadata.json")

stationLats=list()
stationLons=list()
for key in snow_stations:
    stationLats.append(loc_deos[key]['latitude'])
    stationLons.append(loc_deos[key]['longitude'])

fig=plt.figure(figsize=[12,10], dpi=90)
ax = plt.subplot(1,1,1, projection=ccrs.PlateCarree())
ax.set_extent((min_lon, max_lon, min_lat, 40.4))
ax.plot(lon0, lat0,color='k', linewidth=4, marker='o', transform=ccrs.PlateCarree())
im1 = ax.pcolormesh(glon, glat,rain,cmap=cmap_rain, vmin=0, vmax=50, transform = ccrs.PlateCarree())
im2 = ax.pcolormesh(glon, glat,ice,cmap=cmap_ice, vmin=0, vmax=50,transform = ccrs.PlateCarree())
im3 = ax.pcolormesh(glon, glat,sleet,cmap=cmap_sleet, vmin=0, vmax=50,transform = ccrs.PlateCarree())
im4 = ax.pcolormesh(glon, glat,snow,cmap=cmap_snow, vmin=0, vmax=50,transform = ccrs.PlateCarree())

# plot up the title 
title = 'CEMA Precipitation Type & 1000-500mb Thickness Lines '
timestr = local.strftime('%Y-%m-%d %H:%M ') + et
plt.title(title + "\n" + timestr, fontsize=12)

im5 = ax.contour(glon, glat, gridthick,levels=[5450, 5500,5550,5600,5650, 5700,5750, 5800], colors='k',linestyles='--', transform = ccrs.PlateCarree())
im6 = ax.contour(glon, glat, gridthick,levels = [5200, 5250, 5300, 5350,5400], colors='blue',linestyles='--',linewidths=2, transform = ccrs.PlateCarree())

for l in range(0,len(stationLons)):
    text = plt.text(stationLons[l],stationLats[l],snow_stations[l], size=7,weight='bold',verticalalignment='center',
    horizontalalignment='center',transform=ccrs.PlateCarree(),zorder=5)
    text.set_path_effects([path_effects.Stroke(linewidth=2.5, foreground='white'),path_effects.Normal()])
                        
                    
                    
# add contour labels
plt.clabel(im5, fmt='%1.0f')
plt.clabel(im6,fmt='%1.0f')

# add colorbars
cbaxes = fig.add_axes([0.75, 0.15, 0.02, 0.15]) 
cb1 = plt.colorbar(im1, cax = cbaxes, )  
cb1.ax.get_yaxis().labelpad = 6
cb1.ax.set_ylabel('Rainfall  [dBZ]', fontsize=12)
cb1.set_ticks([0, 10,20, 30, 40, 50])

cbaxes = fig.add_axes([0.75, 0.333, 0.02, 0.15]) 
cb2 = plt.colorbar(im2, cax = cbaxes)  
cb2.ax.get_yaxis().labelpad = 12
cb2.ax.set_ylabel('Ice  [dBZ]', fontsize=12)
cb2.set_ticks([0, 10,20, 30, 40, 50])

cbaxes = fig.add_axes([0.75, 0.513, 0.02, 0.15]) 
cb3 = plt.colorbar(im3, cax = cbaxes)  
cb3.ax.get_yaxis().labelpad = 12
cb3.ax.set_ylabel('Sleet  [dBZ]', fontsize=12)
cb3.set_ticks([0, 10,20, 30, 40, 50])

cbaxes = fig.add_axes([0.75, 0.7, 0.02, 0.15]) 
cb4 = plt.colorbar(im4, cax = cbaxes)  
cb4.ax.get_yaxis().labelpad = 12
cb4.ax.set_ylabel('Snowfall  [dBZ]', fontsize=12)
cb4.set_ticks([0, 10,20, 30, 40, 50])

# plot coasts/states/counties/lakes
request = cimgt.GoogleTiles(url="https://cartodb-basemaps-d.global.ssl.fastly.net/light_nolabels/{z}/{x}/{y}.png")
ax.add_image(request, 7, zorder=0, interpolation='none')
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '10m',edgecolor='black', facecolor='none',linewidth=1.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '10m',edgecolor='black', facecolor='none',linewidth=1.5))
ax.add_feature(USCOUNTIES.with_scale('500k'), linewidth=0.5, edgecolor="black")
    
plt.savefig('/Users/james/Downloads/test.png', bbox_inches='tight',dpi=90)
plt.close()


