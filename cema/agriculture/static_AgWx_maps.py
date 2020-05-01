# deos geopandas to be run at home 
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.colors import ListedColormap
import geopandas as gpd
import xarray as xr
import pyproj
from pyproj import Proj
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
import geopandas
from shapely.geometry import box, mapping
import matplotlib.colors as clr
from matplotlib.colors import BoundaryNorm
import matplotlib as mpl
import pandas as pd
import matplotlib.patheffects as path_effects
import time
import matplotlib.image as image
from datetime import datetime, timedelta

# declare paths
shapePaths = "/Users/James/Downloads/mapLayers/"
colorPaths = "/Users/James/Downloads/colorramps/"

def check_crs(crs):
    """Checks if the crs represents a valid grid, projection or ESPG string.
    Examples
    --------
    >>> p = check_crs('+units=m +init=epsg:26915')
    >>> p.srs
    '+units=m +init=epsg:26915 '
    >>> p = check_crs('wrong')
    >>> p is None
    True
    Returns
    -------
    A valid crs if possible, otherwise None
    """
    if isinstance(crs, pyproj.Proj) or isinstance(crs, Grid):
        out = crs
    elif isinstance(crs, dict) or isinstance(crs, string_types):
        try:
            out = pyproj.Proj(crs)
        except RuntimeError:
            try:
                out = pyproj.Proj(init=crs)
            except RuntimeError:
                out = None
    else:
        out = None
    return out
def proj_to_cartopy(proj):
    """Converts a pyproj.Proj to a cartopy.crs.Projection
    Parameters
    ----------
    proj: pyproj.Proj
        the projection to convert
    Returns
    -------
    a cartopy.crs.Projection object
    """

    import cartopy.crs as ccrs

    proj = check_crs(proj)

    #if proj.is_latlong():
        #return ccrs.PlateCarree()

    srs = proj.srs

    km_proj = {'lon_0': 'central_longitude',
               'lat_0': 'central_latitude',
               'x_0': 'false_easting',
               'y_0': 'false_northing',
               'k': 'scale_factor',
               'zone': 'zone',
               }
    km_globe = {'a': 'semimajor_axis',
                'b': 'semiminor_axis',
                }
    km_std = {'lat_1': 'lat_1',
              'lat_2': 'lat_2',
              }
    kw_proj = dict()
    kw_globe = dict()
    kw_std = dict()
    for s in srs.split('+'):
        s = s.split('=')
        if len(s) != 2:
            continue
        k = s[0].strip()
        v = s[1].strip()
        try:
            v = float(v)
        except:
            pass
        if k == 'proj':
            if v == 'tmerc':
                cl = ccrs.TransverseMercator
            if v == 'lcc':
                cl = ccrs.LambertConformal
            if v == 'merc':
                cl = ccrs.Mercator
            if v == 'utm':
                cl = ccrs.UTM
        if k in km_proj:
            kw_proj[km_proj[k]] = v
        if k in km_globe:
            kw_globe[km_globe[k]] = v
        if k in km_std:
            kw_std[km_std[k]] = v

    globe = None
    if kw_globe:
        globe = ccrs.Globe(**kw_globe)
    if kw_std:
        kw_proj['standard_parallels'] = (kw_std['lat_1'], kw_std['lat_2'])

    # mercatoooor
    if cl.__name__ == 'Mercator':
        kw_proj.pop('false_easting', None)
        kw_proj.pop('false_northing', None)

    return cl(globe=globe, **kw_proj)

# read in deos special shapefiles
deos_boundarys = gpd.read_file(shapePaths + 'deoscounties.shp')
bigdeos = gpd.read_file(shapePaths + 'TRISTATE_OVERVIEW.shp')
inland_bays = gpd.read_file(shapePaths + 'InlandBays.shp')
state_outline = gpd.read_file(shapePaths + 'tristateMultiaddedPACo.shp')

# create cartopy instance of obscure projection
c = Proj('+proj=tmerc +lat_0=38 +lon_0=-75.41666666666667 +k=0.999995 +x_0=200000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
oldproj = proj_to_cartopy(c)

# Declare bounds of the data
bounds=(-76.2,38.3,-74.85, 40.3)
# open agwx main dataset
agwx_main = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/DEOSAG.nc")
agwx_main = agwx_main.sel(latitude=slice(bounds[3], bounds[1]), longitude=slice(bounds[0],bounds[2]))

# open NCEP stage IV QC dataset
dsPrec = xr.open_dataset("http://thredds.demac.udel.edu/thredds/dodsC/NCEPIVQC.nc")
dsPrec = dsPrec.sel(lat=slice(bounds[1], bounds[3]), 
                    lon=slice(bounds[0],bounds[2]), 
                    time=slice(datetime.strptime("2014-01-01", "%Y-%m-%d"),
                              date.today()))

dsPrec = dsPrec.reindex(lat=list(reversed(dsPrec.lat)))
dsPrec = dsPrec.rename(name_dict= {'lat' : 'latitude'})
dsPrec = dsPrec.rename(name_dict= {'lon' : 'longitude'})
dsPrec = dsPrec.drop('crs')

sum_dict = dict(zip(['Heating Degree Days', 'Cooling Degree Days','DEOS Precip',  'Energy Density', 
                     'Reference Evapotranspiration', 'Growing Degree Days'],
                    ['HDD', 'CDD',  'dailyprecip', 'energyDens','refET', 'GDD']))
mean_dict = dict(zip(['Mean Temperature', 'Max Temperature', 'Min Temperature', 'Mean Wind Speed', 'Mean Dew Point',
                      'Mean Relative Humidity', 'Max Relative Humidity', 'Min Relative Humidity', 
                      'Mean Soil Temperature', 'Max Soil Temperature', 'Min Soil Temperature',
                      'Mean Volumetric Water Content','Max Volumetric Water Content', 'Min Volumetric Water Content',
                      'Mean Solar', 'Mean Wind Direction', 'Wind Gust', 'Min Wind Chill'],
                     ['meanTemp', 'maxTemp', 'minTemp', 'meanWS','meanDP','meanRH', 'maxRH', 'minRH',
                      'meanST', 'maxST', 'minST','meanVWC', 'maxVWC', 'minVWC','meanSolar', 'meanWD',
                      'dailyGust', 'dailyMinWC']))

datasets = list(sum_dict.keys()) + list(mean_dict.keys()) + ['NCEP Stage IV Precip', 'NCEP Stage IV Precip - DEOS RefET']
var = 'resET'
fig = plt.figure(figsize=(12,12))
ax = fig.add_subplot(111, projection=ccrs.Mercator())
ax.set_extent([-76.15, -75.03, 38.44, 40.26], crs=ccrs.PlateCarree())
for ind in range(0,len(bigdeos)):
        ax.add_geometries([bigdeos['geometry'][ind]], oldproj,
                      facecolor='silver', edgecolor='black')
im=ax.pcolormesh(cl['longitude'].values,cl['latitude'].values,cl.values,cmap=cmap,transform=ccrs.PlateCarree(),zorder=2)
for ind in range(0,len(deos_boundarys)):
    ax.add_geometries([deos_boundarys['geometry'][ind]], ccrs.PlateCarree(),
                      facecolor='none', edgecolor='black', zorder=3, linewidth=1.5)
for ind in range(0,len(inland_bays)):
    ax.add_geometries([inland_bays['geometry'][ind]], oldproj,
                      facecolor='white', edgecolor='black',zorder=3, linewidth=1.5)
ax.add_geometries([state_outline['geometry'][74]], oldproj, facecolor='none', edgecolor='black',zorder=3, linewidth=1.5)
ax.add_geometries([bigdeos['geometry'][121]], oldproj, facecolor='none', edgecolor='black',zorder=3, linewidth=1.5)
plt.title('refET')
#plt.text(-76.13, 38.503, fancyDict[var],horizontalalignment='left',color='white',weight='bold',size=9,zorder=30,transform=ccrs.PlateCarree())
#plt.text(-76.13, 38.473, deos_dateSTR,horizontalalignment='left',color='white',weight='bold',size=9,zorder=30,transform=ccrs.PlateCarree())
plt.colorbar(im)
im1 = image.imread("Downloads/maplayers/deos_logo.png")
plt.figimage(im1, 25, 40 ,zorder=30, alpha=1)
plt.savefig("Downloads/" + nameDict[var] + ".png")
