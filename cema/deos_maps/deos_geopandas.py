import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.colors import ListedColormap
import geopandas as gpd
from shapely.geometry import box
from scipy.interpolate import griddata

## Load in other functions from metpy that couldn't be laded 
def remove_nan_observations(x, y, z):
    r"""Remove all x, y, and z where z is nan.
    Will not destroy original values.
    Parameters
    ----------
    x: array_like
        x coordinate
    y: array_like
        y coordinate
    z: array_like
        observation value
    Returns
    -------
    x, y, z
        List of coordinate observation pairs without
        nan valued observations.
    """
    x_ = x[~np.isnan(z)]
    y_ = y[~np.isnan(z)]
    z_ = z[~np.isnan(z)]

    return x_, y_, z_






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






import rioxarray
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

    if proj.is_latlong():
        return ccrs.PlateCarree()

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


da = xr.DataArray(zi,dims=['lat', 'lon'],coords={'lon': x, 'lat' :y})

da.rio.set_crs("epsg:4326")
da.attrs['units'] = 'Fahrenheit'
da.attrs['standard_name'] = 'Temperature'
da.rio.set_spatial_dims('lon', 'lat')
da.rio.to_raster('/Users/james/Downloads/test.tif')

xds = rioxarray.open_rasterio("/Users/james/Downloads/test.tif")



# read in deos special shapefiles
deos_boundarys = gpd.read_file('Downloads/mapLayers/deoscounties.shp')
bigdeos = gpd.read_file('Downloads/mapLayers/TRISTATE_OVERVIEW.shp')
inland_bays = gpd.read_file('Downloads/mapLayers/InlandBays.shp')

# create cartopy instance of obscure projection
c = Proj('+proj=tmerc +lat_0=38 +lon_0=-75.41666666666667 +k=0.999995 +x_0=200000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
oldproj = proj_to_cartopy(c)

# clip the interpolated data based on the shapefiles
clipped = xds.rio.clip(deos_boundarys.geometry.apply(mapping), xds.rio.crs, drop=True)
cl = clipped.rio.clip(inland_bays.geometry.apply(mapping), oldproj.proj4_init, drop=False, invert=True)



### Call the function make_cmap which returns your colormap

txt_cmap =  pd.read_csv('Downloads/colorramps/at_ramp.txt', header=None,names=['bound', 'r', 'g', 'b', 'a'],delimiter=' ')
raw_rgb = []
for i in range(0,len(txt_cmap)):
    raw_rgb.append(tuple([txt_cmap['r'][i], txt_cmap['g'][i], txt_cmap['b'][i]]))

cmap = clr.LinearSegmentedColormap.from_list('custom blue', raw_rgb, N=700)

fig = plt.figure(figsize=(12,12))
ax = fig.add_subplot(111, projection=ccrs.Mercator())
ax.set_extent([-76.1, -75.02, 38.35, 40.3], crs=ccrs.PlateCarree())
for ind in range(0,len(bigdeos)):
        ax.add_geometries([bigdeos['geometry'][ind]], oldproj,
                      facecolor='gray', edgecolor='black')
im=ax.pcolormesh(cl['x'].values,cl['y'].values,cl.values[0], cmap=cmap, vmin=-30, vmax=120,transform=ccrs.PlateCarree(),zorder=2)
plt.plot(lons,lats,'k.', transform=ccrs.PlateCarree())
for ind in range(0,len(deos_boundarys)):
    ax.add_geometries([deos_boundarys['geometry'][ind]], ccrs.PlateCarree(),
                      facecolor='none', edgecolor='black', zorder=3, linewidth=1.5)
for ind in range(0,len(inland_bays)):
        ax.add_geometries([inland_bays['geometry'][ind]], oldproj,
                      facecolor='white', edgecolor='black',zorder=3, linewidth=1.5)
plt.colorbar(im)


