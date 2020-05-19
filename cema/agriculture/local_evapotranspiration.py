
# deos geopandas to be run at home 
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.colors import ListedColormap
import geopandas as gpd
from scipy.interpolate import griddata
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
## define custom functions
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
def make_cmap(colors, position=None, bit=False):
    '''
    make_cmap takes a list of tuples which contain RGB values. The RGB
    values may either be in 8-bit [0 to 255] (in which bit must be set to
    True when called) or arithmetic [0 to 1] (default). make_cmap returns
    a cmap with equally spaced colors.
    Arrange your tuples so that the first color is the lowest value for the
    colorbar and the last is the highest.
    position contains values from 0 to 1 to dictate the location of each color.
    '''
    import matplotlib as mpl
    import numpy as np
    bit_rgb = np.linspace(0,1,256)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256*4)
    return cmap
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


# read in deos special shapefiles
deos_boundarys = gpd.read_file(shapePaths + 'deoscounties.shp')
bigdeos = gpd.read_file(shapePaths + 'TRISTATE_OVERVIEW.shp')
inland_bays = gpd.read_file(shapePaths + 'InlandBays.shp')
state_outline = gpd.read_file(shapePaths + 'tristateMultiaddedPACo.shp')

# create cartopy instance of obscure projection
c = Proj('+proj=tmerc +lat_0=38 +lon_0=-75.41666666666667 +k=0.999995 +x_0=200000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
oldproj = proj_to_cartopy(c)

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

nowtime = datetime.utcnow()
monthDict = dict(zip([1,2,3,4,5,6,7,8,9,10,11,12], ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'June', 'July', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']))
# create a dict of station IDs with the name of the station
station_dict = {}

for s in list(loc_deos.columns):
    station_dict[s] = loc_deos[s]['station_id']


rev_station_dict = dict(zip(station_dict.values(),station_dict.keys()))



nameDict = dict(zip(['Gage Precipitation (60)','Air Temperature','Dew Point Temperature','Wind Speed','Wind Direction','Barometric Pressure','Relative humidity', '24 Hour Precipitation', 'Peak Wind Gust Speed (60)'], ['precip', 'airT', 'dewP', 'wspeed', 'wdir', 'pres', 'rh', 'dailyprecip', 'gust']))
fancyDict = dict(zip(list(nameDict.keys()), ['1-hr Rain (in)', 'Air Temperature (F)', 'Dew Point (F)', '5-min Wind Speed', 'Wind', 'Pressure (mb)', 'Relative Humidity (%)', '24-hr Rain (in)', '1-hr Peak Wind Gust (mph)']))


lats=list()
lons=list()
et = list()
for key in rev_station_dict:
    stat_path = "http://128.175.28.202/deos_json/daily_summary/" + rev_station_dict[key] + "_" + monthDict[nowtime.month] + "-" + str(nowtime.year) + ".json"
    print(stat_path)
    try:
        et_data = pd.read_json(stat_path)
        # use this for when we are real-time et.append(int(float(et_data[rev_station_dict[key]][str(str(nowtime.year) + "-" + str("{0:0=2d}".format(nowtime.month)) + "-" + str("{0:0=2d}".format(nowtime.day)))]['Reference Evapotrans.']['Value'])))
        et.append(round(float(et_data[rev_station_dict[key]][str('2020-02-17')]['Reference Evapotrans.']['Value']),4))
        lats.append(loc_deos[rev_station_dict[key]]['latitude'])
        lons.append(loc_deos[rev_station_dict[key]]['longitude'])
    except:
        pass


# add in four corners to expand the interpolated grid
lons = lons + list([-76.15,-76.15, -74.98,  -74.98])
lats = lats + list([38.3, 40.3, 38.3, 40.3])
try:
    t1 = float(deos_data[date_deos][2321][var])
except:
    t1 = np.nanmean(et)

try:
    t2 = float(deos_data[date_deos][2980][var])
except:
    t2 = np.nanmean(et)

try:
    t3 = float(deos_data[date_deos][2304][var])
except:
    t3 = np.nanmean(et)

try:
    t4 = float(deos_data[date_deos][2983][var])
except:
    t4 = np.nanmean(et)

et = et + list([t1,t2,t3,t4])

lons=np.array(lons)
lats=np.array(lats)
et = np.array(et)

x = np.linspace(min(lons), max(lons), 750)
y = np.linspace(min(lats), max(lats), 750)
xi,yi = np.meshgrid(x,y)
# interpolate
#zi = griddata((lons,lats),temp,(xi,yi),method='cubic')
# try the idw interpolation scheme
xi, yi = xi.flatten(), yi.flatten()

# Calculate IDW
zi = linear_rbf(lons,lats,et,xi,yi)
zi=zi.reshape((len(x), len(y)))


da = xr.DataArray(zi,dims=['lat', 'lon'],coords={'lon': x, 'lat' :y})
da.rio.set_crs("epsg:4326")
da.attrs['units'] = 'mm.day-1'
da.attrs['standard_name'] = 'Reference ET'
da.rio.set_spatial_dims('lon', 'lat')
da.rio.to_raster('Downloads/ET.tif', overwrite=True)

xds = rioxarray.open_rasterio('Downloads/ET.tif')
# clip the interpolated data based on the shapefiles
clipped = xds.rio.clip(deos_boundarys.geometry.apply(mapping), xds.rio.crs, drop=True)
cl_et = clipped.rio.clip(inland_bays.geometry.apply(mapping), oldproj.proj4_init, drop=False, invert=True)


# now reading in ncep stage IV precip

ncepPrec = xr.open_dataset('http://thredds.demac.udel.edu/thredds/dodsC/ncep_stage_iv.nc')

ncepPrec = ncepPrec.Precipitation_Flux.sel(time=str('2020-02-17'), lat=slice(36,42), lon=slice(-78,-72))
ncepPrec = ncepPrec.sum('time') # should this be mean or sum? 

prec = list()

for val in range(0,len(lons)):
    prec.append(round(float(ncepPrec.sel(lat=lats[val], lon=lons[val], method='nearest').values),4))

zi = linear_rbf(lons,lats,prec,xi,yi)
zi=zi.reshape((len(x), len(y)))

da = xr.DataArray(zi,dims=['lat', 'lon'],coords={'lon': x, 'lat' :y})
da.rio.set_crs("epsg:4326")
da.attrs['units'] = 'mm.day-1'
da.attrs['standard_name'] = 'Precipitation_Flux'
da.rio.set_spatial_dims('lon', 'lat')
da.rio.to_raster('Downloads/PREC.tif', overwrite=True)

xds = rioxarray.open_rasterio('Downloads/PREC.tif')
# clip the interpolated data based on the shapefiles
clipped = xds.rio.clip(deos_boundarys.geometry.apply(mapping), xds.rio.crs, drop=True)
cl_prec = clipped.rio.clip(inland_bays.geometry.apply(mapping), oldproj.proj4_init, drop=False, invert=True)

cl = cl_prec - cl_et
temp = (prec - et) / 10 # convert from mm to cm

#### custom colormap  
import matplotlib as mpl
startcolor = '#8B4513'
midcolor = '#FFFFFF'
endcolor = '#008000'
own_cmap1 = mpl.colors.LinearSegmentedColormap.from_list( 'own2', [startcolor, midcolor, endcolor] )


fig = plt.figure(figsize=(12,12))
ax = fig.add_subplot(111, projection=ccrs.Mercator())
ax.set_extent([-76.15, -75.03, 38.44, 40.26], crs=ccrs.PlateCarree())
for ind in range(0,len(bigdeos)):
        ax.add_geometries([bigdeos['geometry'][ind]], oldproj,
                      facecolor='silver', edgecolor='black')
im=ax.pcolormesh(cl['x'].values,cl['y'].values,cl.values[0],cmap=own_cmap1,transform=ccrs.PlateCarree(),zorder=2, vmin=-10, vmax=10)

for l in range(0,len(lons)):
    if lons[l] != -76.15 and lons[l] != -74.98 and lons[l] != -75.062685 and lons[l] != -75.118033 and lons[l] != -75.247235 and lons[l] != -75.640685 and lons[l] != -75.527755 and lons[l] != -75.118033 and lons[l] != -75.148629 and lons[l] != -75.727202:
            text = plt.text(lons[l],lats[l],str(round(temp[l], 2)), size=9,weight='bold',verticalalignment='center',
            horizontalalignment='center',transform=ccrs.PlateCarree(),zorder=5)
            text.set_path_effects([path_effects.Stroke(linewidth=4, foreground='white'),path_effects.Normal()])
for ind in range(0,len(deos_boundarys)):
    ax.add_geometries([deos_boundarys['geometry'][ind]], ccrs.PlateCarree(),
                      facecolor='none', edgecolor='gray', zorder=3, linewidth=1.5)
for ind in range(0,len(inland_bays)):
    ax.add_geometries([inland_bays['geometry'][ind]], oldproj,
                      facecolor='white', edgecolor='black',zorder=3, linewidth=1.5)
ax.add_geometries([state_outline['geometry'][74]], oldproj, facecolor='none', edgecolor='black',zorder=3, linewidth=1.5)
ax.add_geometries([bigdeos['geometry'][121]], oldproj, facecolor='none', edgecolor='black',zorder=3, linewidth=1.5)
plt.title("(Precipitation - Ref ET) Proof of Concept")
cbar = plt.colorbar(im)
cbar.set_label('mm day-1')
plt.text(-76.13, 38.503, 'P-ET',horizontalalignment='left',color='white',weight='bold',size=9,zorder=30,transform=ccrs.PlateCarree())
plt.text(-76.13, 38.473, '2020-02-17',horizontalalignment='left',color='white',weight='bold',size=9,zorder=30,transform=ccrs.PlateCarree())

im1 = image.imread("Downloads/maplayers/deos_logo.png")
plt.figimage(im1, 350, 150 ,zorder=30, alpha=1)
plt.savefig("Downloads/P_et.png")


